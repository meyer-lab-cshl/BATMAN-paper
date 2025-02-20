import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../batman_functions')
from batman_functions import *


#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset to 9-mer and 10-mer binding selected TCRs only
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data['peptide'].str.len().isin([9,10])].copy()

# Extract 9-mer and 10-mer TCRs
TCR_list_9mer = pd.unique(tcr_pmhc_data.tcr[tcr_pmhc_data['peptide'].str.len()==9]) 

TCR_list_10mer = pd.unique(tcr_pmhc_data.tcr[tcr_pmhc_data['peptide'].str.len()==10]) 

TCR_list = list(TCR_list_9mer) + list(TCR_list_10mer)

#Array to store AUCs: TCR-class_pair-fold
aucs_tcr = np.zeros((len(TCR_list),3,5)) 

#%%'''Run within TCR predictions: unpooled on TCRs'''
for tcr_index in np.arange(len(TCR_list)):
    tcr = TCR_list[tcr_index]
    
    # Subset to TCR-specific data
    tcr_data = tcr_pmhc_data[tcr_pmhc_data.tcr == tcr].copy()

    # Loop over 5 folds
    for fold in np.arange(5):
        train_data = tcr_data[tcr_data.training_fold!=fold].copy()
        test_data = tcr_data[tcr_data.training_fold==fold].copy()
        
        # Add index peptide in training data for all fold
        train_data = pd.concat([train_data,
                                tcr_data[pd.isna(tcr_data.training_fold)]])
    
        # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
        inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(train_data,
                                                    'full', # Asymmetric AA matrix
                                                    'blosum100', #AA matrix prior
                                                    consider_mhc_binding=True,
                                                    steps = 80000,
                                                    seed = 111)
        
        # Get test set mutant2index distance and MHC effect as TCR activation predictor
        test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                                  test_data['peptide'].tolist(),
                                  inferred_aa_matrix, 
                                  inferred_weights.loc[test_data['tcr'],:].to_numpy()
                                  ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                              test_data['BindLevel']
                                                              ].to_numpy())        
        
        
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        
        # Loop over class pairs
        for i in np.arange(0,3):
            fpr, tpr, _ = roc_curve(test_data.activation[test_data.activation !=i],
                                          -test_data.mutant_distance[test_data.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            aucs_tcr[tcr_index,i,fold] = auc(fpr, tpr)
        
            
# Average over folds and save
aucs_all = np.nanmean(aucs_tcr,axis=2) #nanmean helps if a test fold has only 2 classes   

'''##############################################################'''
# All AUCs: Within-TCR
aucs_all=pd.DataFrame(aucs_all,
                      index=TCR_list,columns=['WS','NS','NW'])
aucs_all['type'] = 'Full_unpooled'

aucs_all.to_csv('aucs_Full_unpooled.csv')











