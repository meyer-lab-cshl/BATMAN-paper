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

#%%'''Run within TCR predictions: pooled on 9mer-binding TCRs'''
# Pool within 9-mer-binding TCRs
tcr_pmhc_data_9mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_9mer)].copy()
aucs_9mer = np.zeros((len(TCR_list_9mer),3)) 

#Array to store AUCs: TCR-class_pair-fold
aucs_tcr = np.zeros((len(TCR_list_9mer),3,5)) 
# Loop over 5 folds
for fold in np.arange(5):
    train_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.training_fold!=fold].copy()
    test_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.training_fold==fold].copy()
    
    # Add index peptide in training data for all fold
    train_data = pd.concat([train_data,
                            tcr_pmhc_data_9mer[pd.isna(tcr_pmhc_data_9mer.training_fold)]])

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
    
    # Calculate TCR-specific AUC, averaged over 3 class-pairs
    for tcr_index in np.arange(len(TCR_list_9mer)):
        tcr = TCR_list_9mer[tcr_index]
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        
        # Loop over class pairs
        for i in np.arange(0,3):
            test_data_tcr = test_data[test_data.tcr==tcr].copy()
            fpr[i], tpr[i], _ = roc_curve(test_data_tcr.activation[test_data_tcr.activation !=i],
                                          -test_data_tcr.mutant_distance[test_data_tcr.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            aucs_tcr[tcr_index,i,fold] = auc(fpr[i], tpr[i])
        
        
# Average over folds and save
aucs_9mer = np.nanmean(aucs_tcr,axis=2) #nanmean helps if a test fold has only 2 classes   

'''###############################################################'''

#%%'''Run within TCR predictions: pooled on 10-mer-binding TCRs'''
# Pool within 10-mer-binding TCRs
tcr_pmhc_data_10mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_10mer)].copy()

#Array to store AUCs: TCR-class_pair-fold
aucs_tcr = np.zeros((len(TCR_list_10mer),3,5)) 

# Loop over 5 folds
for fold in np.arange(5):
    train_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.training_fold!=fold].copy()
    test_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.training_fold==fold].copy()
    
    # Add index peptide in training data for all fold
    train_data = pd.concat([train_data,
                            tcr_pmhc_data_10mer[pd.isna(tcr_pmhc_data_10mer.training_fold)]])

    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(train_data,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Get mutant2index distance from test set as TCR activation predictor
    test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              inferred_aa_matrix, 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                          test_data['BindLevel']
                                                          ].to_numpy()) 
    
    # Calculate TCR-specific AUC, averaged over 3 class-pairs
    for tcr_index in np.arange(len(TCR_list_10mer)):
        tcr = TCR_list_10mer[tcr_index]
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        
        # Loop over class pairs
        for i in np.arange(0,3):
            test_data_tcr = test_data[test_data.tcr==tcr].copy()
            fpr[i], tpr[i], _ = roc_curve(test_data_tcr.activation[test_data_tcr.activation !=i],
                                          -test_data_tcr.mutant_distance[test_data_tcr.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            aucs_tcr[tcr_index,i,fold] = auc(fpr[i], tpr[i])
        
        

# Average over folds and save
aucs_10mer = np.nanmean(aucs_tcr,axis=2) #nanmean helps if a test fold has only 2 classes   
    

'''##############################################################'''
# All AUCs: Within-TCR
aucs_all=pd.DataFrame(np.concatenate((aucs_9mer,aucs_10mer),axis=0),
                      index=TCR_list,columns=['WS','NS','NW'])
aucs_all['type'] = 'Full_across'

aucs_all.to_csv('aucs_Full_across.csv')











