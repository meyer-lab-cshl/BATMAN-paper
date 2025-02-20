import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../../batman_functions')
from batman_functions import *


#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset selected index peptides
index_peptide_list = ['FRDYVDRFYKTLRAEQASQE','FEAQKAKANKAVD','SIINFEKL',
                      'TPQDLNTML','VPSVWRSSL','FMNKFIYEI','SLLMWITQC',
                      'IMDQVPFSV','NLVPMVATV','VVVGAVGVGK']

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.index_peptide.isin(index_peptide_list)]


# Subset for TCRs with CDR3a/b sequences
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# For Ova peptide, subset for educated TCRs and OT1
tcr_pmhc_data = tcr_pmhc_data[(tcr_pmhc_data['index_peptide']!='SIINFEKL') | 
                              (tcr_pmhc_data['tcr'].str.contains('Ed')==True) |
                              (tcr_pmhc_data['tcr']=='OT1')]

# empty dataframe to store aucs
auc_loo_within_repertoire = pd.DataFrame(columns=['tcr', 'index_peptide', 'auc'])

#%%'''Run LOO-TCR predictions for BATMAN'''
for index_peptide in index_peptide_list:
    
    # Subset to antigen-specific TCR repertoire
    tcr_data = tcr_pmhc_data[tcr_pmhc_data.index_peptide==index_peptide].copy()
    
    # Rename col
    tcr_data.rename(columns={'index_peptide': 'index'}, inplace=True)
    
    # LOO TCR loop
    for tcr in pd.unique(tcr_data['tcr']):
        
        # Loo TCR data
        train_data = tcr_data[tcr_data.tcr!=tcr].copy()
        test_data = tcr_data[tcr_data.tcr==tcr].copy()
        
        # since it is pan-TCR, discard all TCR names
        train_data.tcr = 'pan-TCR'
        test_data.tcr = 'pan-TCR'
        
        # Run BATMAN to infer AA matrix and pan-TCR weights and MHC effect
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

        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        auc_tcr = np.zeros((1,3))
        
        # Loop over class pairs
        for i in np.arange(0,3):
            fpr[i], tpr[i], _ = roc_curve(test_data.activation[test_data.activation !=i],
                                          -test_data.mutant_distance[test_data.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            auc_tcr[0,i] = auc(fpr[i], tpr[i])
        
        # Record data
        new_data = pd.DataFrame({'tcr':tcr,
                                 'index_peptide':index_peptide,
                                 'auc':np.nanmean(auc_tcr)},#if 2 classes only, avoids nan
                                index=[tcr]) 

        auc_loo_within_repertoire = pd.concat([auc_loo_within_repertoire,
                                               new_data])


# Save data
auc_loo_within_repertoire.to_csv('BATMAN_within_repertoire_LOO_AUCs.csv')








