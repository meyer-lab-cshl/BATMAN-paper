import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../../batman_functions')
from batman_functions import *

# Inputs
tcr_index = int(sys.argv[1]) #Selected TCR

#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Full list of all TCRs
TCR_list = pd.unique(tcr_pmhc_data.tcr)

# Subset data for selected TCR
tcr = TCR_list[tcr_index]
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr]

# AA matrix list
aa_matrix_list = ['hamming','blosum62'] + list(np.char.add(np.repeat('blosum',15),
                                  np.arange(30,105,5).astype(str)))

aa_matrix_list = aa_matrix_list + list(np.char.add(np.repeat('pam',50),
                                  np.arange(10,510,10).astype(str)))

aa_matrix_list = aa_matrix_list + ['dayhoff','gonnet']

aa_matrix_list.remove('blosum95') #absent aa matrix


#Array to store correlation results of different binding pairs and aa matrices
auc_all = np.zeros((len(aa_matrix_list),3)) 

for aa_matrix_index in np.arange(len(aa_matrix_list)):
    
    # select AA matrix
    aa_matrix = aa_matrix_list[aa_matrix_index]
    
    auc_fold = np.zeros((5,3)) #5 Fold, 3 activation classes
    # Loop over 5 folds
    for fold in np.arange(5):
        train_data = tcr_pmhc_data[tcr_pmhc_data.training_fold!=fold].copy()
        test_data = tcr_pmhc_data[tcr_pmhc_data.training_fold==fold].copy()    
        
        # Run BATMAN to infer weights
        inferred_weights = train(train_data,
                                 'weights_only', # Asymmetric AA matrix
                                    aa_matrix, #AA matrix prior
                                    consider_mhc_binding=False,
                                    steps = 50000,
                                    seed = 111)
    
        
        # Get mutant2index distance from test set as TCR activation predictor
        test_distances = peptide2index(test_data['index'].tolist(), 
                                  test_data['peptide'].tolist(),
                                  aa_matrix, 
                                  inferred_weights.to_numpy())
        
        # Measure AUCs
        # Loop over class pairs
        for i in np.arange(0,3):
            fpr, tpr, _ = roc_curve(test_data.activation[test_data.activation !=i],
                                          -test_distances[test_data.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            auc_fold[fold,i] = auc(fpr, tpr)
            
    #Save data
    auc_all[aa_matrix_index,:] = np.nanmean(auc_fold,axis=0)

# Save data
auc_all = pd.DataFrame(auc_all,columns=['WS','NS','NW'])
auc_all['aa_matrix'] = aa_matrix_list
auc_all['tcr'] = tcr

auc_all.to_csv(''.join(['aucs_TCR_AA_matrix/auc_',tcr,'.csv']))








