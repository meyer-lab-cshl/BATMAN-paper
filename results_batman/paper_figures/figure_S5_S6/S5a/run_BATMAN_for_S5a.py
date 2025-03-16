import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os
import scipy.stats

from batman_regression_functions import *

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


#Array to store correlation results of different folds and aa matrices
corr_fold = np.zeros((len(aa_matrix_list),5)) 

for aa_matrix_index in np.arange(len(aa_matrix_list)):
    
    # select AA matrix
    aa_matrix = aa_matrix_list[aa_matrix_index]
    # Loop over 5 folds
    for fold in np.arange(5):
        train_data = tcr_pmhc_data[tcr_pmhc_data.training_fold!=fold].copy()
        test_data = tcr_pmhc_data[tcr_pmhc_data.training_fold==fold].copy()    
        
        # Run BATMAN to infer weights
        inferred_weights = train(tcr_pmhc_data,#Path to file or pandas df with TCR data (see example for format)
                  # Default values of optional parameters
                  aa_matrix=aa_matrix,#Named or user-defined AA matrix used for regularization
                  seed=100,#seed for sampling
                  steps=50000)# number of steps for sampling
    
        
        # Get mutant2index distance from test set as TCR activation predictor
        test_distances = peptide2index(test_data['index'].tolist(), 
                                  test_data['peptide'].tolist(),
                                  aa_matrix, 
                                  inferred_weights)
        
        # Measure correlation
        corr_fold[aa_matrix_index,fold] = np.abs(scipy.stats.spearmanr(test_data['peptide_activity'],
                                                  test_distances).correlation)

# Save data
corr_mean = np.mean(corr_fold,axis=1)

corr_av = pd.DataFrame(corr_mean,columns=['Spearman r'])
corr_av['aa_matrix'] = aa_matrix_list
corr_av['tcr'] = tcr

corr_av.to_csv(''.join(['Corrs_TCR_AA_matrix/corr_',tcr,'.csv']))








