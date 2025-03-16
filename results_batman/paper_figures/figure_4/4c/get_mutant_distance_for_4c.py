import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../../batman_functions')
from batman_functions import *

##### Load data ################################################
#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset selected TCRs
TCR_list = ['a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','868Z11'] 

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]
###############################################################
#Array to save data
mutant_distance = np.zeros((len(TCR_list),3)) # Array to store data
   
for tcr_index in np.arange(len(TCR_list)):
    tcr = TCR_list[tcr_index]    
    
    # Loo TCR data
    train_data = tcr_pmhc_data[tcr_pmhc_data.tcr!=tcr].copy()
    test_data = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()
    
    # since it is pan-TCR, discard all TCR names
    train_data.tcr = 'pan-TCR'
    test_data.tcr = 'pan-TCR'
    
    # Run BATMAN to infer AA matrix and pan-TCR weights
    inferred_weights,inferred_aa_matrix = train(train_data,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=False,
                                                steps = 80000,
                                                seed = 111)
    
    # Get mutant2index distance from test set as TCR activation predictor
    loo_tcr_distances = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              inferred_aa_matrix, 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              )
    # Save mean, median, and std peptide2index distance
    mutant_distance[tcr_index,0] = np.mean(loo_tcr_distances)
    mutant_distance[tcr_index,1] = np.median(loo_tcr_distances)
    mutant_distance[tcr_index,2] = np.std(loo_tcr_distances)
    
    
# Store data       
mutant_distance=pd.DataFrame(mutant_distance,
                             index=TCR_list,
                             columns=['d_mean','d_median','d_std'])        

mutant_distance.to_csv('median_mutant_distances.csv')
