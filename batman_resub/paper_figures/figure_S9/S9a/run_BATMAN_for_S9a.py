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

# Subset selected 9mer-binding TCRs with all positions sampled
tcr_list = ['a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','868Z11']

# Subset data to listed TCRs
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data['tcr'].isin(tcr_list)]

tcr_pmhc_data = tcr_pmhc_data[['tcr','index_peptide','peptide','activation',
                               'mutation_position']]


#%% '''Run within TCR predictions: selected TCRs (unpooled only)'''
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Empty df for storing alll mutant2index distances
peptide_distances = pd.DataFrame(columns=['tcr','index',
                                          'peptide','activation',
                                          'mutation_position',
                                          'd','is_al'])

# Calculate pan-TCR parameters, for LOO TCR mode
for tcr in tcr_list:
    train_data = tcr_pmhc_data[tcr_pmhc_data.tcr!=tcr].copy()
    test_data = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()
    
    train_data['tcr'] = 'pan-TCR' #in LOO TCR mode
    
    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_aa_matrix = train(train_data,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=False,
                                                steps = 80000,
                                                seed = 111)
    
    # Get 1HD mutant distances    
    distances = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              inferred_aa_matrix, 
                              inferred_weights.to_numpy())
    
    test_data['d'] = distances
    
    # Mark round 1 AL peptides
    test_data['is_al'] = 0 #initialize
    
    _,BinEdges=pd.qcut(distances,3,retbins=True)
    
    distance1 = BinEdges[1]
    distance2 = np.median(distances)  
    
    # rank positional weights
    w_ranks = (inferred_weights.to_numpy()).argsort().argsort()
    
    # For every position, get a mutant
    train_indices = np.zeros((9)).astype(int)    
    
    for mutation in np.arange(9): # peptide positions to sample
    
        #get distances for mutants with a particular mutation position
        distances_position = distances[test_data.mutation_position==(mutation+1)] 
        
        # at odd steps, pick judiciously
        if w_ranks[0,mutation]>=5:
            
            # Mutant with distance closest to distance1
            distance_index =  np.argsort(abs(distances_position-distance1))[0]
            
            train_indices[mutation] = test_data[
            test_data.mutation_position==(mutation+1)].index[distance_index]                              

        else: 
            # Mutant with distance closest to distance2
            distance_index =  np.argsort(abs(distances_position-distance2))[0]
            
            train_indices[mutation] = test_data[
            test_data.mutation_position==(mutation+1)].index[distance_index]                              
        
    test_data.loc[train_indices,'is_al'] = 1
    
    # Add to data
    peptide_distances = pd.concat([peptide_distances,test_data])
    
# Export data
peptide_distances.to_csv('peptide_distances_loo.csv')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

