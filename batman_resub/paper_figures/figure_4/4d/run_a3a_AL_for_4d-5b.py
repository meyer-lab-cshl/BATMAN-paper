import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN AL functions file
sys.path.append('../4a-b')
from batman_AL_functions import *


'''Active Learning'''

##### Load data ################################################
epitope_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
epitope_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset selected 9mer-binding TCRs with all positions sampled
TCR_list = ['a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','868Z11'] 

epitope_data = epitope_data[epitope_data.tcr.isin(TCR_list)]

peptide_length = 9
###############################################################
# Run LOO TCR and AL for a3a and record weights  
# LOO TCR train data
train_data = epitope_data[epitope_data.tcr != 'a3a']
train_data.loc[:,'tcr'] = 'pan-TCR'

# LOO TCR test data
test_data_full = epitope_data[(epitope_data.tcr == 'a3a')].copy().reset_index(drop=True)

# identify unmutated peptide and take it out
test_data_index = test_data_full[
 test_data_full.peptide==test_data_full['index'][0]].copy().reset_index(drop=True)

test_data_full = test_data_full[
 test_data_full.peptide!=test_data_full['index'][0]].reset_index(drop=True)

################# Pan TCR #################################################
# Train BATMAN to infer pan-TCR weights and AA distance matrix
inferred_weights_pan,inferred_matrix_pan = train(train_data,'full',
                                                 'blosum100',
                                                 steps = 80000,
                                                 seed = 111)

# Get mutant distance for 1 Hamming
mutant_distance_pan = peptide2index(test_data_full['index'].tolist(), 
                          test_data_full['peptide'].tolist(),
                          inferred_matrix_pan, 
                          inferred_weights_pan.to_numpy())

# Set flat weight prior mean and sd for first cycle
flat_weights_sd = pd.DataFrame(np.zeros((1,peptide_length))) #dummy    
prior_mean_sd = pd.concat([inferred_weights_pan,flat_weights_sd]).to_numpy()

inferred_weights_pan = inferred_weights_pan.to_numpy()

# Store weights to be used as priors
w_al_old = inferred_weights_pan

#########################################################################      
# Active Learning cycle: choose peptides
train_set =  test_data_index.copy() #initiate train data with the unmutated 
    
_,BinEdges=pd.qcut(mutant_distance_pan ,3,retbins=True)

distance1 = BinEdges[1]
distance2 = np.median(mutant_distance_pan )  

# rank positional weights
w_ranks = (w_al_old).argsort().argsort()

# For every position, get a mutant
train_indices = np.zeros((peptide_length)).astype(int)    

for mutation in np.arange(peptide_length): # peptide positions to sample

    #get distances for mutants with a particular mutation position
    distances_position = mutant_distance_pan[test_data_full.mutation_position==(mutation+1)] 
    
    # Pick AL peptides
    if w_ranks[0,mutation]>=5:
        
        # Mutant with distance closest to distance1
        distance_index =  np.argsort(abs(distances_position-distance1))[0]
        
        train_indices[mutation] = test_data_full[
        test_data_full.mutation_position==(mutation+1)].index[distance_index]                              

    else: 
        # Mutant with distance closest to distance2
        distance_index =  np.argsort(abs(distances_position-distance2))[0]
        
        train_indices[mutation] = test_data_full[
        test_data_full.mutation_position==(mutation+1)].index[distance_index]                              
        
# Train and test set
test_set = test_data_full.loc[np.setdiff1d(test_data_full.index,
                                           train_indices)].copy() 
train_set = pd.concat([train_set,test_data_full.loc[train_indices,]])

# update prior mean
prior_mean_sd = pd.concat([pd.DataFrame(w_al_old),flat_weights_sd]).to_numpy()

# run AL
w_al,auc_mean = active_learning_cycle(
    train_set,test_set,
    inferred_matrix_pan,
    prior_mean_sd,
    steps=40000,
    seed=111)     
    
# Save weight profiles
weights = pd.DataFrame(np.concatenate([inferred_weights_pan,w_al],axis=0),
                       index=['Pan-TCR','AL'])

weights.to_csv('a3a_weights.csv')



'''Calculate distances of 1HD and >1 HD mutants with Pan-TCR and AL weights'''
# Read >1 Hamming data
multi_hamming_a3a = pd.read_csv(''.join(['../../../tcr_epitope_datasets/',
                                         'more_than_one_hamming_datasets/MAGE_papers/',
                                          'a3a_multi_hamming_all.csv']))

# Filter 9mers
multi_hamming_a3a = multi_hamming_a3a[
    np.char.str_len(multi_hamming_a3a.peptide.to_numpy().astype(str))==9]

multi_hamming_a3a.rename(columns={'index_peptide': 'index'}, inplace=True)

#Add to 1 Hamming data
single_hamming_a3a = test_set[['tcr','index','activation','peptide']]

all_mutants_a3a = pd.concat([single_hamming_a3a,multi_hamming_a3a])

# Find Hamming, pan-TCR, and AL mutant distances
all_mutants_a3a['d_hamming'] = peptide2index(all_mutants_a3a['index'].tolist(), 
                          all_mutants_a3a['peptide'].tolist(),
                          'hamming', 
                          np.ones((1,9)))

all_mutants_a3a['d_pan'] = peptide2index(all_mutants_a3a['index'].tolist(), 
                          all_mutants_a3a['peptide'].tolist(),
                          inferred_matrix_pan, 
                          inferred_weights_pan)

all_mutants_a3a['d_al'] = peptide2index(all_mutants_a3a['index'].tolist(), 
                          all_mutants_a3a['peptide'].tolist(),
                          inferred_matrix_pan, 
                          w_al)

# Save data
all_mutants_a3a.to_csv('mutant_distances_a3a.csv')














