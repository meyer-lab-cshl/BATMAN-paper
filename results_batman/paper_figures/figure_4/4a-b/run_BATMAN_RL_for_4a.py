import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

import sys

# LOO TCR
rand_seed = int(sys.argv[1]) #random seed
left_out_TCR_index = int(sys.argv[2]) #LOO TCR


np.random.seed(rand_seed)

# Join path to BATMAN functions file
from batman_AL_functions import *

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
#Array to save data
n_al_step = 10 # max number of AL steps
roc_auc_active = np.zeros((n_al_step+1,1)) # Array to store AUCs after AL cycles
# AL cycle (random choice, AL cycles...)
   
# LOO TCR train data
train_data = epitope_data[epitope_data.tcr != TCR_list[left_out_TCR_index]]
train_data.loc[:,'tcr'] = 'pan-TCR'

# LOO TCR test data
test_data_full = epitope_data[(epitope_data.tcr == TCR_list[
                        left_out_TCR_index])].copy().reset_index(drop=True)

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
                                                 seed = 111 + rand_seed)

# Get mutant distance
mutant_distance_pan = peptide2index(test_data_full['index'].tolist(), 
                          test_data_full['peptide'].tolist(),
                          inferred_matrix_pan, 
                          inferred_weights_pan.to_numpy())

  
# Calculate pairwise AUC for LOO-TCR-inferred distances
fpr = dict()
tpr = dict()
auc_al = np.zeros((1,3))
for i in np.arange(0,3):
    fpr[i], tpr[i], _ = roc_curve(test_data_full.activation[test_data_full.activation !=i],
                                  -mutant_distance_pan[test_data_full.activation !=i],
                                  pos_label=max(set(np.arange(0,3)).difference({i})))
    auc_al[0,i] = auc(fpr[i], tpr[i])        
auc_mean = np.nanmean(auc_al) 
roc_auc_active[0,0] = auc_mean

# Set flat weight prior mean and sd for first cycle
flat_weights_sd = pd.DataFrame(np.zeros((1,peptide_length))) #dummy    
prior_mean_sd = pd.concat([inferred_weights_pan,flat_weights_sd]).to_numpy()

inferred_weights_pan = inferred_weights_pan.to_numpy()
# Normalize learned weights
inferred_weights_pan_normalized = (inferred_weights_pan-np.min(inferred_weights_pan))/(
    np.max(inferred_weights_pan)-np.min(inferred_weights_pan))     

print(inferred_weights_pan)

# Store weights to be used as priors
w_al_old = inferred_weights_pan

#########################################################################      
# Active Learning cycles
train_set =  test_data_index.copy() #initiate train data with the unmutated 
# Run AL steps 
for step in np.arange(n_al_step):
    # Get distances based on learned weight w_AL
    distances = peptide2index(test_data_full['index'].tolist(),
                              test_data_full['peptide'].tolist(), 
                              inferred_matrix_pan, 
                              w_al_old)
    
    _,BinEdges=pd.qcut(distances,3,retbins=True)
    
    distance1 = BinEdges[1]
    distance2 = np.median(distances)  
    
    # rank positional weights
    w_ranks = (w_al_old).argsort().argsort()
    
    # For every position, get a mutant
    train_indices = np.zeros((peptide_length)).astype(int)    
    
    for mutation in np.arange(peptide_length): # peptide positions to sample
    
        #get distances for mutants with a particular mutation position
        distances_position = distances[test_data_full.mutation_position==(mutation+1)] 
        
        #for all steps, pick random peptides    
        train_indices[mutation] = np.random.choice(test_data_full[
        test_data_full.mutation_position==(mutation+1)].index)
               
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
        steps=40000+4000*step,
        seed=111+rand_seed)        
    
    # Normalize learned weights
    w_al = (w_al-np.min(w_al))/(np.max(w_al)-np.min(w_al))        
    
    #print(w_al)        
       
    w_al_old = w_al #store new w_al 
    roc_auc_active[step+1,0] = auc_mean
    
    #update full test set to untested peptides
    test_data_full = test_set.copy() 
        

# Store AUC data       
roc_auc_active=pd.DataFrame(roc_auc_active)        

roc_auc_active.to_csv(''.join(['BATMAN_RL_aucs/batman_RL_auc_',
                               TCR_list[left_out_TCR_index],'_',
                               str(rand_seed),'.csv']))
