import numpy as np
import pandas as pd
from active_learning_functions import train,active_learning_cycle
from active_learning_functions import generate_mutant_features,peptide2index

import sys

# LOO TCR
rand_seed = int(sys.argv[1]) #random seed

# Load epitope list
epitope_data = pd.read_csv('data/batman_train_data.csv');

# List of TCR names
TCR_list = ['18A2','868Z11','A23','A6','NYE-S1','TCR2-T','TCR3','TCR6','TCR7',
            'T1','FLT3DY']
peptide_length = 9
np.random.seed(rand_seed)

# Record mutation position
epitope_data['mutation_position']=np.sum(
    np.reshape(np.array(list(''.join(epitope_data['peptide']))) != np.array(list(
                                              ''.join(epitope_data['index']))),
    (len(epitope_data),peptide_length))*np.tile(np.arange(1,peptide_length+1),
                                                (len(epitope_data),1)),
                                                axis=1)                                                
        
#Array to save data
n_al_step = 10 # max number of AL steps
roc_auc_active = np.zeros((n_al_step,len(TCR_list))) # Array to store AUCs after AL cycles
# AL cycle-by-TCR (random choice, AL cycles...)

for left_out_TCR_index in np.arange(len(TCR_list)): # loop over TCR index for test TCR

    print(left_out_TCR_index)    
    
    # LOO TCR train data
    train_data = epitope_data[epitope_data.tcr != TCR_list[left_out_TCR_index]]
    train_data.loc[:,'tcr'] = 'pan-TCR'
    train_data.to_csv(''.join(['tmp_train_data_',str(rand_seed),'.csv']), index = False)
    
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
    inferred_weights_pan,_,inferred_matrix_pan = train(''.join(['tmp_train_data_',str(rand_seed),'.csv']),
                                                                 'full',
                                                                 'blosum100',
                                                                 steps = 80000,
                                                                 seed = 100+rand_seed)
    
    # Get mutant distance
    mutant_distance_pan = peptide2index(test_data_full['index'].tolist(), 
                              test_data_full['peptide'].tolist(),
                              inferred_matrix_pan, 
                              inferred_weights_pan.to_numpy())
    
    # Set flat weight prior mean and sd for first cycle
    flat_weights_sd = pd.DataFrame(np.zeros((1,peptide_length))) #dummy    
    prior_mean_sd = pd.concat([inferred_weights_pan,flat_weights_sd]).to_numpy()
    
    inferred_weights_pan = inferred_weights_pan.to_numpy()
    # Normalize learned weights
    inferred_weights_pan_normalized = (inferred_weights_pan-np.min(inferred_weights_pan))/(
        np.max(inferred_weights_pan)-np.min(inferred_weights_pan))        
    
    
    print(inferred_weights_pan)
    
    w_al_old = inferred_weights_pan
    #############################################################################
    # Array of mutation positions and distances     
    positional_distances = generate_mutant_features(test_data_full['index'].tolist(),
                                         test_data_full['peptide'].tolist(), 
                                         inferred_matrix_pan)
    
    mutation_position = pd.DataFrame(np.where(positional_distances!=0)[1],
                                     columns = ['position'])       
    #########################################################################       
    # Active Learning
    train_set =  test_data_index.copy() #initiate train data with the unmutated 
    # Run AL steps 
    position_to_sample = np.arange(peptide_length) #sample all positions
    for step in np.arange(n_al_step):
        # Array of mutation positions and distances     
        positional_distances = generate_mutant_features(
                                             test_data_full['index'].tolist(),
                                             test_data_full['peptide'].tolist(), 
                                             inferred_matrix_pan)
        
        mutation_position = pd.DataFrame(np.where(positional_distances!=0)[1],
                                         columns = ['position'])
        
        # Get distance based on learned weight w_AL
        distances = peptide2index(test_data_full['index'].tolist(),
                                  test_data_full['peptide'].tolist(), 
                                  inferred_matrix_pan, 
                                  w_al_old)
        
        median_distance = np.median(distances)
        _,BinEdges=pd.qcut(distances,3,retbins=True)
        
        distance1 = BinEdges[1]
        distance2 = median_distance  
        
        # rank positional weights
        w_ranks = (w_al_old).argsort().argsort()
        
        # For every position, get a mutant
        train_indices = np.zeros((peptide_length)).astype(int)    
        
        # Pick positions to select random peptides from
        random_positions = np.random.choice(peptide_length, 3, replace=False)
        for mutation in np.arange(peptide_length): # peptide positions to sample
        
            # get the closest to median distance mutant
            distances_position = distances[mutation_position.position==mutation] 
            
            # at first step, pick judiciously
            if w_ranks[0,mutation]>=5:             
                distance_index =  np.argsort(abs(distances_position-distance1))[0]
                
                train_indices[mutation] = mutation_position[
                mutation_position.position==mutation].index[distance_index]                              

            else:                
                distance_index =  np.argsort(abs(distances_position-distance2))[0]
                
                train_indices[mutation] = mutation_position[
                mutation_position.position==mutation].index[distance_index]
            
            # Pick 3 random peptides at random positions
            if mutation in random_positions:
                train_indices[mutation] = np.random.choice(mutation_position[
                mutation_position.position==mutation].index)
                   
                        
        train_indices = train_indices[position_to_sample]    
        # Train and test set
        test_set = test_data_full.iloc[np.setdiff1d(
                               np.arange(len(test_data_full)),train_indices)].copy() 
        train_set = pd.concat([train_set,test_data_full.iloc[train_indices,]])
        
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
        roc_auc_active[step,left_out_TCR_index] = auc_mean
        
        #update full test set to untested peptides
        test_data_full = test_set.copy() 
        

# Store AUC data       
roc_auc_active=pd.DataFrame(roc_auc_active,columns=TCR_list)        
roc_auc_active.index = (np.arange(n_al_step)).astype(str) 

roc_auc_active.to_csv(''.join(['data/auc_active_random_mixed_',str(rand_seed),'.csv']))

   


