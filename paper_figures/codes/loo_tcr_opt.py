import numpy as np
import pandas as pd
from active_learning_functions import train,active_learning_cycle,peptide2index
import sys

# LOO TCR
left_out_TCR_index = int(sys.argv[1]) #0 to 10

# Load epitope list
epitope_data = pd.read_csv('data/batman_train_data.csv');

# List of TCR names
TCR_list = ['18A2','868Z11','A23','A6','NYE-S1','TCR2-T','TCR3','TCR6','TCR7',
            'T1','FLT3DY']
peptide_length = 9
np.random.seed(111)

#Array to save data
n_al_step = 45 # max number of AL steps (= number of peptides)
roc_auc_active = pd.DataFrame(index=range(n_al_step),columns=["auc","peptide"]);

# LOO TCR train data
train_data = epitope_data[epitope_data.tcr != TCR_list[left_out_TCR_index]]
train_data.loc[:,'tcr'] = 'pan-TCR'
train_data.to_csv(''.join(['../data/tmp_train_data_',
                        TCR_list[left_out_TCR_index],'.csv']), index = False)

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
inferred_weights_pan,_,inferred_matrix_pan = train(''.join(['../data/tmp_train_data_',
                                             TCR_list[left_out_TCR_index],'.csv']),
                                                             'full',
                                                             'blosum100',
                                                             steps = 80000,
                                                             seed = 100)

# Get mutant2index distance for new TCR
mutant_distance_pan = peptide2index(test_data_full['index'].tolist(), 
                          test_data_full['peptide'].tolist(),
                          inferred_matrix_pan, 
                          inferred_weights_pan.to_numpy())

    
# Set flat weight prior mean and sd for first cycle
flat_weights_sd = pd.DataFrame(np.zeros((1,peptide_length))) #dummy    
prior_mean_sd = pd.concat([inferred_weights_pan,flat_weights_sd]).to_numpy()

inferred_weights_pan = inferred_weights_pan.to_numpy()

w_al_old = inferred_weights_pan

## Step 0 of OPT learning: pick median of each class #######
train_set =  test_data_index.copy() #initiate train data with the unmutated 

# In step 0 we pick 3 samples: 1 SB, 1 WB, and 1 NB, with median distances
train_indices = np.zeros((3))

for activation in np.arange(3):
    distances = mutant_distance_pan[test_data_full.activation==activation]
    indices = test_data_full[test_data_full.activation==activation].index.to_numpy()
    median_peptide = np.argsort(distances)[len(distances)//2]
    train_indices[activation] = indices[median_peptide]    

# Train and test set
test_set = test_data_full.iloc[np.setdiff1d(test_data_full.index,
                                           train_indices)].copy() 
train_set = pd.concat([train_set,test_data_full.iloc[train_indices,]])

# update prior mean
prior_mean_sd = pd.concat([pd.DataFrame(w_al_old),flat_weights_sd]).to_numpy()

# run AL
w_al,auc_mean = active_learning_cycle(
    train_set,test_set,
    inferred_matrix_pan,
    prior_mean_sd,
    steps=20000,
    seed=111)        

# Normalize learned weights
w_al = (w_al-np.min(w_al))/(np.max(w_al)-np.min(w_al))        

auc_best = auc_mean

roc_auc_active.iloc[0:4,0] = auc_best
roc_auc_active.iloc[0,1] = np.unique(train_set['index'])
roc_auc_active.iloc[1:4,1] = test_data_full.iloc[train_indices,1]

# save data    
roc_auc_active.to_csv(''.join(['../data/AL_OPT/',TCR_list[left_out_TCR_index],'_opt.csv']),
                                       index = True)

for cycle in np.arange(45):
    # Find remaining indices of unknown peptides
    remaining_indices = np.setdiff1d(np.arange(len(test_data_full)),
                                     train_indices)
        
    # Test one peptide at a time to find the best one    
    for test_peptide in remaining_indices:
        
        # variable train and test set
        train_set = test_data_full.iloc[np.append(train_indices,test_peptide),]   
        test_set = test_data_full.iloc[np.setdiff1d(test_data_full.index,
                                np.append(train_indices,test_peptide))].copy()
        
        w_al,auc_mean = active_learning_cycle(train_set,test_set,
                                              inferred_matrix_pan,
                                              prior_mean_sd,
                                              steps=20000,
                                              seed=100)
        # Check if peptide is best and update
        if auc_mean>auc_best:
            best_peptide = test_peptide
            auc_best = auc_mean
            best_weights = w_al
        print([left_out_TCR_index,cycle,test_peptide,auc_best])
        
    # Store best auc and peptide  
    roc_auc_active.iloc[cycle+4,0] = auc_best
    roc_auc_active.iloc[cycle+4,1] = test_data_full.iloc[best_peptide,1]

    # save data    
    roc_auc_active.to_csv(''.join(['../data/AL_OPT/',TCR_list[left_out_TCR_index],'_opt.csv']),
                                       index = True) 
    
    # Normalize best learned weights
    best_weights = (best_weights-np.min(best_weights))/(
                                     np.max(best_weights)-np.min(best_weights))
    # reset prior
    prior_mean_sd = pd.concat([pd.DataFrame(best_weights),
                               flat_weights_sd]).to_numpy()

    #No improvement in adding peptide
    if roc_auc_active.iloc[cycle+4,0]==roc_auc_active.iloc[cycle+3,0]: 
        break
    
    # Add best peptide to train set
    train_indices = np.append(train_indices,best_peptide)

# save data    
roc_auc_active.to_csv(''.join(['../data/AL_OPT/',TCR_list[left_out_TCR_index],'_opt.csv']),
                                       index = True)    























