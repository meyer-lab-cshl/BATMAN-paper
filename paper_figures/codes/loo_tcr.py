import numpy as np
import pandas as pd
from active_learning_functions import train,active_learning_cycle
from active_learning_functions import generate_mutant_features,peptide2index

# Load epitope list
epitope_data = pd.read_csv('../data/batman_train_data.csv');

# List of TCR names
TCR_list = ['18A2','868Z11','A23','A6','NYE-S1','TCR2-T','TCR3','TCR6','TCR7',
            'T1','FLT3DY']
peptide_length = 9
np.random.seed(100)

for left_out_TCR_index in np.arange(len(TCR_list)): # loop over TCR index for test TCR

    print(left_out_TCR_index)    
    
    # LOO TCR train data
    train_data = epitope_data[epitope_data.tcr != TCR_list[left_out_TCR_index]]
    train_data.loc[:,'tcr'] = 'pan-TCR'
    train_data.to_csv('tmp_train_data.csv', index = False)
    
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
    inferred_weights_pan,_,inferred_matrix_pan = train('tmp_train_data.csv',
                                                                 'full',
                                                                 'blosum100',
                                                                 steps = 80000,
                                                                 seed = 100)
    
    # Save matrices and weights
    inferred_weights_pan.to_csv(''.join(['../data/loo_distances/w_',
                                         TCR_list[left_out_TCR_index],'.csv']))
    inferred_matrix_pan.to_csv(''.join(['../data/loo_distances/aa_',
                                         TCR_list[left_out_TCR_index],'.csv']))
    
    
    