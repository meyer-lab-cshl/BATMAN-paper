# Code to run BATMAN in both within-TCR (5 fold cv) and cross-TCR mode for 11 TCRs

# Import functions
import numpy as np
import pandas as pd
from original_batman_functions import train, peptide2index

# Load full TCR peptide data for 11 TCRs
full_data = pd.read_csv('fig_2_training_files/full_training_data.csv');

'''
5 fold CV in within TCR
'''
# Create dataframe to store BATMAN scores
BATMAN_scores = full_data.copy()
BATMAN_scores['within_tcr_score'] = 'NA'

for fold in np.arange(0,5):
    # Create BATMAN training data
    train_data = full_data[full_data.fold!=fold]
    train_data = train_data.rename(columns={"index_peptide":"index"})
    train_data.to_csv('fig_2_training_files/train_data.csv',
                               index = False)
    
    # Train BATMAN
    inferred_weight,inferred_matrix = train('fig_2_training_files/train_data.csv',
                                            'full','blosum100', steps = 80000,
                                            seed = 111)
    # Get BATMAN score for test data
    
    test_data = full_data[full_data.fold==fold]
    BATMAN_scores.loc[
        full_data.fold==fold,'within_tcr_score'] = peptide2index(
            test_data['index_peptide'].tolist(),
                           test_data['peptide'].tolist(),
                           inferred_matrix,
                           inferred_weight.loc[test_data.tcr].to_numpy())
    

'''
LOO TCR
'''
# Create dataframe to store BATMAN scores
BATMAN_scores['cross_tcr_score'] = 'NA'
tcr_names = pd.unique(full_data.tcr)

for TCR in tcr_names:
    # Create BATMAN training data
    train_data = full_data[full_data.tcr!=TCR]
    test_data = full_data[full_data.tcr==TCR]
    
    train_data = train_data.rename(columns={"index_peptide":"index"})
    train_data.loc[:,'tcr'] = TCR #Rename to have pan-TCR weight
    train_data.to_csv('fig_2_training_files/train_data.csv',
                               index = False)
    
    # Train BATMAN
    inferred_weight,inferred_matrix = train('fig_2_training_files/train_data.csv',
                                            'full','blosum100', steps = 80000,
                                            seed = 111)
    # Get BATMAN score for test data
    BATMAN_scores.loc[
        full_data.tcr==TCR,'cross_tcr_score'] = peptide2index(
            test_data['index_peptide'].tolist(),
                           test_data['peptide'].tolist(),
                           inferred_matrix,
                           inferred_weight.loc[test_data.tcr].to_numpy())
        
# Save BATMAN scores
BATMAN_scores.to_csv('fig_2_training_files/BATMAN_scores_full_matrix.csv',
                           index = False)
    
    
    
    
    
    
    
    
    
    
    
    
