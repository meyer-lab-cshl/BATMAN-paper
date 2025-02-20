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

# List of TCRs with matched 1HD and >1 HD data
tcr_list = ['868Z11','a3a','EWW-TCR','TIL1383I','c259','TCR-1E6',
            'MBP-TCR','A3-05', 'A3-10']

# Subset data to listed TCRs
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data['tcr'].isin(tcr_list)]

# Load >1 HD data for all TCRs
tcr_pmhc_data_multihamming = pd.DataFrame(columns=['tcr', 'peptide',
                                                   'index_peptide','activation'])

for tcr in ['868Z11','EWW-TCR','TIL1383I','c259','TCR-1E6','MBP-TCR']:
    
    # Open file
    tcr_data = pd.read_csv(''.join(['../../../tcr_epitope_datasets/',
              'more_than_one_hamming_datasets/',tcr,'_papers/',tcr,
              '_multi_hamming_all.csv']))
    
    # Add data
    tcr_pmhc_data_multihamming = pd.concat([tcr_pmhc_data_multihamming,
                                            tcr_data])

# MAGE TCRs
tcr_data = pd.read_csv(''.join(['../../../tcr_epitope_datasets/',
          'more_than_one_hamming_datasets/MAGE_papers/',
          'a3a_multi_hamming_all.csv']))

# Add a3a data
tcr_pmhc_data_multihamming = pd.concat([tcr_pmhc_data_multihamming,
                                        tcr_data])

# Read data for A3-05 and A3-10 TCRs, plus more a3a data
tcr_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'more_than_one_hamming_datasets/MAGE_papers/',
          'MAGE_self_peptides_Immunity_paper.xlsx']))

tcr_data_a3a = tcr_data.iloc[:,[0,1]].copy()
tcr_data_a3a.columns = ['peptide','activation']
tcr_data_a3a['tcr'] = 'a3a'

tcr_data_A305 = tcr_data.iloc[:,[0,2]].copy()
tcr_data_A305.columns = ['peptide','activation']
tcr_data_A305['tcr'] = 'A3-05'

tcr_data_A310 = tcr_data.iloc[:,[0,3]].copy()
tcr_data_A310.columns = ['peptide','activation']
tcr_data_A310['tcr'] = 'A3-10'

tcr_data = pd.concat([tcr_data_a3a,tcr_data_A305,tcr_data_A310])
tcr_data['index_peptide'] = 'EVDPIGHLY'


# Add data
tcr_pmhc_data_multihamming = pd.concat([tcr_pmhc_data_multihamming,
                                        tcr_data])

# Add 1HD data
tcr_pmhc_data_multihamming = pd.concat([tcr_pmhc_data_multihamming,
                                        tcr_pmhc_data[['tcr','peptide',
                                                       'index_peptide',
                                                       'activation']]])

# Remove duplicates
tcr_pmhc_data_multihamming = tcr_pmhc_data_multihamming.drop_duplicates().reset_index(drop=True)

# Remove instances when index peptide and peptide lengths do not match
tcr_pmhc_data_multihamming = tcr_pmhc_data_multihamming[
    np.char.str_len(tcr_pmhc_data_multihamming.peptide.to_numpy().astype(str)
                    )==np.char.str_len(
                        tcr_pmhc_data_multihamming.index_peptide.to_numpy().astype(str)
                                    )]                     
                        
#%% '''Run within TCR predictions: selected TCRs (unpooled only)'''
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Empty df for storing alll mutant2index distances
peptide_distances = pd.DataFrame(columns=['tcr', 'peptide',
                                                   'index_peptide','activation',
                                                   'd_hamming','d'])

# Calculate TCR-specific parameters, individually for each TCR
for tcr in tcr_list:
    tcr_pmhc_data_subset = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()    
    
    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_aa_matrix = train(tcr_pmhc_data_subset,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=False,
                                                steps = 80000,
                                                seed = 111)
    
    # Get 1HD and >1 HD mutant distances
    
    # Subset TCR specific data
    tcr_data = tcr_pmhc_data_multihamming[tcr_pmhc_data_multihamming.tcr==tcr].copy()

    # Find Hamming and inferred mutant distances
    peptide_length = np.unique(np.char.str_len(
                                      tcr_data.peptide.to_numpy().astype(str)))[0]
    
    tcr_data['d_hamming'] = peptide2index(tcr_data['index_peptide'].tolist(), 
                              tcr_data['peptide'].tolist(),
                              'hamming', 
                              np.ones((1,peptide_length)))

    tcr_data['d'] = peptide2index(tcr_data['index_peptide'].tolist(), 
                              tcr_data['peptide'].tolist(),
                              inferred_aa_matrix, 
                              inferred_weights.to_numpy())
    
    # Add to data
    peptide_distances = pd.concat([peptide_distances,tcr_data])
    
# Export data
#peptide_distances.to_csv('peptide_distances_multihamming.csv')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

