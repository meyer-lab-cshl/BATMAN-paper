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

# Subset selected TCRs with NLV index peptide
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.index_peptide=='NLVPMVATV']

tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

'''##############################################################'''
#%% '''Run within TCR predictions: NLV-TCRs (all pooled)'''
pan_tcr_data = tcr_pmhc_data.copy()

pan_tcr_data['tcr'] = 'Pan-TCR'
# Calculate pan-TCR parameters
# Run BATMAN to infer AA matrix and pan-TCR weights and MHC effect
inferred_weights,inferred_aa_matrix = train(pan_tcr_data,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Read NLV-expanded TCR activity data
nlv_activity = pd.read_csv('NLV_mutant_activity.csv',index_col=0)

# Get mutant distance
nlv_activity['pan-TCR'] = peptide2index('NLVPMVATV', 
                          nlv_activity['peptide'].tolist(),
                          inferred_aa_matrix, 
                          inferred_weights.to_numpy())


#%% '''Run within TCR predictions: NLV-TCRs, individual TCRs'''

for tcr in ['NLV2','NLV3','TCR2','TCR52-10','TCR82-14']:
    tcr_data = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()
    # Calculate pan-TCR parameters
    # Run BATMAN to infer AA matrix and pan-TCR weights and MHC effect
    inferred_weights,inferred_aa_matrix = train(tcr_data,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=False,
                                                steps = 80000,
                                                seed = 111)
    
    
    # Get mutant distances
    nlv_activity[tcr] = peptide2index('NLVPMVATV', 
                              nlv_activity['peptide'].tolist(),
                              inferred_aa_matrix, 
                              inferred_weights.to_numpy())




# Export inferred weights
nlv_activity.to_csv('nlv_activity_distance.csv')   
    










