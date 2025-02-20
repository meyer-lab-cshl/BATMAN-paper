import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../batman_functions')
from batman_functions import *


#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset selected TCRs
TCR_list = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 


tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

'''##############################################################'''
#%%'''Run within TCR predictions: pooled on all 9mer-binding TCRs'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_full.csv') 

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data,
                                            'symm', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_symm.csv') 

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights = train(tcr_pmhc_data,'weights_only', # No AA matrix inferred
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_BLOSUM100.csv')

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights = train(tcr_pmhc_data,'weights_only', # No AA matrix inferred
                                            'pam10', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_PAM10.csv') 

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights = train(tcr_pmhc_data,'weights_only', # No AA matrix inferred
                                            'hamming', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_HAMMING.csv') 

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights = train(tcr_pmhc_data,'weights_only', # No AA matrix inferred
                                            'dayhoff', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_DAYHOFF.csv')

'''##############################################################'''

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights = train(tcr_pmhc_data,'weights_only', # No AA matrix inferred
                                            'gonnet', #AA matrix prior
                                            consider_mhc_binding=False,
                                            steps = 80000,
                                            seed = 111)

# Export data
inferred_weights.to_csv('inferred_weights/inferred_weights_GONNET.csv') 

 











