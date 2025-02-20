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
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset selected TCRs with PDB structures
TCR_list = pd.read_excel('PDB_ID_list.xlsx')['tcr']
   
'''##############################################################'''
#%% '''Run within TCR predictions: selected TCRs (unpooled only)'''

# Calculate TCR-specific parameters, individually for each TCR
for tcr in TCR_list:
    tcr_pmhc_data_subset = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()    
    
    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data_subset,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Export inferred weights
    inferred_weights.to_csv(''.join(['inferred_weights/inferred_weights_',
                                        tcr,'.csv']))   
        










