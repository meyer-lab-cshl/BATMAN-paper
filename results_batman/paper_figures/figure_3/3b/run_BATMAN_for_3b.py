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

# Calculate TCR-specific parameters, individually for each MHCII TCR
# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)

# Export inferred weights
inferred_weights.to_csv('inferred_weights_NLV_TCRs.csv')   
    










