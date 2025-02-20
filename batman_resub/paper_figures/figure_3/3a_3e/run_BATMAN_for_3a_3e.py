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

# Subset selected TCRs
TCR_list_9mer = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

#TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']

TCR_list_mhcii = ['TCR-F5','TCR-3598-2','MBP-TCR','B3K508']

TCR_list = TCR_list_9mer + TCR_list_mhcii

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

'''##############################################################'''
#%%'''Run within TCR predictions: pooled on all 9mer-binding TCRs'''
tcr_pmhc_data_9mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_9mer)].copy()

# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data_9mer,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)

# Store inferred weights
inferred_params_9mer = pd.concat([inferred_weights,inferred_mhc_effect],axis=1)

# export AA matrix
inferred_aa_matrix.to_csv('inferred_parameters/inferred_aa_matrix_9_mers.csv')

# Run BATMAN in Pan-TCR mode
tcr_pmhc_data_9mer['tcr'] = 'pan-TCR'

inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data_9mer,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)

# Store inferred weights
inferred_params_9mer_pan = pd.concat([inferred_weights,inferred_mhc_effect],axis=1)
inferred_params_9mer = pd.concat([inferred_params_9mer,inferred_params_9mer_pan])

# Export data
inferred_params_9mer.to_csv('inferred_parameters/inferred_weights_9_mers.csv')
   
'''##############################################################'''
#%% '''Run within TCR predictions: MHCII TCRs (unpooled only)'''

# Calculate TCR-specific parameters, individually for each MHCII TCR
for tcr_index in np.arange(len(TCR_list_mhcii)):
    tcr = TCR_list_mhcii[tcr_index]
    tcr_pmhc_data_mhcii = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()
    
    
    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_aa_matrix,inferred_mhc_effect = train(tcr_pmhc_data_mhcii,
                                                'full', # Asymmetric AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Export inferred weights
    inferred_params_tcr = pd.concat([inferred_weights,inferred_mhc_effect],axis=1)
    inferred_params_tcr.to_csv(''.join(['inferred_parameters/inferred_weights_',
                                        tcr,'.csv']))   
        










