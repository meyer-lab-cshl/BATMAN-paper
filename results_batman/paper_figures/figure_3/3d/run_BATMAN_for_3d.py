import numpy as np
import pandas as pd

# Join path to BATMAN functions file
from batman_functions_for_posterior_analysis import *


#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset SIINFEKL TCRs
tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.index_peptide=='SIINFEKL']

tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)


   
'''##############################################################'''
#%% '''Run within TCR predictions: selected TCRs (pooled only)'''
# Calculate TCR-specific parameters
# Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
inferred_weights,inferred_weights_sd,_,_ = train(tcr_pmhc_data,
                                            'full', # Asymmetric AA matrix
                                            'blosum100', #AA matrix prior
                                            consider_mhc_binding=True,
                                            steps = 80000,
                                            seed = 111)


# Export inferred weights and sd of weights
inferred_weights.to_csv('inferred_weights_SIINFEKL.csv')   
inferred_weights_sd.to_csv('inferred_weights_sd_SIINFEKL.csv')   

    






