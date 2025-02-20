import numpy as np
import pandas as pd

# Import BATMAN database
tcr_pmhc_data = pd.read_excel(''.join(['../../../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset selected TCRs
TCR_list_9mer = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']

TCR_list_mhcii = ['TCR-F5','TCR-3598-2','MBP-TCR','B3K508']

TCR_list = TCR_list_9mer + TCR_list_10mer + TCR_list_mhcii

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

# Subset to TCRs with seq
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# Format data according to TCRPrediction format
# https://github.com/kyoheikoyama/TCRPrediction/tree/main
tcr_pmhc_data = tcr_pmhc_data[['peptide','cdr3a','cdr3b']]

tcr_pmhc_data.columns = ['peptide','tcra','tcrb']
tcr_pmhc_data['sign'] = 0 #dummy

# Export data
tcr_pmhc_data.to_csv('tcrprediction_input_peptides.csv',index=None)