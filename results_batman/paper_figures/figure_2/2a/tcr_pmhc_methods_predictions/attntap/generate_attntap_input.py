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

# Subset for TCRs with available CDR3b
tcr_pmhc_data = tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# Columns to format for Attn-TAP
tcr_pmhc_data = tcr_pmhc_data[['cdr3b','peptide','activation']]
tcr_pmhc_data.activation[tcr_pmhc_data.activation==2] = 1

tcr_pmhc_data.columns = ['tcr','antigen','label']

# Write peptide seqs to output files
tcr_pmhc_data.to_csv('attntap_input.csv',index=None)






