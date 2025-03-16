import numpy as np
import pandas as pd

# Import BATMAN database
tcr_pmhc_data = pd.read_excel(''.join(['../../../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset selected TCRs
TCR_list_9mer = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']


TCR_list = TCR_list_9mer + TCR_list_10mer

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

# Available MHCs for NetTepi: 
# https://services.healthtech.dtu.dk/services/NetTepi-1.0/MHC_allele_names.txt
nettepi_mhc_list = ['HLA-A*02:01', 'HLA-B*27:05', 'HLA-A*01:01',
       'HLA-B*07:02', 'HLA-A*24:02', 'HLA-A*11:01']

# Write peptide seqs to output files separately for 9mer and 10mer
for mhc in nettepi_mhc_list:
    #9mers
    peptides = tcr_pmhc_data.peptide[(tcr_pmhc_data.mhc==mhc) & 
                                     (tcr_pmhc_data.tcr.isin(TCR_list_9mer))]
    if len(peptides)>0: #non-empty
        peptides.to_csv(''.join(['nettepi_input_9mer_peptides_',
                             mhc.replace(':', '').replace('*', ''),
                             '.csv']),index=None)
    #10mers    
    peptides = tcr_pmhc_data.peptide[(tcr_pmhc_data.mhc==mhc) & 
                                     (tcr_pmhc_data.tcr.isin(TCR_list_10mer))]
    if len(peptides)>0: #non-empty
        peptides.to_csv(''.join(['nettepi_input_10mer_peptides_',
                             mhc.replace(':', '').replace('*', ''),
                             '.csv']),index=None)

