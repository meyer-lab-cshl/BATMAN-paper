import numpy as np
import pandas as pd

# Import BATMAN database
tcr_pmhc_data = pd.read_excel(''.join(['../../../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset selected TCRs
TCR_list_9mer = ['a3a',
                 #'TIL1383I', # CDR3b too long
                 'TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']

TCR_list_mhcii = ['TCR-F5','TCR-3598-2','MBP-TCR','B3K508']

TCR_list = TCR_list_9mer + TCR_list_10mer

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

# Drop TCRs without CDR3b sequence (minimal requirement for epitcr)
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# Format mhc data to conform to epitcr input
# https://github.com/ddiem-ri-4D/epiTCR/blob/main/data/hlaCovertPeudoSeq/HLAWithPseudoSeq.csv.zip

tcr_pmhc_data = tcr_pmhc_data.replace(
    ['HLA-DRB1*11:01', 'HLA-DRB1*04:01', 'I-Ab', 'HLA-A*02:01','HLA-B*81:01',
     'HLA-B*27:05', 'HLA-A*01:01', 'HLA-B*07:02', 'HLA-A*11:01'],
    ['DRB1*11:01:01:01','DRB1*04:01:01:01','H2-AB', 'A*02:01:01:01','B*81:01:01:01',
     'B*27:05:02:01','A*01:01:01:01','B*07:02:01:01','A*11:01:01:01']
    )

# Add binding category (used for AUC calculation)
tcr_pmhc_data['binder'] = tcr_pmhc_data['activation']
tcr_pmhc_data.binder[tcr_pmhc_data.binder==2] = 1

# Subset data
tcr_pmhc_data = tcr_pmhc_data[['cdr3b','peptide','mhc','binder'
                                    ]].copy().reset_index(drop=True)

tcr_pmhc_data.columns = ['CDR3b','epitope','HLA','binder']


# Get HLA seqs from https://github.com/ddiem-ri-4D/epiTCR/blob/main/data/hlaCovertPeudoSeq/HLAWithPseudoSeq.csv.zip
hla_seq = pd.read_csv('HLAWithPseudoSeq.csv',index_col=0)


# Add mhc pseudosequences
tcr_pmhc_data['MHC'] = hla_seq.seq[tcr_pmhc_data['HLA']].to_numpy()

# Write peptide seqs to input files 
tcr_pmhc_data.to_csv('epitcr_input_peptides_mhci.csv',index=None)

