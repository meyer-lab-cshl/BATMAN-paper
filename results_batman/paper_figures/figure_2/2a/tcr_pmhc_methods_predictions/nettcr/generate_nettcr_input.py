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

# Drop TCRs without TRAV/ sequence (minimal requirement for nettcr)
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.trav)==False]

# Format data to conform to pmtnet input
# https://pmtnet-omni-document.readthedocs.io/en/latest/input_format/index.html

# Add TRAV, TRBV before gene names if available
tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trav)==False,
              'trav'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trav)==False,
                            'trav'].astype(str).apply(lambda x: 'TRAV' + x)

tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbv)==False,
              'trbv'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbv)==False,
                            'trbv'].astype(str).apply(lambda x: 'TRBV' + x)                                          

# Replace . by -
tcr_pmhc_data.trav[tcr_pmhc_data.trav=='TRAV2.3'] = 'TRAV2-3'

                                         
# Add allele *01 where allele names are absent
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trav)==False) & 
                 (tcr_pmhc_data.trav.str.contains('\*')==False),
             'trav'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trav)==False) & 
                            (tcr_pmhc_data.trav.str.contains('\*')==False),
                                'trav'].astype(str).apply(lambda x: x + '*01')
                                         
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbv)==False) & 
                 (tcr_pmhc_data.trbv.str.contains('\*')==False),
             'trbv'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbv)==False) & 
                            (tcr_pmhc_data.trbv.str.contains('\*')==False),
                                'trbv'].astype(str).apply(lambda x: x + '*01')
                                         
# Take out C from cdr3ab sequences
tcr_pmhc_data.loc[:,'cdr3a'] = tcr_pmhc_data.loc[:,'cdr3a'].astype(str).apply(
                                    lambda x: x.replace('C', ''))   

tcr_pmhc_data.loc[:,'cdr3b'] = tcr_pmhc_data.loc[:,'cdr3b'].astype(str).apply(
                                    lambda x: x.replace('C', ''))

# Add CDR12a and CDR12b sequences
cdr12a = pd.read_excel('IMGT_CDR12_in_TRAV.xlsx',index_col=0)
cdr12b = pd.read_excel('IMGT_CDR12_in_TRBV.xlsx',index_col=0)

tcr_pmhc_data['cdr1a'] = cdr12a.cdr1a[tcr_pmhc_data['trav']].to_numpy()
tcr_pmhc_data['cdr2a'] = cdr12a.cdr2a[tcr_pmhc_data['trav']].to_numpy()

tcr_pmhc_data['cdr1b'] = cdr12b.cdr1b[tcr_pmhc_data['trbv']].to_numpy()
tcr_pmhc_data['cdr2b'] = cdr12b.cdr2b[tcr_pmhc_data['trbv']].to_numpy()

# Subset data
tcr_pmhc_data = tcr_pmhc_data[['peptide','cdr1a','cdr2a','cdr3a',
                               'cdr1b','cdr2b','cdr3b']].reset_index(drop=True)
# Save data
tcr_pmhc_data.to_csv('nettcr_input_peptides.csv')

