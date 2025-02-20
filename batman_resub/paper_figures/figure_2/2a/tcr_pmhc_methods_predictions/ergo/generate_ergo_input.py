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

# Drop TCRs without CDR3b sequence (minimal requirement for ergo)
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

# Format data to conform to ergo input
# https://github.com/IdoSpringer/ERGO-II/blob/master/example.csv

# Add TRAV, TRBV etc before gene names if available
tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trav)==False,
              'trav'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trav)==False,
                            'trav'].astype(str).apply(lambda x: 'TRAV' + x)

tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbv)==False,
              'trbv'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbv)==False,
                            'trbv'].astype(str).apply(lambda x: 'TRBV' + x)

tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.traj)==False,
              'traj'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.traj)==False,
                            'traj'].astype(str).apply(lambda x: 'TRAJ' + x)

tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbj)==False,
              'trbj'] = tcr_pmhc_data.loc[pd.isna(tcr_pmhc_data.trbj)==False,
                            'trbj'].astype(str).apply(lambda x: 'TRBJ' + x)                                          

# Replace . by -
tcr_pmhc_data.trav[tcr_pmhc_data.trav=='TRAV2.3'] = 'TRAV2-3'
tcr_pmhc_data.trbj[tcr_pmhc_data.trbj=='TRBJ2.3'] = 'TRBJ2-3'


# delete allele *01 where allele names are present
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trav)==False) & 
                 (tcr_pmhc_data.trav.str.contains('\*')==True),
             'trav'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trav)==False) & 
                            (tcr_pmhc_data.trav.str.contains('\*')==True),
                                'trav'].astype(str).apply(lambda x: x.split('*')[0])
                                         
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbv)==False) & 
                 (tcr_pmhc_data.trbv.str.contains('\*')==True),
             'trbv'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbv)==False) & 
                            (tcr_pmhc_data.trbv.str.contains('\*')==True),
                                'trbv'].astype(str).apply(lambda x: x.split('*')[0])
                                         
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.traj)==False) & 
                 (tcr_pmhc_data.traj.str.contains('\*')==True),
             'traj'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.traj)==False) & 
                            (tcr_pmhc_data.traj.str.contains('\*')==True),
                                'traj'].astype(str).apply(lambda x: x.split('*')[0])
                                         
tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbj)==False) & 
                 (tcr_pmhc_data.trbj.str.contains('\*')==True),
             'trbj'] = tcr_pmhc_data.loc[(pd.isna(tcr_pmhc_data.trbj)==False) & 
                            (tcr_pmhc_data.trbj.str.contains('\*')==True),
                                'trbj'].astype(str).apply(lambda x: x.split('*')[0])                                         
                                         
# Additional TRAV/TRBV changes to conform input
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['TRAJ24 (31)', 'TRBV 5-1', 'TRBJ 2-7'],
    ['TRAJ24', 'TRBV5-1', 'TRBJ2-7']
    )

# MHC format
tcr_pmhc_data.loc[:,'mhc'] = tcr_pmhc_data.loc[:,
                                'mhc'].astype(str).apply(lambda x: x.split(':')[0])                                         
   
# T cell type column: 
tcr_pmhc_data['T-Cell-Type'] = 'CD8'
tcr_pmhc_data.loc[tcr_pmhc_data.tcr.isin(['TCR-F5','TCR-3598-2','B3K508']),
              'T-Cell-Type'] = 'CD4'
# Subset data
tcr_pmhc_data = tcr_pmhc_data[['cdr3a','cdr3b','trav','traj','trbv','trbj',
                               'T-Cell-Type','peptide','mhc']]
tcr_pmhc_data.columns = ['TRA','TRB','TRAV','TRAJ','TRBV','TRBJ','T-Cell-Type',
                         'Peptide','MHC']
# Export data
tcr_pmhc_data.to_csv('ergo_input_peptides.csv')

