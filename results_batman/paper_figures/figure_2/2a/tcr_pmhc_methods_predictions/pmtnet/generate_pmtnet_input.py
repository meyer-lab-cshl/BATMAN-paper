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

# Drop TCRs without CDR3b sequence (minimal requirement for pmtnet)
tcr_pmhc_data=tcr_pmhc_data[pd.isna(tcr_pmhc_data.cdr3b)==False]

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
                                         
# Additional TRAV/TRBV changes to conform input
# Available alleles: https://github.com/Yuqiu-Yang/pMTnet_Omni_Document/tree/main/pMTnet_Omni_Document/validation_data
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['TRAV24-1*01','TRAV25*01 F','TRAV5*00(1661.5)','TRAV2-3*01'],
    ['TRAV24*01','TRAV25*01','TRAV5*01','TRAV2*01']
    )

tcr_pmhc_data = tcr_pmhc_data.replace(
    ['TRBV2-1*01','TRBV27*01 F','TRBV20-1*00(1904.2)'],
    ['TRBV2*01','TRBV27*01','TRBV20-1*01']
    )

# End CDR3 seqs in F if F is there
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['CAASKGNNRIFFG','CAWSLGGVAETLYFG'],
    ['CAASKGNNRIFF','CAWSLGGVAETLYF']
    )

# Mouse MHC format
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['I-Ab'],
    ['H-2-IAb']
    ) 

# Species name formats
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['Human'],
    ['human']
    ) 

# Rename columns: 
# https://pmtnet-omni-document.readthedocs.io/en/latest/quick_start.html
tcr_pmhc_data['pMHC_SPECIES'] = tcr_pmhc_data['tcr_source_organism'] 
tcr_pmhc_data['TCR_SPECIES'] = tcr_pmhc_data['tcr_source_organism'] 
tcr_pmhc_data['va'] = tcr_pmhc_data['trav']
tcr_pmhc_data['vb'] = tcr_pmhc_data['trbv']

# Subset data
tcr_pmhc_data = tcr_pmhc_data[['va','cdr3a','vb','cdr3b','peptide','mhc',
                               'pMHC_SPECIES','TCR_SPECIES']]

# Write peptide seqs to input files in 7 batches (max 460 TCR-pMHC pairs)
df_batches = np.array_split(tcr_pmhc_data, 7)

for batch in np.arange(7):
    df = df_batches[batch].reset_index(drop=True)
    df.to_csv(''.join(['pmtnet_input_peptides_',str(batch),'.csv']))

# Separate data without trav/trbv info for pMTnet v1 
tcr_pmhc_data = tcr_pmhc_data[pd.isna(tcr_pmhc_data.va)==True].reset_index(drop=True)
tcr_pmhc_data = tcr_pmhc_data[['cdr3b','peptide','mhc']]
tcr_pmhc_data.columns = ['CDR3','Antigen','HLA']
# MHC format
tcr_pmhc_data = tcr_pmhc_data.replace(
    ['HLA-A*02:01'],
    ['A*02:01']
    )


tcr_pmhc_data.to_csv('pmtnet_input_peptides_7.csv')

