import numpy as np
import pandas as pd

# Read CMV data
cmv_data = pd.read_excel("Affinity_prediction_cmv.xlsx",
                            sheet_name = "Mean",index_col=0)
tcr_info = pd.read_excel("Affinity_prediction_cmv.xlsx",
                            sheet_name = "TCR sequences")
peptide_info = pd.read_excel("Affinity_prediction_cmv.xlsx",
                            sheet_name = "peptides",index_col=0)

cmv_data['peptide'] = peptide_info.loc[cmv_data.index,'Peptide']

# file to save all data
tcr_data_all = pd.DataFrame(columns=['tcr_name', 'va', 'vb', 'cdr3a', 'cdr3b',
                                     'trav', 'traj', 'trbv',
       'trbd', 'trbj', 'assay', 'tcr_source_organism', 'index_peptide',
       'index_peptide_activity', 'mhc', 'pmid', 'peptide_type',
       'peptide','peptide_activity'])

# Separate TCRs
TCR_list = list(tcr_info['TCR id'])

for TCR in TCR_list:
    # Record TCR, MHC, organism etc info
    tcr_seq = tcr_info[tcr_info['TCR id']==TCR]
    tcr_data = pd.DataFrame({
    'tcr_name':[TCR],
    'va':list(tcr_seq['full_alpha_aa']),
    'vb':list(tcr_seq['full_beta_aa']),
    'cdr3a':list(tcr_seq['cdr3_a_aa']),
    'cdr3b':list(tcr_seq['cdr3_b_aa']),
    'trav':list(tcr_seq['TRAV']),
    'traj':list(tcr_seq['TRAJ']),
    'trbv':list(tcr_seq['TRBV']),
    'trbd':list(tcr_seq['TRBD']),
    'trbj':list(tcr_seq['TRBJ']),
    'assay':["NFAT luminescence"],
    'tcr_source_organism':["human"],
    'index_peptide':list(tcr_seq['Epitope']),
    'index_peptide_activity':list(cmv_data.loc[
        cmv_data.peptide=="NLVPMVATV",
        TCR.replace("-","_")[4:]]), #edit TCR name format to match data
    'mhc':["HLA-A*02:01"],
    'pmid':["39151427"],
    'peptide_type':["viral"]
    })

    # Get mutant nfat for the TCR
    mutant_nfat = cmv_data[['peptide',
                            TCR.replace("-","_")[4:]]].reset_index(drop=True) #edit TCR name format to match data
    mutant_nfat.columns = ['peptide','peptide_activity']
    
    # Merge all data columns
    tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_nfat))
                            ].reset_index(drop=True)
    tcr_data = pd.concat([tcr_data,mutant_nfat],axis=1)
    
    # add to data
    tcr_data_all = pd.concat([tcr_data_all,tcr_data])

# Export data
tcr_data_all.to_excel('epitope_data_cmv.xlsx', index=False)













