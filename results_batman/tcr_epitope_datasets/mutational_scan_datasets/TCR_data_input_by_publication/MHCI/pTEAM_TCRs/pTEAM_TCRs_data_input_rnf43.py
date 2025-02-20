import numpy as np
import pandas as pd

# Read neoantigen data
neoantigen_data = pd.read_excel("Affinity_prediction_rnf43_repertoire.xlsx",
                            sheet_name = "Individual APL screening")
tcr_info = pd.read_excel("Affinity_prediction_rnf43_repertoire.xlsx",
                            sheet_name = "Sequences",skiprows=1)

# file to save all data
tcr_data_all = pd.DataFrame(columns=['tcr_name', 'va', 'vb', 'cdr3a', 'cdr3b',
                                     'trav', 'traj', 'trbv',
       'trbd', 'trbj', 'assay', 'tcr_source_organism', 'index_peptide',
       'index_peptide_activity', 'mhc', 'pmid', 'peptide_type',
       'peptide','peptide_activity'])
# Separate TCRs
TCR_list = list(tcr_info.TCR)

for TCR in TCR_list:
    # Record TCR, MHC, organism etc info
    tcr_seq = tcr_info[tcr_info.TCR==TCR]
    tcr_data = pd.DataFrame({
    'tcr_name':[TCR],
    'va':["NA"],
    'vb':["NA"],
    'cdr3a':list(tcr_seq['CDR3α']),
    'cdr3b':list(tcr_seq['CDR3β']),
    'trav':list(tcr_seq['TRAV']),
    'traj':list(tcr_seq['TRAJ']),
    'trbv':list(tcr_seq['TRBV']),
    'trbd':list(tcr_seq['TRBD']),
    'trbj':list(tcr_seq['TRBJ']),
    'assay':["NFAT luminescence"],
    'tcr_source_organism':["human"],
    'index_peptide':["VPSVWRSSL"],
    'index_peptide_activity':list(neoantigen_data.loc[
        neoantigen_data.Peptide=="VPSVWRSSL",TCR]),
    'mhc':list(tcr_seq['Restriction']),
    'pmid':["39151427"],
    'peptide_type':["neoantigen"]
    })

    # Get mutant nfat for the TCR
    mutant_nfat = neoantigen_data[['Peptide',TCR]]
    mutant_nfat.columns = ['peptide','peptide_activity']
    
    # Merge all data columns
    tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_nfat))
                            ].reset_index(drop=True)
    tcr_data = pd.concat([tcr_data,mutant_nfat],axis=1)
    
    # add to data
    tcr_data_all = pd.concat([tcr_data_all,tcr_data])

# Export data
tcr_data_all.to_excel('epitope_data_rnf43.xlsx', index=False)













