import numpy as np
import pandas as pd

'''TCR 3598-2'''
# Get mutant FC data
mutant_fc_raw = pd.read_excel("jpeg2excel.xlsx",index_col=0)

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['MBP-TCR'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["T cell proliferation"],
'tcr_source_organism':["mouse"],
'index_peptide':["ASQKRPSQRSK"],
'mhc':["I-Au"],
'pmid':["10490973"],
'peptide_type':["autoimmunity"]
})   

index_peptide = 'ASQKRPSQRSK'
aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_fc = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(11):
    for aa in np.arange(20):
        
        # Create mutant seq
        peptide = list(index_peptide)
        peptide[position] = aa_list[aa]
        peptide = ''.join(peptide)
        
        # activity for mutant seq in data
        peptide_activity = mutant_fc_raw.iloc[aa+1,position]
        
        # Add data
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_fc = pd.concat([mutant_fc,data_add])

# Average index peptide activity
index_activity1 = mutant_fc_raw.iloc[0,:].to_numpy()
index_activity2 = mutant_fc.peptide_activity[mutant_fc.peptide==index_peptide].to_numpy()
index_activity_mean = np.mean(np.concatenate([index_activity1,index_activity2]))

mutant_fc.peptide_activity[mutant_fc.peptide==index_peptide] = index_activity_mean

# Remove duplicates
mutant_fc = mutant_fc.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_fc))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_fc],axis=1)

# Export data
tcr_data.to_excel('epitope_data_MBP-TCR.xlsx', index=False)













