import numpy as np
import pandas as pd

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['a3a'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["TCR-MAP"],
'tcr_source_organism':["human"],
'index_peptide':["EVDPIGHLY"],
'index_peptide_activity':[1],
'mhc':["HLA-A*01:01"],
'pmid':["38956325"],
'peptide_type':["cancer antigen"]
})

# Get mutant IFNg replicate
mutant_fc_raw = pd.read_excel("41587_2024_2248_MOESM10_ESM.xlsx",
                            sheet_name = "Scren results for Figure 6b",
                            skiprows=1)
index_peptide = "EVDPIGHLY"

mutant_fc = pd.DataFrame(columns=['peptide','peptide_activity'])
for row in np.arange(len(mutant_fc_raw)):
    peptide = list(index_peptide)
    aa = mutant_fc_raw.iloc[row,0]
    position = mutant_fc_raw.iloc[row,1]
    if type(position)==int: #only sample positions 1-9
        peptide[position-1] = aa
        peptide = ''.join(peptide)
        peptide_activity = mutant_fc_raw.iloc[row,2]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_fc = pd.concat([mutant_fc,data_add])

# Remove duplicates
mutant_fc = mutant_fc.drop_duplicates().reset_index(drop=True)

# Index peptide activity
mutant_fc.iloc[np.isnan(mutant_fc.peptide_activity),1] = 1

# minmax normalize data
mutant_fc.peptide_activity = (
    mutant_fc.peptide_activity-np.min(mutant_fc.peptide_activity))/(
        np.max(mutant_fc.peptide_activity)-np.min(mutant_fc.peptide_activity))

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_fc))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_fc],axis=1)

# Export data
tcr_data.to_excel('epitope_data_a3a.xlsx', index=False)













