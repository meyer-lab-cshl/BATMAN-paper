import numpy as np
import pandas as pd

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['47BE7'],
'va':[""],
'vb':[""],
'cdr3a':["CAVSNYNVLYF"],
'cdr3b':["CASSQEPGGYAEQFF"],
'trav':["7-1"],
'traj':["21"],
'trbv':["2"],
'trbd':["2"],
'trbj':["2-1"],
'assay':["IFNg"],
'tcr_source_organism':["mouse"],
'index_peptide':["YGFRNVVHI"],
'index_peptide_activity':[1],
'mhc':["H2-Db"],
'pmid	':["36778273"],
'peptide_type':["neoantigen"]
})



# Get mutant IFNg replicate
mutant_ifng_raw = pd.read_excel("41467_2024_46367_MOESM5_ESM.xlsx",
                            sheet_name = "Fig. 4F",
                            index_col=0,skiprows=2)
index_peptide = "YGFRNVVHI"
aa_order = mutant_ifng_raw.index
mutant_ifng = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(9):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = np.nanmean(mutant_ifng_raw.iloc[aa,
                                            3*position:3*position+3])
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_ifng = pd.concat([mutant_ifng,data_add])

mutant_ifng = mutant_ifng.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_ifng))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_ifng],axis=1)

# Export data
tcr_data.to_excel('epitope_data_47BE7.xlsx', index=False)




















