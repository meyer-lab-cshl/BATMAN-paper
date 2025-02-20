import numpy as np
import pandas as pd

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['TIL1383I'],
'va':["MTLSTLSLAKTTQPISMDSYEGQEVNITCSHNNIATNDYITWYQQFPSQGPRFIIQGYKTKVTNEVASLFIPADRKSSTLSLPRVSLSDTAVYYCLVALNYGGSQGNLIFGKGTKLSVKPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS"],
'vb':["MGHMDAGITQSPRHKVTETGTPVTLRCHQTENHRYMYWYRQDPGHGLRLIHYSYGVKDTDKGEVSDGYSVSRSKTEDFLLTLESATSSQTSVYFCAISPTEEGGLIFPGNTIYFGEGSWLTVVEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD"],
'cdr3a':["CLVALNYGGSQGNLIF"],
'cdr3b':["CAISPTEEGGLIFPGNTIYF"],
'trav':["4*01"],
'traj':["42*01"],
'trbv':["10-3*04"],
'trbd':["NA"],
'trbj':["1-3*01"],
'assay':["ELISA"],
'tcr_source_organism':["human"],
'index_peptide':["YMDGTMSQV"],
'index_peptide_activity':[3592.14],
'mhc':["HLA-A*02:01"],
'pmid':["36424374"],
'peptide_type':["cancer antigen"]
})



# Get mutant IFNg replicate
mutant_il2_raw = pd.read_excel("TIL1383I-IL2.xlsx",index_col=0,
                            skiprows=1)
index_peptide = "YMDGTMSQV"
mutant_il2 = pd.DataFrame(columns=['peptide','peptide_activity'])

# Extract mutation positions and mutated AAs
aa_order = mutant_il2_raw.columns
positions = mutant_il2_raw.index
positions = np.array([s[1] for s in positions]).astype(int)


for position_index in np.arange(len(positions)):
    for aa_index in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[positions[position_index]-1] = aa_order[aa_index]
        peptide = ''.join(peptide)
        peptide_activity = mutant_il2_raw.iloc[position_index,aa_index]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_il2 = pd.concat([mutant_il2,data_add])

mutant_il2 = mutant_il2.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_il2))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_il2],axis=1)

# Export data
tcr_data.to_excel('epitope_data_TIL1383I.xlsx', index=False)




















