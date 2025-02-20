import numpy as np
import pandas as pd

index_peptide = 'NLVPMVATV'

# Get mutant IFNg replicate
mutant_data = pd.read_excel("TScan_figs5_lrg.jpg.xlsx",
                            index_col=0)

aa_order = mutant_data.index

mutant_activity = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(9):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = mutant_data.iloc[aa,position]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_activity = pd.concat([mutant_activity,data_add])

mutant_activity = mutant_activity.drop_duplicates().reset_index(drop=True)

# Export data
mutant_activity.to_csv('NLV_mutant_activity.csv')




















