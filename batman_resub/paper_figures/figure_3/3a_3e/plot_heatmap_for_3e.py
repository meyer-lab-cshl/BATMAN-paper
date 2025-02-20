import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting order for substituted amino acids (interfacial hydrophobicity order)
# Ref: https://www.annualreviews.org/content/journals/10.1146/annurev.biophys.28.1.319
amino_acid_list = ['W','F','Y','L','I','C','M','G','V','S','T','A','N','P',
                   'Q','R','H','K','D','E']

# Load amino acid distance matrix
AA_distance_matrix = pd.read_csv('inferred_parameters/inferred_aa_matrix_9_mers.csv',
                                 index_col=0)

AA_distance_matrix = AA_distance_matrix.loc[amino_acid_list,amino_acid_list]

# Load BLOSUM100 matrix
# Load amino acid distance matrix
blosum100 = pd.read_csv('inferred_parameters/blosum100.csv',
                                 index_col=0)

blosum100 = blosum100.loc[amino_acid_list,amino_acid_list]

# plot heatmap

# 2.3 is the average ration with blosum in hyperparameters
sns.set(font_scale=1.2)
plot = sns.heatmap(AA_distance_matrix/(2.3*blosum100), cmap="RdBu",
            yticklabels=amino_acid_list, 
            xticklabels=amino_acid_list,square=1,center=1,vmin=0,vmax=2) 
plt.yticks(rotation = 0)
plt.xticks(rotation = 0)
plt.show()

plt.savefig('aa_matrix_for_3e.pdf', format='pdf')

# Save data
df = AA_distance_matrix/(2.3*blosum100)
df.to_csv('raw_data_fig_3e.csv')


