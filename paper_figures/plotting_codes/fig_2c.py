# Code to generate AA distance matrices

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Default order for substituted amino acids (alphabetic order)
amino_acid_list = np.array(['A','C','D','E','F','G','H','I','K','L','M','N',
                            'P','Q','R','S','T','V','W','Y'])

# Plotting AA distance heatmap for symmetric matrix
# Load amino acid distance matrix
AA_distance_matrix = pd.read_excel('../data/AA_matrix_ratio_symm.xlsx').to_numpy()

# order of AAs for good clustering
AA_order=[5,12,0,16,15,11,2,3,8,14,13,6,19,18,4,9,7,17,10,1]

#Reorder AA matrix and AA list
AA_distance_matrix_clustered = np.zeros((20,20))
for aa1 in np.arange(0,20,1):
    for aa2 in np.arange(0,20,1):    
        AA_distance_matrix_clustered[aa1,aa2] = AA_distance_matrix[
                                                  AA_order[aa1],AA_order[aa2]]

amino_acid_list = amino_acid_list[AA_order]

# Create mask array to plot only the lower triangular part
mask = np.triu((AA_distance_matrix)) + np.identity(20)

# plot triangular heatmap
sns.set(font_scale=1.2)
plot = sns.heatmap(AA_distance_matrix_clustered, cmap="Blues",
            yticklabels=amino_acid_list, 
            xticklabels=amino_acid_list,mask=mask,square=1) 
plt.yticks(rotation = 0)
plt.xticks(rotation = 0)
plt.show()

#plt.savefig('../figures/fig2/aa_matrix_symm.pdf', format='pdf')

# Default order for substituted amino acids (alphabetic order)
amino_acid_list = np.array(['A','C','D','E','F','G','H','I','K','L','M','N',
                            'P','Q','R','S','T','V','W','Y'])

# Plotting AA distance heatmap for antisymmetric part of matrix
# Load amino acid distance matrix
AA_distance_matrix = pd.read_excel('../data/AA_matrix_antisymm.xlsx').to_numpy()

# order of AAs for good clustering
AA_order=[5,12,0,16,15,11,2,3,8,14,13,6,19,18,4,9,7,17,10,1]

#Reorder AA matrix and AA list
AA_distance_matrix_clustered = np.zeros((20,20))
for aa1 in np.arange(0,20,1):
    for aa2 in np.arange(0,20,1):    
        AA_distance_matrix_clustered[aa1,aa2] = AA_distance_matrix[
                                                  AA_order[aa1],AA_order[aa2]]

amino_acid_list = amino_acid_list[AA_order]

# Create mask array to plot only the lower triangular part
mask = np.tril((AA_distance_matrix)) + np.identity(20)

# plot triangular heatmap
sns.set(font_scale=1.2)
plot = sns.heatmap(AA_distance_matrix_clustered, cmap="coolwarm_r",
            yticklabels=amino_acid_list, 
            xticklabels=amino_acid_list,mask=mask,square=1,vmin=-9,vmax=9) 
plt.yticks(rotation = 0)
plt.xticks(rotation = 0)
plt.show()

plt.savefig('../figures/fig2/aa_matrix_inferred.pdf', format='pdf')

