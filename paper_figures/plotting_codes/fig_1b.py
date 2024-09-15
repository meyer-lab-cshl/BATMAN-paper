# Code to generate mutational scan heatmaps for one TCR and the AA distance matrix

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting order for substituted amino acids (alphabetic order)
amino_acid_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']

# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('../data/TCR_epitope_database.xlsx');

# Select TCR and index to plot
TCR = 'NLV2'
index_peptide = list('NLVPMVATV')
    
# Create array to store mutant peptide activity data for all mutations
heatmap_data = np.zeros((len(amino_acid_list), len(index_peptide)))

# Load mutant sequences
mutant_sequences = mutant_activity_data[
    mutant_activity_data.tcr_name==TCR]['peptide'].values

# Load mutant activity data
mutant_activities = mutant_activity_data[
    mutant_activity_data.tcr_name==TCR]['peptide_activity'].values

# Normalize mutant activities to 1
mutant_activities = mutant_activities/max(mutant_activities)
        
# Store mutant activity data for mutants
for mutated_amino_acid in np.arange(0,len(amino_acid_list),1):
    for position in np.arange(0,len(index_peptide),1):
        mutant = ['A','A','A','A','A','A','A','A','A'] #initialize
        mutant[:] = index_peptide #initialize
        mutant[position] = amino_acid_list[mutated_amino_acid] #put mutation
        mutant = ''.join(mutant) #convert list to string
        #Find mutant in sequence list and record activity
        heatmap_data[mutated_amino_acid,position] = mutant_activities[np.where(
            mutant_sequences == mutant)]    
    
           
# Plotting peptide activity heatmap
fig, ax = plt.subplots(figsize=(5, 10))
sns.set(font_scale=3) 
plot = sns.heatmap(heatmap_data, cmap='plasma', linewidths=0.3,
                   vmin=0, vmax=1, xticklabels=0,
                   yticklabels = amino_acid_list) 
plt.yticks(rotation = 0)
plt.show() 

plt.savefig('../figures/fig1/fig1b/%s.pdf' % TCR, format='pdf')

# Plotting AA distance heatmap
# Load amino acid distance matrix
AA_distance_matrix = pd.read_excel('../data/aa_matrix.xlsx').to_numpy()
# Dendrogram to cluster AAs based on their distances
plot = sns.clustermap(AA_distance_matrix, cmap="Blues_r",
           square=1,yticklabels=0, xticklabels=0) 

plot.ax_col_dendrogram.set_visible(False)
plot.ax_row_dendrogram.set_visible(False)
plt.show()

plt.savefig('../figures/fig1/fig1b/aa_matrix.pdf', format='pdf')

        
    