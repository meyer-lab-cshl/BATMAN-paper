# Code to generate mutational scan heatmaps for one TCR

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting order for substituted amino acids (interfacial hydrophobicity order)
# Ref: https://www.annualreviews.org/content/journals/10.1146/annurev.biophys.28.1.319
amino_acid_list = ['W','F','Y','L','I','C','M','G','V','S','T','A','N','P',
                   'Q','R','H','K','D','E']

# Load full mutant epitope activity data for MHCI 
mutant_activity_data = pd.read_excel(
    ''.join(['../../../tcr_epitope_datasets/mutational_scan_datasets/',
             'TCR_pMHCI_mutational_scan_database.xlsx']))

# Select TCR and index to plot
TCR = 'NLV2'
index_peptide = list('NLVPMVATV')
    
# Create array to store mutant peptide activity data for all mutations
heatmap_data = np.zeros((len(amino_acid_list), len(index_peptide)))

mutant_data = mutant_activity_data[
    mutant_activity_data.tcr==TCR][['peptide','peptide_activity']] #1 Hamming mutants
mutant_data = mutant_data.set_index('peptide')
        
# Normalize mutant activities to 1
mutant_data['peptide_activity'] = mutant_data['peptide_activity']/max(
    mutant_data['peptide_activity'])
        
# Store mutant activity data for mutants
for mutated_amino_acid in np.arange(0,len(amino_acid_list),1):
    for position in np.arange(0,len(index_peptide),1):
        mutant = list(index_peptide) #initialize
        mutant[position] = amino_acid_list[mutated_amino_acid] #put mutation
        mutant = ''.join(mutant) #convert list to string
        #Find mutant in sequence list and record activity
        heatmap_data[mutated_amino_acid,position] = mutant_data.loc[mutant,
                                                                    'peptide_activity']    
    
           
# Plotting peptide activity heatmap
fig, ax = plt.subplots(figsize=(12, 5))
sns.set(font_scale=3) 
plot = sns.heatmap(heatmap_data.T, cmap='plasma', linewidths=0.3,
                   vmin=0, vmax=1, yticklabels=list(index_peptide),
                   xticklabels = amino_acid_list) 
plt.yticks(rotation = 0)
plt.show() 

plt.savefig('fig1b_heatmap_schematic.pdf', format='pdf')
