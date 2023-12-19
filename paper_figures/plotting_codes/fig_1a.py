# Code to generate mutational scan heatmaps for all TCRs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting order for substituted amino acids (alphabetic order)
amino_acid_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']

# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('../data/TCR_epitope_database.xlsx');

# List of unique index peptides available in the dataset
index_peptide_list = np.unique(mutant_activity_data['index_peptide'].values)

# Get the heatmap of mutated peptide activity for for all TCRs for a given epitope
for index_peptide in index_peptide_list:
    
    # List TCR names available in the dataset for the selected index peptide
    TCR_name_list = np.unique(mutant_activity_data[
        mutant_activity_data.index_peptide==index_peptide]['tcr_name'].values)
    
    # Create array to store mutant peptide activity data for TCRs in above list
    heatmap_data = np.zeros((len(index_peptide)*19+1,len(TCR_name_list)))
    
    # Create mask array to store which mutations are missing
    mask = np.zeros((len(index_peptide)*19+1,len(TCR_name_list))) + 1
    
    # Store activity data for mutant peptides for a given TCR
    for TCR_index in np.arange(0,len(TCR_name_list),1):
        
        TCR = TCR_name_list[TCR_index] # Name of TCR
        mutant_sequences = mutant_activity_data[
            mutant_activity_data.tcr_name==TCR]['peptide'].values #1 Hamming mutants
        
        mutant_activities = mutant_activity_data[
            mutant_activity_data.tcr_name==TCR]['peptide_activity'].values
        
        # Normalize mutant activities to 1
        mutant_activities = mutant_activities/max(mutant_activities)
        
        # Store mutant activity data for mutants
        for mutant_index in np.arange(0,len(mutant_sequences),1):
            
            mutant = mutant_sequences[mutant_index] #mutant sequence
            
            # Separate the index-peptide to put it at the end
            if mutant != index_peptide:
                
                # Location of mutation
                mutation_position = np.argwhere(np.char.compare_chararrays(list(
                    index_peptide),list(mutant),'!=','true'))[0][0]
                
                #index of mutated AA in AA list without the WT AA                
                mutated_AA_index = list(np.setdiff1d(
                    np.array(amino_acid_list),
                    [index_peptide[mutation_position]])
                    ).index(mutant[mutation_position])
                
                # Store mutant activity data using a linear index
                heatmap_data[mutation_position*19 + mutated_AA_index, 
                             TCR_index] = mutant_activities[mutant_index]
                
                # Store 0 in mask to plot the data
                mask[mutation_position*19 + mutated_AA_index, 
                             TCR_index] = 0
                
            
            # Put the index peptide at the end   
            if mutant == index_peptide:
                heatmap_data[-1, TCR_index] = mutant_activities[mutant_index]
                mask[-1, TCR_index] = 0
    
    # Plotting
    if len(TCR_name_list) > 1: #cluster TCRs if there are multiple of them
    
        # Dendrogram to cluster TCRs based on similar specificities
        plot = sns.clustermap(heatmap_data, cmap='plasma', row_cluster=False, 
                   figsize=[len(TCR_name_list)*0.5,len(index_peptide)*19+1],
                   vmin=0, vmax=1, yticklabels=0, xticklabels=0, mask=mask) 
        plot.ax_col_dendrogram.set_visible(False)       
     
        plt.show() 
        
        plt.savefig('../figures/fig1/fig1a/%s.pdf' % index_peptide, format='pdf')
        
    if len(TCR_name_list) == 1: # Just one TCR available for the index peptide
    
        # Duplicate data to do a fake clustering to have the same style as above
        heatmap_data = np.repeat(heatmap_data, 2, axis=1)
        mask = np.repeat(mask, 2, axis=1)
        
        # Dendrogram to cluster TCRs based on similar specificities
        plot = sns.clustermap(heatmap_data, cmap='plasma', row_cluster=False, 
                   figsize=[len(TCR_name_list),len(index_peptide)*19+1],
                   vmin=0, vmax=1, yticklabels=0, xticklabels=0, mask=mask) 
        
        plot.ax_col_dendrogram.set_visible(False)
        
        plt.show() 
        
        plt.savefig('../figures/fig1/fig1a/%s.pdf' % index_peptide, format='pdf')
        
    
       
        
