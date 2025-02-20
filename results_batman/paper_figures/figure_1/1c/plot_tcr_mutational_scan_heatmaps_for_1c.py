'''Code to generate TCR activation mutational scan heatmaps'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting order for substituted amino acids (interfacial hydrophobicity order)
# Ref: https://www.annualreviews.org/content/journals/10.1146/annurev.biophys.28.1.319
amino_acid_list = ['W','F','Y','L','I','C','M','G','V','S','T','A','N','P',
                   'Q','R','H','K','D','E']

# Load full mutant epitope activity data for MHCI and MHCII
mutant_activity_data_mhci = pd.read_excel(
    ''.join(['../../../tcr_epitope_datasets/mutational_scan_datasets/',
             'TCR_pMHCI_mutational_scan_database.xlsx']))

mutant_activity_data_mhcii = pd.read_excel(
    ''.join(['../../../tcr_epitope_datasets/mutational_scan_datasets/',
             'TCR_pMHCII_mutational_scan_database.xlsx']))

mutant_activity_data = pd.concat([mutant_activity_data_mhci,
                                  mutant_activity_data_mhcii]).reset_index(drop=True)

# List of unique index peptides available in the dataset
index_peptide_list = np.unique(mutant_activity_data['index_peptide'].values)

# Get the heatmap of mutated peptide activity for for all TCRs for a given epitope
for index_peptide in index_peptide_list:
    
    # List TCR names available in the dataset for the selected index peptide
    TCR_name_list = np.unique(mutant_activity_data[
        mutant_activity_data.index_peptide==index_peptide]['tcr'].values)
    
    # Create array to store mutant peptide activity data for TCRs in above list
    heatmap_data = np.zeros((len(index_peptide)*19+1,len(TCR_name_list)))
    
    # Create mask array to store which mutations are missing
    mask = np.ones((len(index_peptide)*19+1,len(TCR_name_list)))
    
    # Store activity data for mutant peptides for a given TCR
    for TCR_index in np.arange(len(TCR_name_list)):
        
        TCR = TCR_name_list[TCR_index] # Name of TCR
        mutant_data = mutant_activity_data[
            mutant_activity_data.tcr==TCR][['peptide','peptide_activity']] #1 Hamming mutants
        mutant_data = mutant_data.set_index('peptide')
                
        # Normalize mutant activities to 1
        mutant_data['peptide_activity'] = mutant_data['peptide_activity']/max(
            mutant_data['peptide_activity'])
        
        # ordered mutant sequences
        mutant_list = np.full((19,len(index_peptide)),index_peptide)
        for p in np.arange(len(index_peptide)):
            # List all mutant AA at peptide position p
            mutant_aa = [a for a in amino_acid_list if a!=index_peptide[p]]
            for aa in np.arange(19):
                # Create mutant sequence
                mutant = list(index_peptide)
                mutant[p] = mutant_aa[aa]
                mutant = ''.join(mutant)
                mutant_list[aa,p] = mutant
        # Make a linear array        
        mutant_list = np.reshape(mutant_list,(19*len(index_peptide),1),order='F')
        # Add index peptide at beginning
        mutant_list = np.insert(mutant_list, 0, index_peptide)
        
        # Reindex mutant activity data according to the order 
        # (returns NaN for missing data)
        mutant_data = mutant_data.reindex(list(mutant_list))
        
        # Record data
        heatmap_data[:,TCR_index] = mutant_data['peptide_activity']
        
        # create mask where data is absent
        mask[:,TCR_index] = np.isnan(heatmap_data)[:,0]
        

    # Fill NaN values with 0 (they won't be plotted)
    heatmap_data = np.nan_to_num(heatmap_data, nan=0)
    
    # Plotting
    if len(TCR_name_list) > 1: #cluster TCRs if there are multiple of them
    
        # Dendrogram to cluster TCRs based on similar specificities
        plot = sns.clustermap(heatmap_data, cmap='plasma', row_cluster=False, 
                   figsize=[len(TCR_name_list)*0.5,len(index_peptide)*19+1],
                   vmin=0, vmax=1, yticklabels=0, xticklabels=0, mask=mask) 
        plot.ax_col_dendrogram.set_visible(False)       
     
        plt.show() 
        
        plt.savefig('heatmaps_by_index_peptide/%s.pdf' % index_peptide, format='pdf')
        
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
        
        plt.savefig('heatmaps_by_index_peptide/%s.pdf' % index_peptide, format='pdf')
        
    
       
        
