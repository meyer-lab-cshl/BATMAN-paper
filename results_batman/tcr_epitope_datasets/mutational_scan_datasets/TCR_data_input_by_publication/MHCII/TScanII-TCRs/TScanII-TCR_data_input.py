import numpy as np
import pandas as pd

'''TCR F5, F13, and F24'''
tcr_name_list = ['F5','F13','F24']
cdr3b_list = ["CASSGLAGGMDEQFF","CASSRRTSGGTDTQYF","CASSRLAGGMDEQFF"]
trbv_list = ["2-1","2-2","2-1"]

# Get mutant FC data
mutant_fc_raw = pd.read_excel("mmc1.xlsx", sheet_name="Table S1B",skiprows=1)

# Initialize df
all_data=pd.DataFrame()

for tcr in np.arange(3):

    # Record TCR, MHC, organism etc info
    tcr_data = pd.DataFrame({
    'tcr_name':[tcr_name_list[tcr]],
    'va':["NA"],
    'vb':["NA"],
    'cdr3a':["CAFKAAGNKLTF"],
    'cdr3b':[cdr3b_list[tcr]],
    'trav':["24-1"],
    'traj':["NA"],
    'trbv':[trbv_list[tcr]],
    'trbd':["NA"],
    'trbj':["NA"],
    'assay':["TScan II"],
    'tcr_source_organism':["human"],
    'index_peptide':["FRDYVDRFYKTLRAEQASQE"],
    'mhc':["DRB1∗11:01"],
    'pmid':["38016469"],
    'peptide_type':["viral"]
    })   
    
    index_peptide_extended = 'SILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKA'
    aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                     'S','T','V','W','Y']
    mutant_fc = pd.DataFrame(columns=['peptide','peptide_activity'])
    
    for position in np.arange(12,32,1):
        for aa in aa_list:
            
            # Create mutant seq
            peptide = list(index_peptide_extended)
            peptide[position] = aa
            peptide = ''.join(peptide)
            
            # search for mutant seq in data
            peptide_activity = mutant_fc_raw[mutant_fc_raw['Peptide Sequence']==peptide]
            peptide_activity = np.mean(peptide_activity.iloc[:,2+tcr].to_numpy())
            
            # Add data
            data_add = pd.DataFrame({'peptide':[peptide[12:32]],
                                     'peptide_activity':[peptide_activity]})
            mutant_fc = pd.concat([mutant_fc,data_add])
        
    # Remove duplicates
    mutant_fc = mutant_fc.drop_duplicates().reset_index(drop=True)
    
    # minmax normalize data
    mutant_fc.peptide_activity = (
        mutant_fc.peptide_activity-np.min(mutant_fc.peptide_activity))/(
            np.max(mutant_fc.peptide_activity)-np.min(mutant_fc.peptide_activity))

    # Merge all data columns
    tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_fc))
                            ].reset_index(drop=True)
    tcr_data = pd.concat([tcr_data,mutant_fc],axis=1)
    
    # add to full data
    all_data = pd.concat([all_data,tcr_data])

'''TCR 3598-2'''
# Get mutant FC data
mutant_fc_raw = pd.read_excel("mmc2.xlsx", sheet_name="Table S2B",skiprows=1)



# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['TCR-3598-2'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CALRDSGGGADGLTF"],
'cdr3b':["CASSVMTGLNTEAFF"],
'trav':["9-2"],
'traj':["45"],
'trbv':["2"],
'trbd':["NA"],
'trbj':["1-1"],
'assay':["TScan II"],
'tcr_source_organism':["human"],
'index_peptide':["LPVPGVLLKEFTVSGNILTI"],
'mhc':["DRB1∗04:01"],
'pmid':["38016469"],
'peptide_type':["cancer antigen"]
})   

index_peptide_extended = 'AMPFATPMEAELARRSLAQDAPPLPVPGVLLKEFTVSGNILTIRLTAADHRQLQLS'
aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_fc = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(23,43,1):
    for aa in aa_list:
        
        # Create mutant seq
        peptide = list(index_peptide_extended)
        peptide[position] = aa
        peptide = ''.join(peptide)
        
        # search for mutant seq in data
        peptide_activity = mutant_fc_raw[mutant_fc_raw['Peptide Sequence']==peptide]
        peptide_activity = np.mean(peptide_activity.iloc[:,2].to_numpy())
        
        # Add data
        data_add = pd.DataFrame({'peptide':[peptide[23:43]],
                                 'peptide_activity':[peptide_activity]})
        mutant_fc = pd.concat([mutant_fc,data_add])
    
# Remove duplicates
mutant_fc = mutant_fc.drop_duplicates().reset_index(drop=True)

# minmax normalize data
mutant_fc.peptide_activity = (
    mutant_fc.peptide_activity-np.min(mutant_fc.peptide_activity))/(
        np.max(mutant_fc.peptide_activity)-np.min(mutant_fc.peptide_activity))

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_fc))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_fc],axis=1)

# add to full data
all_data = pd.concat([all_data,tcr_data])

# Export data
all_data.to_excel('epitope_data_tscan_ii.xlsx', index=False)













