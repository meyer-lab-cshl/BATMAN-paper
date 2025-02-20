import numpy as np
import pandas as pd

'''TCR A3V'''

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['A3V'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CALSEAGTYKYIF"],
'cdr3b':["CASSVAGGGQETQY"],
'trav':["19"],
'traj':["40"],
'trbv':["9"],
'trbd':["1"],
'trbj':["2-5"],
'assay':["NFAT-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["VVVGAVGVGK"],
'mhc':["HLA-A*03:01"],
'pmid':["39287991"],
'peptide_type':["cancer antigen"]
})

# Get mutant activity raw data
mutant_activity_raw = pd.read_excel("JCI175790.sdval.xlsx",
                            sheet_name = "Figure 2C",
                            skiprows=9)
# Get cysteine substitution values
mutant_activity_C = mutant_activity_raw.iloc[0:10,15].to_numpy().T

# Extract data from the 18X10 grid
mutant_activity_raw = mutant_activity_raw.iloc[0:18,2:12].to_numpy()

# Combine all mutant AAs
mutant_activity_raw = np.insert(mutant_activity_raw, 1, 
                                mutant_activity_C, axis=0)

index_peptide = "VVVGAVGVGK"

# Order of AAs for data
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
            'T','V','W','Y']

# Record data
mutant_activity = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(np.shape(mutant_activity_raw)[1]):
    wt_aa = index_peptide[position]
    aa_list = [a for a in aa_order if a!=wt_aa] #List of substitution AAs
    
    for iaa in np.arange(len(aa_list)):
        aa = aa_list[iaa]
        peptide = list(index_peptide)
        peptide[position] = aa
        peptide = ''.join(peptide)
        peptide_activity = mutant_activity_raw[iaa,position]
        data_add = pd.DataFrame({'peptide':[peptide],
                             'peptide_activity':[peptide_activity]})
        mutant_activity = pd.concat([mutant_activity,data_add])

# Add index peptide activity
data_add = pd.DataFrame({'peptide':[index_peptide],
                      'peptide_activity':[95.9]})
mutant_activity = pd.concat([mutant_activity,data_add])

# Min max activity for normalization 
min_activity = 2.43
max_activity = 92.64    


# minmax normalize data
mutant_activity.peptide_activity = (
    mutant_activity.peptide_activity-min_activity)/(
        max_activity - min_activity)
        
mutant_activity.peptide_activity[mutant_activity.peptide_activity<0] = 0        

mutant_activity = mutant_activity.reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_activity))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_activity],axis=1)

tcr_data_all = tcr_data.copy()



'''TCR A11Va'''
# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['A11Va'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CAVNPPDTGFQKLV"],
'cdr3b':["CASSLSFRQGLREQYF"],
'trav':["12-1"],
'traj':["8"],
'trbv':["28"],
'trbd':["1"],
'trbj':["2-7"],
'assay':["NFAT-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["VVVGAVGVGK"],
'mhc':["HLA-A*11:01"],
'pmid':["39287991"],
'peptide_type':["cancer antigen"]
})

# Get mutant activity raw data
mutant_activity_raw = pd.read_excel("JCI175790.sdval.xlsx",
                            sheet_name = "Figure 2C",
                            skiprows=9)

# Get cysteine substitution values
mutant_activity_C = mutant_activity_raw.iloc[30:40,15].to_numpy().T

# Extract data from the 18X10 grid
mutant_activity_raw = mutant_activity_raw.iloc[30:48,2:12].to_numpy()

# Combine all mutant AAs
mutant_activity_raw = np.insert(mutant_activity_raw, 1, 
                                mutant_activity_C, axis=0)

index_peptide = "VVVGAVGVGK"

# Order of AAs for data
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
            'T','V','W','Y']

# Record data
mutant_activity = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(np.shape(mutant_activity_raw)[1]):
    wt_aa = index_peptide[position]
    aa_list = [a for a in aa_order if a!=wt_aa] #List of substitution AAs
    
    for iaa in np.arange(len(aa_list)):
        aa = aa_list[iaa]
        peptide = list(index_peptide)
        peptide[position] = aa
        peptide = ''.join(peptide)
        peptide_activity = mutant_activity_raw[iaa,position]
        data_add = pd.DataFrame({'peptide':[peptide],
                             'peptide_activity':[peptide_activity]})
        mutant_activity = pd.concat([mutant_activity,data_add])

# Add index peptide activity
data_add = pd.DataFrame({'peptide':[index_peptide],
                      'peptide_activity':[87]})
mutant_activity = pd.concat([mutant_activity,data_add])

# Min max activity for normalization 
min_activity = 2.74
max_activity = 92.567    


# minmax normalize data
mutant_activity.peptide_activity = (
    mutant_activity.peptide_activity-min_activity)/(
        max_activity - min_activity)
        
mutant_activity.peptide_activity[mutant_activity.peptide_activity<0] = 0        

mutant_activity = mutant_activity.reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_activity))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_activity],axis=1)

tcr_data_all = pd.concat([tcr_data_all,tcr_data]).reset_index(drop=True)


'''TCR A11Vb'''

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['A11Vb'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CAVDKDGGYQKVTF"],
'cdr3b':["CSASPRAGQLSSYNSPLHF"],
'trav':["39"],
'traj':["13"],
'trbv':["20-1"],
'trbd':["1"],
'trbj':["1-6"],
'assay':["NFAT-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["VVVGAVGVGK"],
'mhc':["HLA-A*11:01"],
'pmid':["39287991"],
'peptide_type':["cancer antigen"]
})

# Get mutant activity raw data
mutant_activity_raw = pd.read_excel("JCI175790.sdval.xlsx",
                            sheet_name = "Figure 2C",
                            skiprows=9)

# Get cysteine substitution values
mutant_activity_C = mutant_activity_raw.iloc[60:70,15].to_numpy().T

# Extract data from the 18X10 grid
mutant_activity_raw = mutant_activity_raw.iloc[60:78,2:12].to_numpy()

# Combine all mutant AAs
mutant_activity_raw = np.insert(mutant_activity_raw, 1, 
                                mutant_activity_C, axis=0)

index_peptide = "VVVGAVGVGK"

# Order of AAs for data
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
            'T','V','W','Y']

# Record data
mutant_activity = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(np.shape(mutant_activity_raw)[1]):
    wt_aa = index_peptide[position]
    aa_list = [a for a in aa_order if a!=wt_aa] #List of substitution AAs
    
    for iaa in np.arange(len(aa_list)):
        aa = aa_list[iaa]
        peptide = list(index_peptide)
        peptide[position] = aa
        peptide = ''.join(peptide)
        peptide_activity = mutant_activity_raw[iaa,position]
        data_add = pd.DataFrame({'peptide':[peptide],
                             'peptide_activity':[peptide_activity]})
        mutant_activity = pd.concat([mutant_activity,data_add])

# Add index peptide activity
data_add = pd.DataFrame({'peptide':[index_peptide],
                      'peptide_activity':[86.9]})
mutant_activity = pd.concat([mutant_activity,data_add])

# Min max activity for normalization 
min_activity = 2.43
max_activity = 97.38571429

# minmax normalize data
mutant_activity.peptide_activity = (
    mutant_activity.peptide_activity-min_activity)/(
        max_activity - min_activity)
        
mutant_activity.peptide_activity[mutant_activity.peptide_activity<0] = 0        

mutant_activity = mutant_activity.reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_activity))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_activity],axis=1)

tcr_data_all = pd.concat([tcr_data_all,tcr_data]).reset_index(drop=True)


'''TCR A11Vc'''

# Record TCR, MHC, organism etc info
tcr_data = pd.DataFrame({
'tcr_name':['A11Vc'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CATDPGGFKTIF"],
'cdr3b':["CASSLYGGSISYEQYF"],
'trav':["17"],
'traj':["9"],
'trbv':["11-2"],
'trbd':["1"],
'trbj':["2-7"],
'assay':["NFAT-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["VVVGAVGVGK"],
'mhc':["HLA-A*11:01"],
'pmid':["39287991"],
'peptide_type':["cancer antigen"]
})

# Get mutant activity raw data
mutant_activity_raw = pd.read_excel("JCI175790.sdval.xlsx",
                            sheet_name = "Figure 2C",
                            skiprows=9)

# Get cysteine substitution values
mutant_activity_C = mutant_activity_raw.iloc[90:100,15].to_numpy().T

# Extract data from the 18X10 grid
mutant_activity_raw = mutant_activity_raw.iloc[90:108,2:12].to_numpy()

# Combine all mutant AAs
mutant_activity_raw = np.insert(mutant_activity_raw, 1, 
                                mutant_activity_C, axis=0)

index_peptide = "VVVGAVGVGK"

# Order of AAs for data
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
            'T','V','W','Y']

# Record data
mutant_activity = pd.DataFrame(columns=['peptide','peptide_activity'])

for position in np.arange(np.shape(mutant_activity_raw)[1]):
    wt_aa = index_peptide[position]
    aa_list = [a for a in aa_order if a!=wt_aa] #List of substitution AAs
    
    for iaa in np.arange(len(aa_list)):
        aa = aa_list[iaa]
        peptide = list(index_peptide)
        peptide[position] = aa
        peptide = ''.join(peptide)
        peptide_activity = mutant_activity_raw[iaa,position]
        data_add = pd.DataFrame({'peptide':[peptide],
                             'peptide_activity':[peptide_activity]})
        mutant_activity = pd.concat([mutant_activity,data_add])

# Add index peptide activity
data_add = pd.DataFrame({'peptide':[index_peptide],
                      'peptide_activity':[87.5]})
mutant_activity = pd.concat([mutant_activity,data_add])

# Min max activity for normalization 
min_activity = 4.33
max_activity = 90.26667

# minmax normalize data
mutant_activity.peptide_activity = (
    mutant_activity.peptide_activity-min_activity)/(
        max_activity - min_activity)
        
mutant_activity.peptide_activity[mutant_activity.peptide_activity<0] = 0        

mutant_activity = mutant_activity.reset_index(drop=True)

# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_activity))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_activity],axis=1)

tcr_data_all = pd.concat([tcr_data_all,tcr_data]).reset_index(drop=True)


# Export data
tcr_data_all.to_excel('epitope_data_KRAS_TCRs.xlsx', index=False)













