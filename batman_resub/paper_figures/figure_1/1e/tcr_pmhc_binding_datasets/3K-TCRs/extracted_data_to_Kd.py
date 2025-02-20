import numpy as np
import pandas as pd

# Read raw data
raw_data = pd.read_excel('data_extracted_from_bar_plots.xlsx')

# TCR WT Kd and DG list from Table 1
wt_data = pd.DataFrame(data={'tcr':['B3K506', 'B3K508', 'YAe5-62.8',
                             'TCR75-1', 'TCR2W1S12-20.4'],
                      'Kd':[7E-6,2.9E-5,8E-6,8E-6,1.3E-5],
                      'DG':[-7,-6.2,-6.9,-7.0,-6.6]})

# Initialize data record with index peptide
index_peptide = "FEAQKAKANKAVD"

all_data = wt_data[['tcr','Kd']]
all_data['index_peptide'] = index_peptide
all_data['peptide'] = index_peptide
all_data['Kd_is_greater'] = 0 #Kd is equal to record, and not >

wt_data=wt_data.set_index('tcr')

# Record the extracted data
for i in np.arange(len(raw_data)):
    peptide_data = raw_data.iloc[i,:].copy()
    
    # Construct mutant peptide sequence
    peptide=list(index_peptide)
    peptide[peptide_data.position+1] = peptide_data.mutant_aa
    peptide = ''.join(peptide)
    
    if np.isnan(peptide_data.Kd_Table_1): #Kd not found in Table 1
        if np.isnan(peptide_data.DDG_extracted): #Extracted DDG not available
            # Calculate Kd from DDG found from bar plot
            # DDG from bar height
            DDG = (peptide_data.ymax-peptide_data.ymin)/(
                                          peptide_data.y_unit-peptide_data.ymin)
            # Kd from DDG and WT DG
            Kd = np.exp((DDG + wt_data.loc[peptide_data.tcr,'DG'])/0.592)
            Kd_is_greater = abs(np.isnan(peptide_data.DDG_is_greater).astype(int)-1)
        else: #DDG data extracted is available
            DDG = peptide_data.DDG_extracted
            # Kd from DDG and WT DG
            Kd = np.exp((DDG + wt_data.loc[peptide_data.tcr,'DG'])/0.592)
            Kd_is_greater = 1 #Always 1 if DDG is directly recorded
        
    else: #Kd found in Table 1 and recorded
        Kd = peptide_data.Kd_Table_1
        Kd_is_greater = abs(np.isnan(peptide_data.Kd_is_greater).astype(int)-1)
     
    # Add new data
    data_add = pd.DataFrame({'tcr':[peptide_data.tcr],
                             'index_peptide':[index_peptide],
                             'peptide':[peptide],
                             'Kd':[Kd],
                             'Kd_is_greater':[Kd_is_greater]})
    all_data = pd.concat([all_data,data_add])

# Reorder columns and save data
all_data = all_data[['tcr','index_peptide','peptide','Kd','Kd_is_greater']]    
all_data.to_excel('kd_3K_TCRs.xlsx', index=False)   
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
