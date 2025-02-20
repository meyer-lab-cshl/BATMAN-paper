import numpy as np
import pandas as pd

all_data = pd.DataFrame(columns=['peptide',
                                 'peptide_activity','activation'])


# Gather all logFC data
for p1 in np.arange(4,8,1):
    for p2 in np.arange(p1+1,9,1):
        new_data = pd.read_excel(''.join(['data_from_heatmap_',str(p1),
                                        str(p2),'.xlsx']))
            
        #Assign TCR activation
        new_data['activation'] = 1 #initialize
        new_data.loc[new_data.peptide_activity<0.1,'activation']=0 #Non-binder
        new_data.loc[new_data.peptide_activity>=0.5,'activation']=2 #strong-binder
        
        # append data
        all_data = pd.concat([all_data,new_data])
        
# Gather flow data
flow_data = pd.read_excel('flow_data.xlsx')
#Assign TCR activation
flow_data['activation'] = 1 #initialize
flow_data.loc[flow_data.peptide_activity<=0.03,'activation']=0 #Non-binder
flow_data.loc[flow_data.peptide_activity>=1.25,'activation']=2 #strong-binder

# append data
all_data = pd.concat([all_data,flow_data])
        
all_data.reset_index(drop=True, inplace=True)

# Add constant columns
all_data['tcr'] = 'EWW-TCR'
all_data['index_peptide'] = 'EWWRSGGFSF'

all_data = all_data.drop('peptide_activity',axis=1)

all_data.to_csv('EWW-TCR_multi_hamming_all.csv',index=False)


















