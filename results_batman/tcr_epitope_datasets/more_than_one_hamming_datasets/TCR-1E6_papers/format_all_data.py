import numpy as np
import pandas as pd

# Gather PMID: 22102287 data 
all_data = pd.read_excel('peptide_response_1E6_data.xlsx')
all_data.columns = ['Number','peptide','activity']

# Change Null to 0
all_data.activity[all_data.activity=='Null'] = '0.0'
all_data['activity'] = all_data['activity'].astype(float)

# Get index peptide activity
index_activity = np.mean(all_data.activity[all_data.peptide=='ALWGPDPAAA'])

# MinMax normalize
min_activity = np.sort(all_data.activity)[1] #Ignore 0
all_data.activity = (all_data.activity - min_activity)/(np.max(all_data.activity)-min_activity)


#Assign TCR activation
all_data['activation'] = 1 #initialize
all_data.loc[all_data.activity<0.1,'activation']=0 #Non-binder
all_data.loc[all_data.activity>=0.5,'activation']=2 #strong-binder

#Take out an instance where the index is marked WB
all_data = all_data.drop(all_data.index[58])

# Add constant columns
all_data['tcr'] = 'TCR-1E6'
all_data['index_peptide'] = 'ALWGPDPAAA'

all_data = all_data.drop('activity',axis=1)
all_data = all_data.drop('Number',axis=1)

all_data = all_data.drop_duplicates()
all_data.reset_index(drop=True, inplace=True)

all_data.to_csv('TCR-1E6_multi_hamming_all.csv',index=False)


















