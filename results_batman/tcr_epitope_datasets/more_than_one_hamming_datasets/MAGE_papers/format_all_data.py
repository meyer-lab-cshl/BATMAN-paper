import numpy as np
import pandas as pd

# Gather PMID: 23926201 data from two tables
data1 = pd.read_excel('5-197ra103_sm.xlsx',sheet_name='page 10')
data1 = data1.iloc[0:33,[10,13]]
data1.columns = ['activity','peptide']

data2 = pd.read_excel('5-197ra103_sm.xlsx',sheet_name='page 14')
data2 = data2.iloc[0:18,[10,11]]
data2.columns = ['activity','peptide']

data1 = pd.concat([data1,data2])

data1['activation'] = 1
data1.activation[data1.activity=='Yes'] = 2
data1.activation[data1.activity=='×'] = 0

data1 = data1.drop('activity',axis=1)

# Add data from PMID:35389803
data2 = pd.read_excel('science.abl5282_sm.xlsx',sheet_name='page 38',
                      skiprows=1)
data2 = pd.DataFrame(data2.iloc[:,2])
data2.columns = ['peptide']
data2['activation'] = 0 #initialize as NB

sb_set = ['EVDPIGHLY','ESDPIVAQY','ETDPVNHMV','EVDPIGHVY']
data2.activation[data2['peptide'].isin(sb_set)] = 2
#Add remaining SB data in Fig 6D
data2 = pd.concat([data2,
                   pd.DataFrame({'peptide':sb_set,
                                 'activation':[2,2,2,2]})])

# Add additional SBs from PMID: 38956325
sb_set = ['ETDPLTFNF','ETDPIEQVY']
all_data = pd.concat([data1,data2,
                   pd.DataFrame({'peptide':sb_set,
                                 'activation':[2,2]})])

all_data = all_data.drop_duplicates()
all_data.reset_index(drop=True, inplace=True)

# EVDPIRHYY has labels 0 and 1, make it WB
all_data=all_data.drop(49)

# Add constant columns
all_data['tcr'] = 'a3a'
all_data['index_peptide'] = 'EVDPIGHLY'


all_data.to_csv('a3a_multi_hamming_all.csv',index=False)


















