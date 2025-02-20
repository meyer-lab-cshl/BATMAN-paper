import numpy as np
import pandas as pd

# Gather PMID: 23926201 data from two tables
all_data = pd.read_excel('aav0860_table_s5.xlsx',sheet_name='Figure 6d',
                         skiprows=3)
all_data = all_data.iloc[:,[0,4]]
all_data.columns = ['peptide','activity']

# Normalize with max of 1 Hamming
all_data.activity = all_data.activity/249.702057613169

#Assign TCR activation
all_data['activation'] = 1 #initialize
all_data.loc[all_data.activity<0.1,'activation']=0 #Non-binder
all_data.loc[all_data.activity>=0.5,'activation']=2 #strong-binder

# Add constant columns
all_data['tcr'] = '868Z11'
all_data['index_peptide'] = 'SLYNTVATL'

all_data = all_data.drop('activity',axis=1)

all_data.to_csv('868Z11_multi_hamming_all.csv',index=False)


















