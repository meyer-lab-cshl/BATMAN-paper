import numpy as np
import pandas as pd

# Gather PMID: 36424374 data from two tables
data1 = pd.read_excel('multiple_mutation_IL.jpg.xlsx',sheet_name='page-1_table-1',
                         skiprows=0)
data1 = data1.iloc[:,[1,3]]
data1.columns = ['peptide','activity']

data2 = pd.read_excel('multiple_mutation_IL.jpg.xlsx',sheet_name='page-1_table-2',
                         skiprows=0)
data2 = data2.iloc[:,[1,3]]
data2.columns = ['peptide','activity']

all_data = pd.concat([data1,data2])

# Normalize with max of 1 Hamming, and by scaling of index peptide activity
all_data.activity = all_data.activity*(3592.14/100)*(1/4196.05)

#Assign TCR activation
all_data['activation'] = 1 #initialize
all_data.loc[all_data.activity<0.1,'activation']=0 #Non-binder
all_data.loc[all_data.activity>=0.5,'activation']=2 #strong-binder

# Add constant columns
all_data['tcr'] = 'TIL1383I'
all_data['index_peptide'] = 'YMDGTMSQV'

all_data = all_data.drop('activity',axis=1)

all_data.to_csv('TIL1383I_multi_hamming_all.csv',index=False)


















