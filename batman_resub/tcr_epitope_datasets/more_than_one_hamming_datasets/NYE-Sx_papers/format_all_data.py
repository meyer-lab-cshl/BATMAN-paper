import numpy as np
import pandas as pd

data1 = pd.read_excel('NYE_SLL_S1 ScHLA.xlsx')
data2 = pd.read_excel('NYE_SLL_S2 ScHLA.xlsx')
data3 = pd.read_excel('NYE_SLL_S3 ScHLA.xlsx')
data= pd.concat([data1,data2,data3])

data = pd.DataFrame(pd.unique(data['Peptide']))
data.columns = ['Peptide']
data['mhc'] = 'HLA-A*02:01'
data=data[['mhc','Peptide']]
data.to_csv('all_data.csv')

