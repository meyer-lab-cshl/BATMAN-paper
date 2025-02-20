# Reads TCRAI training data and outputs the number of non-binder and binder
import pandas as pd
import numpy as np

# Load full mutant epitope activity data
mutant_activity_data1 = pd.read_csv('public_TCRs.csv').astype('str');
mutant_activity_data2 = pd.read_csv('CNN-prediction-with-REGN-pilot-version-2.csv').astype('str');

mutant_activity_data = pd.concat([mutant_activity_data1, mutant_activity_data2]).drop_duplicates()


# initialize number of binder and non-binder (discard random ones)
non_binder = len(mutant_activity_data[mutant_activity_data.binds=='0'])
binder = len(mutant_activity_data[mutant_activity_data.binds=='1'])

# unique TCRs and epitopes
TCR = len(np.unique(mutant_activity_data[['TRB_v_gene',
                                'TRB_j_gene',
                                'TRA_v_gene',
                                'TRA_j_gene',
                                'TRB_cdr3',
                                'TRA_cdr3']].apply(''.join, axis=1)))

epitope = len(np.unique(mutant_activity_data.id))
  
print(TCR,non_binder,binder,epitope)

