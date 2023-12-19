# Reads pTEAM training data and outputs the number of non-binder and binder
import pandas as pd
import numpy as np

#Sheet 1
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('media-1.xlsx',sheet_name='Normalized data');
mutant_activity = np.array(mutant_activity_data)[1:154,1:16]

# initialize number of binder and non-binder
non_binder = len(mutant_activity[mutant_activity<=46.9])
binder = len(mutant_activity[mutant_activity>46.9])

#Sheet 2
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('media-2.xlsx',sheet_name='Normalized data');
mutant_activity = np.array(mutant_activity_data)[1:154,1:22]

# initialize number of binder and non-binder
non_binder = non_binder + len(mutant_activity[mutant_activity<=46.9])
binder = binder + len(mutant_activity[mutant_activity>46.9])

#Sheet 3
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('media-3.xlsx',sheet_name='Normalized by PC');
mutant_activity = np.array(mutant_activity_data)[0:134,4:11]

# initialize number of binder and non-binder
non_binder = non_binder + len(mutant_activity[mutant_activity<=46.9])
binder = binder + len(mutant_activity[mutant_activity>46.9])
    
print(non_binder,binder)