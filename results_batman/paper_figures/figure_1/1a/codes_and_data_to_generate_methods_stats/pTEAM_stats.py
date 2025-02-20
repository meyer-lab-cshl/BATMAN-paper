# Reads pTEAM training data and outputs the number of non-binder and binder
import pandas as pd
import numpy as np

#Sheet 1
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pteam_S1.xlsx',sheet_name='Normalized data',
                                     skiprows=1,index_col=0)
mutant_activity = np.array(mutant_activity_data)

# initialize number of binder and non-binder with stats for the Ova-repertoire
non_binder = len(mutant_activity[mutant_activity<46.9])
binder = len(mutant_activity[mutant_activity>=46.9])

#Sheet 2
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pteam_S2.xlsx',sheet_name='Normalized data',
                                     skiprows=1,index_col=0)
mutant_activity = np.array(mutant_activity_data)

# add number of binder and non-binder
non_binder = non_binder + len(mutant_activity[mutant_activity<46.9])
binder = binder + len(mutant_activity[mutant_activity>=46.9])

#Sheet 3
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pteam_S3.xlsx',sheet_name='Normalized by PC');
mutant_activity = np.array(mutant_activity_data)[:,4:11]

# add number of binder and non-binder
non_binder = non_binder + len(mutant_activity[mutant_activity< 66.09])
binder = binder + len(mutant_activity[mutant_activity>=  66.09])

#Sheet 4
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pteam_S4.xlsx',sheet_name='Mean',
                                     index_col=0)
mutant_activity = np.array(mutant_activity_data)

# add number of binder and non-binder
non_binder = non_binder + len(mutant_activity[mutant_activity< 40.0])
binder = binder + len(mutant_activity[mutant_activity>=  40.0])

    
print(non_binder,binder)