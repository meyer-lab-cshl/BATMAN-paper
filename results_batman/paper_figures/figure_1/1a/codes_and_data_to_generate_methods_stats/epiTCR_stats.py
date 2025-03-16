# Reads epiTCR training data and outputs the number of non-binder and binder
import pandas as pd
import numpy as np

# Load full mutant epitope activity data
mutant_activity_data = pd.read_csv('full-training-with-categories.csv');


# initialize number of binder and non-binder (discard random ones)
non_binder = len(mutant_activity_data[mutant_activity_data.binder==0])
binder = len(mutant_activity_data[mutant_activity_data.binder==1])

# unique TCRs and epitopes
TCR = len(np.unique(mutant_activity_data.CDR3b))
epitope = len(np.unique(mutant_activity_data.epitope))
  
print(TCR,non_binder,binder,epitope)

