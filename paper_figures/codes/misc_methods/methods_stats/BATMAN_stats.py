# Reads BATMAN training data and outputs the number of NB and WB+SB
import numpy as np
import pandas as pd

# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('BATMAN_data.xlsx');

# initialize number of binder and non-binder
binder = 0;
non_binder = 0;
# List TCR names available in the dataset for the selected index peptide
TCR_name_list = np.unique(mutant_activity_data['tcr_name'].values)

# binding category for mutant peptides for a given TCR 
for TCR_index in np.arange(0,len(TCR_name_list),1):
    TCR = TCR_name_list[TCR_index] # Name of TCR
    
    mutant_activities = mutant_activity_data[
        mutant_activity_data.tcr_name==TCR]['peptide_activity'].values
    
    # Normalize mutant activities to 1
    mutant_activities = mutant_activities/max(mutant_activities)
    non_binder = non_binder + len(mutant_activities[mutant_activities<0.1])
    binder = binder + len(mutant_activities[mutant_activities>=0.1])
    
print(non_binder,binder)