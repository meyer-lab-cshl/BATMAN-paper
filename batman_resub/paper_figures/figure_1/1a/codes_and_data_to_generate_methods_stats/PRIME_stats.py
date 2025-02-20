# Reads PRIME 2.0 training data and outputs the number of non-binder and binder
import pandas as pd

#Sheet 1
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('1-s2.0-S2405471222004707-mmc6.xlsx',
                                     sheet_name='TableS4');


# initialize number of binder and non-binder (discard random ones)
non_binder = len(mutant_activity_data[
                          (mutant_activity_data.Random==0) &
                          (mutant_activity_data.Immunogenicity==0)])
binder = len(mutant_activity_data[
                          (mutant_activity_data.Random==0) &
                          (mutant_activity_data.Immunogenicity==1)])
   
print(non_binder,binder)