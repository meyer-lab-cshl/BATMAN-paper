# Reads pTEAM training data and outputs the number of non-binder and binder
import pandas as pd

#Sheet 1
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pcbi.1003266.s001.xlsx',sheet_name='Sheet1');


# initialize number of binder and non-binder
non_binder = len(mutant_activity_data[
                          mutant_activity_data.Immunogenicity == 'non-immunogenic'])
binder = len(mutant_activity_data[
                          mutant_activity_data.Immunogenicity == 'immunogenic'])

#Sheet 2
# Load full mutant epitope activity data
mutant_activity_data = pd.read_excel('pcbi.1003266.s002.xlsx',sheet_name='Sheet1');


# initialize number of binder and non-binder
non_binder = non_binder + len(mutant_activity_data[
                          mutant_activity_data.epitope_status == 'non-epitope'])
binder = binder + len(mutant_activity_data[
                          mutant_activity_data.epitope_status == 'epitope'])

   
print(non_binder,binder)