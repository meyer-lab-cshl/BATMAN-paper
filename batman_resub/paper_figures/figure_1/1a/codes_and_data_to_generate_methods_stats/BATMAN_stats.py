# Reads BATMAN training data and outputs the number of NB and WB+SB
import numpy as np
import pandas as pd
import sys
# Add the BATMAN function directory to the system path
sys.path.append('../../../../batman_functions')
from batman_functions import *

# Load full mutant epitope activity data for MHCI and MHCII
mutant_activity_data_mhci = pd.read_excel(
    ''.join(['../../../../tcr_epitope_datasets/mutational_scan_datasets/',
             'TCR_pMHCI_mutational_scan_database.xlsx']))

mutant_activity_data_mhcii = pd.read_excel(
    ''.join(['../../../../tcr_epitope_datasets/mutational_scan_datasets/',
             'TCR_pMHCII_mutational_scan_database.xlsx']))

mutant_activity_data = pd.concat([mutant_activity_data_mhci,
                                  mutant_activity_data_mhcii]).reset_index(drop=True)

# Normalize peptide activity and add peptide activation labels
mutant_activity_data = assign_activation_category(mutant_activity_data)
 
# Get statistics   
print(len(pd.unique(mutant_activity_data['tcr'])), 
      len(pd.unique(mutant_activity_data['peptide'])),
      len(mutant_activity_data[mutant_activity_data.activation==0]),
      len(mutant_activity_data[mutant_activity_data.activation>0]))