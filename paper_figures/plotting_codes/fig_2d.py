import numpy as np
import pandas as pd
from original_batman_functions import train, peptide2index
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns

# Load epitope list
epitope_data = pd.read_excel('../data/TCR_epitope_database.xlsx');

# Filter out TCRs without any available sequence
TCRs_included = ['18A2','A23','T1','TCR3','TCR6','TCR7','FLT3DY','868Z11','A6',
                 'NYE-S1','TCR2-T']
epitope_data = epitope_data[epitope_data.tcr_name.isin(TCRs_included)==True]

# Filter TCRs binding to 9-mers
peptide_length = 9
sequence_length=np.char.str_len(np.array(epitope_data.peptide).astype(str))
epitope_data = epitope_data[sequence_length==peptide_length]

# Get unique TCRs in data
TCR_names, TCR_index = np.unique(epitope_data['tcr_name'], return_inverse=True)

# Normalize peptide activities and indicate peptide activity level
for TCR in TCR_names:
    activity = epitope_data.loc[epitope_data.tcr_name==TCR,'peptide_activity']
    max_activity = max(activity)
    epitope_data.loc[epitope_data.tcr_name==TCR,
                     'peptide_activity'] = activity/max_activity
    
peptide_activity = epitope_data.peptide_activity
peptide_activity_level = np.zeros((len(peptide_activity),1))+1 #initialize

peptide_activity_level[peptide_activity<0.1]=0 #Non-binder
peptide_activity_level[peptide_activity>=0.5]=2 #strong-binder

epitope_data['peptide_activity_level'] = peptide_activity_level

# Create and save a file for training BATMAN
batman_train_data = pd.DataFrame(
    data={'tcr':epitope_data['tcr_name'],
          'peptide':epitope_data['peptide'],
          'index':epitope_data['index_peptide'],
          'activation':epitope_data['peptide_activity_level'].astype(int)})
 
batman_train_data.to_csv('batman_train_data.csv', index = False)


# Train BATMAN to infer TCR-specific weights


_,inferred_matrix_full = train('batman_train_data.csv','full',
                           'blosum100', steps = 80000, seed = 100)

# Plotting
# Order of substituted amino acids for good clusering
AA_list = np.array(['I','V','L','F','C','M','A','W','G','T','S','Y','P','H',
                    'N','D','Q','E','K','R'])

AA_distance_matrix_full = inferred_matrix_full.loc[AA_list,AA_list]

AA_matrix = pd.read_csv('../data/blosum100.csv',index_col=0)
blosum = AA_matrix.loc[AA_list,AA_list]

# plot heatmaps
sns.set_theme(style="whitegrid")
sns.set(font_scale=2.3)
plt.figure(figsize=(15, 15))
#plt.suptitle('Inferred AA substitution distance matrix from BATMAN')

sns.heatmap((AA_distance_matrix_full/2.984)/blosum,cmap="RdBu",
            # divided the pooled mean of AA matrix So that ratio with blosum=1 
            yticklabels=AA_list,
            xticklabels=AA_list,center=1,vmin=0,vmax=1.6,
            square=1)
plt.yticks(rotation=0);
plt.xlabel('To AA');
plt.ylabel('From AA');
#plt.title('Inferred symmetric matrix');

plt.show()
plt.savefig('aa_matrix_inferred.pdf', format='pdf')

