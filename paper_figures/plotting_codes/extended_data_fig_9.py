import numpy as np
import pandas as pd
from original_batman_functions import generate_mutant_features,peptide2index
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

# Load epitope list
epitope_data = pd.read_csv('../data/batman_train_data_with_tcr_names.csv');

# List of TCR names
TCR_list = ['18A2','TCR3','A23','T1','868Z11','TCR2-T','TCR7','A6',
            'FLT3DY','NYE-S1','TCR6']
peptide_length = 9
np.random.seed(100)

# Read data on AUC improvement by 1 AL step
auc_improve = pd.read_csv('../data/AL_OPT/auc_improvement_by_al.csv',
                          index_col=1).iloc[:,1:4]

auc_improve = auc_improve.loc[TCR_list,:]

max_d = np.zeros((1,11)) # array to store max peptide2index distance among mutants

plt.figure(figsize=(15,15))

for tcr in np.arange(11):
    
    loo_tcr = TCR_list[tcr]
    
    # LOO TCR data
    loo_data = epitope_data[epitope_data.tcr == loo_tcr].copy().reset_index(drop=True)
    
    # read LOO weights and matrix
    w = pd.read_csv(''.join(['../data/loo_distances/w_',
                                         loo_tcr,'.csv']),
                    index_col=0).to_numpy()
    #w=np.ones((1,9))
    
    aa_mat = pd.read_csv(''.join(['../data/loo_distances/aa_',
                                         loo_tcr,'.csv']),
                    index_col=0)
    
    # Get distances
    loo_distances = peptide2index(loo_data['index'].to_list(), 
                      loo_data['peptide'].to_list(), 
                      aa_mat, 
                      w)
    
    max_d[0,tcr] = np.median(loo_distances)
        
    ax = plt.subplot(6,2,tcr+1)
    
    # All possible mutants
    plt.plot(loo_distances[loo_data.activation==0],
             np.random.normal(0,0.1,len(loo_distances[loo_data.activation==0])),
             'bo',markersize=10,alpha=0.25)
    
    plt.plot(loo_distances[loo_data.activation==1],
             np.random.normal(1,0.1,len(loo_distances[loo_data.activation==1])),
             'go',markersize=10,alpha=0.25)
    
    plt.plot(loo_distances[loo_data.activation==2],
             np.random.normal(2,0.1,len(loo_distances[loo_data.activation==2])),
             'ro',markersize=10,alpha=0.25)
    
    # OPT data
    # Load top OPT peptides
    opt_data = pd.read_csv(''.join(['../data/AL_OPT/',loo_tcr,'_opt.csv'])).iloc[0:7,:];
    
    opt_full = loo_data[loo_data.peptide.isin(opt_data['peptide'].tolist())==True]
    
    opt_distances = peptide2index(opt_full['index'].to_list(), 
                                  opt_full['peptide'].tolist(), 
                                  aa_mat, 
                                  w)
     
    
    plt.plot(opt_distances[opt_full.activation==0],
             np.random.normal(0,0.01,len(opt_full[opt_full.activation==0])),
             'ko',markersize=10)
    
    plt.plot(opt_distances[opt_full.activation==1],
             np.random.normal(1,0.01,len(opt_full[opt_full.activation==1])),
             'ko',markersize=10)
    
    plt.plot(opt_distances[opt_full.activation==2],
             np.random.normal(2,0.01,len(opt_full[loo_data.activation==2])),
             'ko',markersize=10)
    
    plt.xlim([-0.35 , 6])
    plt.yticks([0,1,2],labels=["Non","Weak","Strong"])
    
    ax.set_title(loo_tcr)
    
    
ax = plt.subplot(6,2,12)

plt.scatter(max_d.T,auc_improve.iloc[:,2],s=60,
         c = np.array(['#088F8F','#808080','#5D3FD3','#E4D00A','#0096FF',
                       '#FAA0A9','#50C878','#b30000','#0047AB','#EC5800',
                       '#bc80bd']),alpha=0.9)
'''
plt.ylim([0.64 , 0.88])
plt.yticks([0.65,0.75,0.85])
plt.xticks([0.8,1.0,1.2])
'''
ax.set_aspect(2)



plt.tight_layout()
plt.savefig("../figures/extended_data_fig9/mutant_distances.pdf",
            format="pdf", bbox_inches="tight")




