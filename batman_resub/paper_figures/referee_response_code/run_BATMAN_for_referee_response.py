import numpy as np
import pandas as pd
import arviz as az
import pymc as pm
import matplotlib.pyplot as plt


'''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)

# Subset 9-mer-binding benchmark TCRs
TCR_list = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

'''Generate mutant peptide features'''
#AA names
amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']

# Load AA matrix
aa_distance_matrix = pd.read_csv('blosum100.csv',index_col=0)

# save row and column AA orders for row to column AA substitutions
from_AA = np.array(aa_distance_matrix.index).astype(str)
to_AA = np.array(aa_distance_matrix.columns).astype(str)

# Convert to np matrix for easy indexing
aa_distance_matrix = aa_distance_matrix.to_numpy()

peptide_length = 9

index_peptide_repeated = np.array(list("".join(tcr_pmhc_data['index'])))
mutant_peptides_joined = np.array(list("".join(tcr_pmhc_data['peptide'])))

# Compare mutant seqs to index peptide seq to find mutation locations
is_mutation_location = np.compare_chararrays(list(index_peptide_repeated),
                        list(mutant_peptides_joined),
                        "!=",'True')

# Locations of WT AAs in from_AA array
from_AA_index = np.nonzero(index_peptide_repeated[is_mutation_location,None] ==
                           from_AA)[1]
# Locations of mutated AAs in to_AA array
to_AA_index = np.nonzero(mutant_peptides_joined[is_mutation_location,None] ==
                           to_AA)[1]

# Collect distance matrix elements corresponding to mismatches
aa_distance = np.zeros(len(index_peptide_repeated)) #initialize
aa_distance[is_mutation_location] = aa_distance_matrix[from_AA_index,
                                                       to_AA_index]

# Reshape AA distance array to dim #mutants-by-peptide_length
peptide_feature = np.reshape(aa_distance, 
                         (len(tcr_pmhc_data),peptide_length))


'''featurize TCR-pMHC data'''        
# Assign indices to unique TCRs in data
TCR_names, TCR_index = np.unique(tcr_pmhc_data['tcr'], return_inverse=True)
n_tcr = len(np.unique(TCR_index))

# Peptide activation categories
peptide_activation_category = tcr_pmhc_data['activation'].to_numpy()

# store MHC binding data as binary feature vector
bindlevel = tcr_pmhc_data['BindLevel']
peptide_bindlevel = np.zeros((len(tcr_pmhc_data),3))
peptide_bindlevel[bindlevel=='NB',0] = 1
peptide_bindlevel[bindlevel=='WB',1] = 1
peptide_bindlevel[bindlevel=='SB',2] = 1

#full AA matrix, indicate a single number to index AA substitution, 
#and reshape flattened array
# Indices of index and mutant AAs in AA list (flattened)

index_aa = np.where(np.array(list(''.join(index_peptide_repeated)))[:, None] == 
                    np.array(amino_acid_list)[None, :])[1]
peptide_aa = np.where(np.array(list(''.join(mutant_peptides_joined)))[:, None] == 
                      np.array(amino_acid_list)[None, :])[1]
aa_subs = (20*index_aa + peptide_aa)
aa_subs[index_aa==peptide_aa]=0 #assign 0 to positions with unchanged AA 
aa_change_index = aa_subs.reshape(np.shape(peptide_feature))

# Number of TCR activation levels
n_level = 3

# Renumber AA change index matrix with unique AA change indices
unique_indices,_, aa_change_index_unique = np.unique(aa_change_index, 
                               return_index=True,
                               return_inverse=True);
aa_change_index_unique = aa_change_index_unique.reshape(
                                                 np.shape(peptide_feature))


'''Build Bayesian classifier'''
with pm.Model() as peptide_classifier_model: 
    
    # TCR-indepedent common amino acid distance matrix multiplier flattened        
    aa_distance_multiplier = pm.math.concatenate([[0],
                                      pm.Normal("BLOSUM_ratio",
                            mu=pm.Normal("mu", mu=0, sigma=0.5),
                            sigma=pm.Exponential("sigma",lam=1),
                            shape=len(unique_indices)-1)], 
                                      axis=0)
    #0 at beginning put for No AA substituion
    
    if np.sum(peptide_bindlevel)!=0: #If MHC binding info is present
        # parameters for incorporating MHC-binding information            
        
        mhc_effect_add = pm.Beta("mhc_effect_add", # MHC feature effect, pooled across TCRs
                           alpha = pm.Gamma("alpha_add", alpha=2, beta=2),
                           beta = pm.Gamma("beta_add", alpha=2, beta=2),
                           shape = (n_tcr,3))
                    
    # positional weights, pooled over TCRs and positions
    weights = pm.Beta("positional_weights",
                      alpha = pm.Gamma("alpha_w", alpha=2, beta=2),
                      beta = pm.Gamma("beta_w", alpha=2, beta=2),
                      shape=(n_tcr,peptide_length))        
    
    # Hyperprior parameters for TCR-specific intercept        
    mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
    sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
    normal = pm.Normal("normal", mu=0, sigma=1, shape=n_tcr) 
    
    intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)
    
    #Full predictor
    # baseline
    eta = - pm.math.sum(
        weights[TCR_index,:]*peptide_feature*
        (1 + aa_distance_multiplier[aa_change_index_unique]),
        axis=1) #positional weight * D *(1+multiplier) summed over positions 
    
    if np.sum(peptide_bindlevel)!=0: # pMHC BindLevel info available
        eta = eta - pm.math.sum(peptide_bindlevel*mhc_effect_add[TCR_index,:],
                                                    axis=1)
   
       
    eta = eta + intercepts[TCR_index]
    # Binomial Regression
    # Generate cutpoints
    cutpoints=pm.Normal("cutpoints", 
                        mu=0.0, sigma=2.0, shape=[n_level-1],
     transform=pm.distributions.transforms.univariate_ordered,
     initval=np.linspace(0.1,0.5,n_level-1))
    
    peptide_activation_category_obs = pm.OrderedLogistic(
                                    "peptide_activation_category_obs",
                                    eta=eta,cutpoints=cutpoints,
                                    observed=peptide_activation_category,
                                    compute_p=False)
    
'''ADVI sampling and saving trace'''
with peptide_classifier_model:
    advi = pm.ADVI()
    
# Make tracker callback for stats
tracker = pm.callbacks.Tracker(
mean=advi.approx.mean.eval,  # callable that returns mean
std=advi.approx.std.eval)  # callable that returns std

# ADVI fit to posterior
approx = advi.fit(80000, callbacks=[tracker],
                  progressbar=True)

def inv_logit(p):
    return np.exp(p) / (1 + np.exp(p))

# Save inferred parameter traces (returned as logodds which we logit inverse)
traces_mean=inv_logit(np.array(tracker["mean"]))
traces_std=np.array(tracker["std"])
trace_elbo = approx.hist #ELBO trace

'''Plotting'''

#Weights index for all 14 TCRs binding 9 mers in benchmark
plt.rcParams.update({'font.size': 28})

plt.subplot(2,1,1)
plt.plot(np.arange(len(traces_mean[:,443:569])),
        traces_mean[:,443:570], '-k')
plt.ylim(0,1)
plt.xlabel('ADVI steps')
plt.ylabel('Mean of inferred positional weight')
plt.title('Inferred positional weight traces')

plt.subplot(2,1,2)
plt.plot(np.arange(len(traces_mean[:,443:569])),
        trace_elbo, '-k')
#plt.ylim(0,1)
plt.xlabel('ADVI steps')
plt.ylabel('ELBO loss')
plt.title('ELBO trace')


        


    






