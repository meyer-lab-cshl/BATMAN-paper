##### Functions for different types of pooled inference
def pooled_inference_weights_only(TCR_index, # TCR index for each peptide
                                  peptide_feature, # features for each peptide
                                  peptide_binding_category,#binding category
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np
    
    # Infer number of TCRs and peptide length
    n_tcr = len(np.unique(TCR_index))
    peptide_length = np.shape(peptide_feature)[1]
    
    # Build Bayesian classifier
    with pm.Model() as peptide_classifier_model:        
        
        # TCR-specific parameters
        
        # positional weights, pooled over TCRs and positions
        weights = pm.Beta("weights",
                          alpha = pm.Gamma("alpha_w", mu=1.0, sigma=0.5),
                          beta = pm.Gamma("beta_w", mu=5.0, sigma=1.0),
                          shape=(n_tcr,peptide_length))        
        
        # Hyperprior parameters for TCR-specific intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1, shape=n_tcr) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)  
             
        
        # Full Predictor
        eta = intercepts[TCR_index] - pm.math.sum(
            peptide_feature*weights[TCR_index,:],axis=1)
        
        # Binomial Regression
        # Generate cutpoints
        cutpoints=pm.Normal("cutpoints", 
                            mu=0.0, sigma=2.0, shape=2,
         transform=pm.distributions.transforms.univariate_ordered,
         initval=[0.2,0.4])
        
        peptide_binding_category_obs = pm.OrderedLogistic(
                                        "peptide_binding_category_obs",
                                        eta=eta,cutpoints=cutpoints,
                                        observed=peptide_binding_category,
                                        compute_p=False)
        
    # Sampling with approximate posterior
    with peptide_classifier_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000))
         
    # Extract position-dependent weights of TCRs         
    inferred_weights=np.reshape(inferred_params.iloc[
        (n_tcr+3):(n_tcr+3+n_tcr*peptide_length),0].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
    
    return inferred_weights
####################################################################
def pooled_inference(TCR_index, # TCR index for each peptide
                                  peptide_feature, # features for each peptide
                                  peptide_binding_category,#binding category
                                  aa_change_index,#indexing which AA to which
                                  mode, # symmetric (default) or full
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np
    
    # Infer number of TCRs and peptide length
    n_tcr = len(np.unique(TCR_index))
    peptide_length = np.shape(peptide_feature)[1]
    
    # Renumber with unique aa_change_index
    _,unique_indices, aa_change = np.unique(aa_change_index, 
                                   return_index=True,
                                   return_inverse=True);
    
    # Create AA change matrix with peptide features
    aa_change_matrix = np.zeros((np.shape(peptide_feature)[0],
                                 len(unique_indices)))
    
    for peptide in np.arange(0,np.shape(peptide_feature)[0],1):
        aa_change_matrix[peptide,aa_change[peptide]] = sum(
            peptide_feature[peptide,]) # since there is only one non-zero element
        
    
    # Binarize peptide feature matrix
    peptide_feature[peptide_feature!=0]=1
    
    # Build Bayesian classifier
    with pm.Model() as peptide_classifier_model: 
        
        # TCR-indepedent common amino acid distance matrix flattened        
        aa_distance = pm.Normal("aa_distance",
                                mu=pm.Normal("mu", mu=0, sigma=0.5),
                                sigma=pm.Exponential("sigma",lam=1),
                                shape=(len(unique_indices),1))
        
        # TCR-specific parameters
        
        # positional weights, pooled over TCRs and positions
        weights = pm.Beta("weights",
                          alpha = pm.Gamma("alpha_w", mu=1.0, sigma=0.5),
                          beta = pm.Gamma("beta_w", mu=5.0, sigma=1.0),
                          shape=(n_tcr,peptide_length))        
        
        # Hyperprior parameters for TCR-specific intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1, shape=n_tcr) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)  
             
        
        # Full Predictor
        eta = intercepts[TCR_index] - pm.math.sum(
            peptide_feature*weights[TCR_index,:],axis=1)*(pm.math.dot(
            aa_change_matrix,(1+aa_distance)[:,0]))
                    
        
        # Binomial Regression
        # Generate cutpoints
        cutpoints=pm.Normal("cutpoints", 
                            mu=0.0, sigma=2.0, shape=2,
         transform=pm.distributions.transforms.univariate_ordered,
         initval=[0.2,0.4])
        
        peptide_binding_category_obs = pm.OrderedLogistic(
                                        "peptide_binding_category_obs",
                                        eta=eta,cutpoints=cutpoints,
                                        observed=peptide_binding_category,
                                        compute_p=False)
        
    # Sampling with approximate posterior
    with peptide_classifier_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000))
         
    # Extract position-dependent weights of TCRs         
    inferred_weights=np.reshape(inferred_params.iloc[
        (n_tcr+len(unique_indices)+5):
        (n_tcr+len(unique_indices)+5+n_tcr*peptide_length),0].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
        
    # Extract inferred AA distance factors 
    inferred_aa_distance = inferred_params.iloc[
                                        1:(1+len(unique_indices)),0].to_numpy()
    aa_indices = aa_change_index[unique_indices]
    
    # Reconstruct full AA factor matrix (to be multiplied with BLOSUM100)
    # Initialize with mean value
    inferred_aa_matrix = np.zeros((20,20))+(1+inferred_params.iloc[0,0])
    
    if mode!='full': #Default
        for aa1 in np.arange(0,19,1):
            for aa2 in np.arange(0,19,1):
               if 20*min(aa1,aa2)+max(aa1,aa2) in aa_indices:
                   inferred_aa_matrix[aa1,aa2] = 1+inferred_aa_distance[
                       np.where(aa_indices==20*min(aa1,aa2)+max(aa1,aa2))[0][0]]                   
                   # Takes care of AA reindexing done at the beginning
      
    if mode=='full':
        for aa1 in np.arange(0,19,1):
            for aa2 in np.arange(0,19,1):
               if 20*aa1+aa2 in aa_indices:
                   inferred_aa_matrix[aa1,aa2] = 1+inferred_aa_distance[
                       np.where(aa_indices==20*aa1+aa2)[0][0]]                   
                   # Takes care of AA reindexing done at the beginning
 
    return inferred_weights,inferred_aa_matrix
####################################################################


    



#################### Actual Script ########################

import numpy as np
import pandas as pd

# Load mutant data
epitope_data = pd.read_excel('../data/training_folds.xlsx');

# Filter data based on a criteria (peptide length or index peptide)

# Remove rows with NaN in aa change index
epitope_data = epitope_data[epitope_data.aa_change_symm.isnull().values==0]
epitope_data = epitope_data[epitope_data.aa_change_full.isnull().values==0]

# Select 9-mers (or 8-mers) only (compulsory)
peptide_length = 9

epitope_data = epitope_data[np.char.str_len(
    epitope_data.peptide_list.astype(str).tolist())==peptide_length]


# filter based on index peptide (optional)
index_peptide_selected ='TPQDLNTML'
epitope_data = epitope_data[epitope_data.index_peptide==index_peptide_selected]


# Renumber TCR index for the selected subset
unique_array, unique_indices = np.unique(
                                 epitope_data.TCR_index, return_inverse=True);
epitope_data.TCR_index = unique_indices.tolist()

# Store unique TCR names in the order they appear
TCR_names = pd.unique(epitope_data.TCR_name)

#AA names
amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
'''
################## Within TCR ##############################
# Loop over 5 folds and store inferred weights

for fold in np.arange(0,5,1):
    
    epitope_data_train = epitope_data[epitope_data.training_folds!=fold]
    
    # Build Bayesian Inference and Sample
    # Weights and symmetric/full AA matrix
    # Size: TCR-by-position
    inferred_weights,inferred_aa_matrix = pooled_inference(
                                epitope_data_train.TCR_index.to_numpy(),
                                epitope_data_train.iloc[0:len(epitope_data),
                                              7:(7+peptide_length)].to_numpy(),
                        epitope_data_train.peptide_binding_category.to_numpy(),
                        epitope_data_train.aa_change_full.to_numpy(),
                        mode='full',seed=100, steps=100000) 
    
    
    # Weights only
    # Size: TCR-by-position
    inferred_weights = pooled_inference_weights_only(
                                epitope_data_train.TCR_index.to_numpy(),
                                epitope_data_train.iloc[0:len(epitope_data),
                                              7:(7+peptide_length)].to_numpy(),
                        epitope_data_train.peptide_binding_category.to_numpy(),
                                  seed=100, steps=100000) 
    
    
    # If you run into NaN, try a new seed
    # Add TCR and AA name
    inferred_weights = pd.DataFrame(data=inferred_weights, 
                                    index=TCR_names)
    inferred_aa_matrix = pd.DataFrame(data=inferred_aa_matrix, 
                                    index=amino_acid_list,
                                    columns=amino_acid_list)
    # Save data
    inferred_weights.to_csv(''.join([
        '../data/within_tcr_pooled/full_aa_matrix/weights_all_',
        str(peptide_length),'_mer_fold_',str(fold),'.csv'])) 
    
    inferred_aa_matrix.to_csv(''.join([
        '../data/within_tcr_pooled/full_aa_matrix/aa_matrix_all_',
        str(peptide_length),'_mer_fold_',str(fold),'.csv']))
#######################################################################    
'''   
################### LOO TCRs ##########################################
# Loop over left out TCRs

for TCR in TCR_names:
    print(TCR)
    epitope_data_train = epitope_data[epitope_data.TCR_name!=TCR]#LOO TCR
    
    # Discard TCR indices, since we infer a common weight profile
    epitope_data_train.iloc[0:len(epitope_data_train),5] = 0
    
    
    # Build Bayesian Inference and Sample
    # Weights and symmetric AA matrix
    inferred_weights,inferred_aa_matrix = pooled_inference(
                                epitope_data_train.TCR_index.to_numpy(),
                                epitope_data_train.iloc[0:len(epitope_data),
                                              7:(7+peptide_length)].to_numpy(),
                        epitope_data_train.peptide_binding_category.to_numpy(),
                        epitope_data_train.aa_change_full.to_numpy(),
                        mode='symm',seed=100, steps=50000)
    
    inferred_weights = pd.DataFrame(data=inferred_weights)
    
    inferred_aa_matrix = pd.DataFrame(data=inferred_aa_matrix, 
                                    index=amino_acid_list,
                                    columns=amino_acid_list)
    
    inferred_weights.to_csv(''.join([
        '../data/loo_tcr/symm_aa_matrix/weights_',
        index_peptide_selected,'_',TCR,'.csv'])) 
    
    inferred_aa_matrix.to_csv(''.join([
        '../data/loo_tcr/symm_aa_matrix/aa_matrix_',
        index_peptide_selected,'_',TCR,'.csv']))
    
    
    # Weights and symmetric AA matrix
    inferred_weights,inferred_aa_matrix = pooled_inference(
                                epitope_data_train.TCR_index.to_numpy(),
                                epitope_data_train.iloc[0:len(epitope_data),
                                              7:(7+peptide_length)].to_numpy(),
                        epitope_data_train.peptide_binding_category.to_numpy(),
                        epitope_data_train.aa_change_full.to_numpy(),
                        mode='full',seed=100, steps=50000) 
    
    inferred_weights = pd.DataFrame(data=inferred_weights)
    
    inferred_aa_matrix = pd.DataFrame(data=inferred_aa_matrix, 
                                    index=amino_acid_list,
                                    columns=amino_acid_list)
    
    inferred_weights.to_csv(''.join([
        '../data/loo_tcr/full_aa_matrix/weights_',
        index_peptide_selected,'_',TCR,'.csv'])) 
    
    inferred_aa_matrix.to_csv(''.join([
        '../data/loo_tcr/full_aa_matrix/aa_matrix_',
        index_peptide_selected,'_',TCR,'.csv']))
    
    
    # Weights only
    # Size: TCR-by-position
    inferred_weights = pooled_inference_weights_only(
                                epitope_data_train.TCR_index.to_numpy(),
                                epitope_data_train.iloc[0:len(epitope_data),
                                              7:(7+peptide_length)].to_numpy(),
                        epitope_data_train.peptide_binding_category.to_numpy(),
                                  seed=100, steps=50000) 
    
    inferred_weights = pd.DataFrame(data=inferred_weights)
       
    inferred_weights.to_csv(''.join([
        '../data/loo_tcr/weights_only/weights_',
        index_peptide_selected,'_',TCR,'.csv'])) 
    
    






