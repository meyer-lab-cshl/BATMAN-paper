import os
# Get the path of the current module
current_module_path = os.path.dirname(__file__)
# Define the path to your data directory
data_directory = current_module_path

# Access your data files using relative paths
#dataset1_path = os.path.join(data_directory, 'dataset1.csv')

##############################################################
'''
Function to generate features of a list of mutant peptides, using their
position-dependent distances from the index peptide,
based on a given AA distance matrix.

Index peptide can be a string (common for all mutants) or a list (potentially
different for all mutants) with length equal to that of mutant list. All input 
sequences must have equal length.

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

'''

def generate_mutant_features(index_peptide, # Single index peptide or list 
                 mutant_peptide_list, #List of mutant peptide sequences
                 aa_matrix): #Named or user-defined AA matrix   
    
    import pandas as pd
    import numpy as np
    import sys
    import os
    
    #AA names
    amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                     'S','T','V','W','Y']
    
    # Check if mutant peptide list is a list
    if (isinstance(mutant_peptide_list,list) == False):
        sys.exit("Mutant peptides must be supplied as a list")
        
    # Check if there is only one index or as many as mutants
    if ((isinstance(index_peptide, str)==False) and 
        (len(mutant_peptide_list)!=len(index_peptide))):
        sys.exit("Supply either a unique index peptide or one for each mutant")
     
    # Check if index peptides form a list
    if ((isinstance(index_peptide, str)==False) and 
        ((isinstance(index_peptide,list) == False))):
        sys.exit("More than one index peptide must be supplied as a list")
        
    
    # Check if all mutants and index peptides have the same length
    if isinstance(index_peptide, str): #One index peptide provided
        if len(np.unique(np.char.str_len([index_peptide]+mutant_peptide_list)))!=1:
            sys.exit("All index and mutant sequences must have equal length")
    else: #Index peptide list provided
        if len(np.unique(np.char.str_len(index_peptide+mutant_peptide_list)))!=1:
            sys.exit("All index and mutant sequences must have equal length")
    
    
    
    # Check if input sequences contain non-AA and wildcard characters
    if isinstance(index_peptide, str): #One index peptide provided
        if set(list("".join(([index_peptide] +
                             mutant_peptide_list)))).issubset(
                                 set(amino_acid_list)) == False:
           sys.exit("Discard sequences with non-AA characters, including X and *")
    else: #Index peptide list provided
        if set(list("".join((index_peptide +
                             mutant_peptide_list)))).issubset(
                                 set(amino_acid_list)) == False:
           sys.exit("Discard sequences with non-AA characters, including X and *")        
    
    
    
    # If AA matrix provided, check size is 20-by-20 and it has row and column names
    if isinstance(aa_matrix, str)==False: #Name of AA matrix not provided
        if ((isinstance(aa_matrix, pd.DataFrame) == False) or
            (aa_matrix.shape != (20,20)) or
            (set(aa_matrix.index) != set(amino_acid_list)) or
            (set(aa_matrix.columns) != set(amino_acid_list))):
            sys.exit("".join(["Custom AA matrix must be a 20-by-20 Pandas ",
                        "Dataframe object with AAs as row and column names."]))
        else: #If all goes well, load custom AA matrix
            aa_distance_matrix = aa_matrix
    
    # Load named AA distance matrix if the name is provided
    if isinstance(aa_matrix, str): #Name of AA matrix provided
        # Check if input AA distance matrix name exists, if it is provided
        filename = "".join([aa_matrix,'.csv'])
        if os.path.isfile(os.path.join(data_directory,'AA_matrices',filename))==False:
            sys.exit("".join(['AA matrix ', aa_matrix, ' does not exist. ',
                              'Check spelling and use all lower cases. ',
                       'You can also define your custom matrix as a 20-by-20 ',
                        'pandas DataFrame with AAs as row and column names.']))
        else: #load stored AA matrix
            aa_distance_matrix = pd.read_csv(os.path.join(data_directory,
                                                          'AA_matrices',
                                                          filename),
                                             index_col=0);    
    
    
    # save row and column AA orders for row to column AA substitutions
    from_AA = np.array(aa_distance_matrix.index).astype(str)
    to_AA = np.array(aa_distance_matrix.columns).astype(str)
    
    # Convert to np matrix for easy indexing
    aa_distance_matrix = aa_distance_matrix.to_numpy()
    
    #infer peptide length
    if isinstance(index_peptide, str): #One index peptide provided
        peptide_length = len(index_peptide)
    else: #List of index peptides provided
        peptide_length = len(index_peptide[0]) 
    
    # To locate mutation positions together in all mutant sequences, join them all
    if isinstance(index_peptide, str): #One index peptide provided
        index_peptide_repeated = np.array(list("".join(list(np.repeat(index_peptide,
                                            len(mutant_peptide_list))))))
    else: #List of index peptide provided
        index_peptide_repeated = np.array(list("".join(index_peptide)))
        
    mutant_peptides_joined = np.array(list("".join(mutant_peptide_list)))
    
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
    aa_distance = np.reshape(aa_distance, 
                             (len(mutant_peptide_list),peptide_length))
    
    
    return aa_distance

##########################################################################

'''
Function to generate peptide to index distance by multiplying positional weight
profile with features of a list of mutant peptides, created using their
position-dependent distances from the index peptide,
based on a given AA distance matrix.

Index peptide can be a string (common for all mutants) or a list (potentially
different for all mutants) with length equal to that of mutant list. All input 
sequences must have equal length.

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

Weight profile should be array of size peptide_length-by-1 (common for all) or 
peptide_length-by-#mutants (different weights for different mutants)

'''
def peptide2index(index_peptide, # Single index peptide or list 
                 mutant_peptide_list, #List of mutant peptide sequences
                 aa_matrix, #Named or user-defined AA matrix
                 weight_profile): #array of weights
    
    import numpy as np
    import sys
    
    
    # Create mutant features
    mutant_features = generate_mutant_features(index_peptide, 
                                               mutant_peptide_list, 
                                               aa_matrix)
    
    # Check if weight profile is an array with desired size
    peptide_length = len(mutant_peptide_list[0])    
    if ((weight_profile.dtype not in ('int','float')) or
        (np.shape(weight_profile) not in ((1,peptide_length),
                                (len(mutant_peptide_list),peptide_length)))):
        sys.exit("".join(["Weight profile must be a numerical array of shape ",
                  "(1,peptide length) or (number of mutants,peptide length)"]))
    
    # Calculate peptide-to-index distances by multiplication
    if (np.shape(weight_profile)==(1,peptide_length)):
        distance = mutant_features.dot(weight_profile.transpose())[:,0]
    else:
        distance = np.multiply(mutant_features,weight_profile).sum(axis=1)
    
    return distance
#############################################################################
'''
Function that outputs sampled positional weights after pooled inference with peptide
features and TCR activation categories
'''


def pooled_inference_weights_only(peptide_feature, # features for each peptide
                                  peptide_activity,#activation category
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np
    
      
    peptide_length = np.shape(peptide_feature)[1]
    
    with pm.Model() as peptide_regression_model:
        # positional weights, pooled over TCRs and positions
        weights = pm.Beta("weights",
                          alpha = pm.Gamma("alpha_w", mu=1.0, sigma=0.5),
                          beta = pm.Gamma("beta_w", mu=5.0, sigma=1.0),
                          shape=(1,peptide_length))        
        
        # Hyperprior parameters for intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)
        
        #Full predictor
        # baseline
        mu = pm.math.sum(peptide_feature*weights[0,:],axis=1) 
        #positional weight * peptide feature summed over position
            
        mu = mu + intercepts 
        
        sigma = pm.Exponential("sigma", lam=1)
            
        # Gaussian Regression
        peptide_activity_obs = pm.Normal("peptide_activity_obs",
                                        mu=mu,
                                        sigma=sigma,
                                        observed=peptide_activity)
        
    # Sampling with approximate posterior
    with peptide_regression_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000,
                                                             random_seed=seed))
     
 
    # Extract position-dependent weights of TCRs
    # Find start index of inferred weight parameters
    weight_index = np.argwhere(inferred_params.index=='weights[0, 0]')[0,0]
        
    # extract and reshape TCR-specific positional weights
    inferred_weights=np.reshape(inferred_params.iloc[
        weight_index:weight_index+peptide_length,
        0].to_numpy(),newshape=(1,peptide_length),order='C') #TCR-by-position
                
    return inferred_weights 
    
#######################################################################

'''
Function to train BATMAN, using TCR activation data provided in a csv file.
Refer to the example file to see format of the input TCR activation datafile.

TCR activation levels must start from 0 (weakest) and are integers 
(no missing level allowed)

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

Runs in one of 3 modes: weights_only, full, or symm, based on if only weights
are inferred (AA matrix being the one specified) or if both weights and AA matrices
are inferred (symmetric or full AA matrix)

(Optional) NetMHCPan pMHC binding info can be input as one of: SB, WB, NB
(Optional) pMHC features under columns 'pmhc_feature_*' for regression

Also takes #steps for sampling and a seed argument for reproduction 

'''

def train(filename,#Path to file or pandas df with TCR data (see example for format)
          # Default values of optional parameters
          aa_matrix='BLOSUM100',#Named or user-defined AA matrix used for regularization
          seed=100,#seed for sampling
          steps=20000):# number of steps for sampling


    import pandas as pd
    import numpy as np
    import os.path, sys
    
    '''Check input'''
    
    # If input is not a pandas dataframe
    if isinstance(filename, pd.DataFrame)==False:
        # Check if file exists and is csv
        if os.path.isfile(filename)==False:
            sys.exit("File or DataFrame does not exist. Check filename and directory.")
        
        if filename.endswith('.csv')==False:
            sys.exit("Input must be csv file or pandas DataFrame. See example input for details.")
    
        # Read file
        peptide_data = pd.read_csv(filename)
        
    else: # If input is a pandas DataFrame
        peptide_data = filename.copy() # Read DataFrame    
    
    # Check that input has correct headers for relevant columns
    if {'peptide_activity', 'index', 'peptide', 'tcr'}.issubset(
            set(peptide_data.columns))==False:
        sys.exit(''.join(["Input must have these 4 headers:", 
                 "'peptide_activity', 'index', 'peptide', 'tcr'.",
                 "Check spelling and have all headers in lower case."]))    
    
    
    '''featurize peptide data'''        
    # Assign indices to unique TCRs in data
    
    # Make list of index peptide
    index_list = peptide_data['index'].tolist()
    # Make list of mutant peptide
    peptide_list = peptide_data['peptide'].tolist()
    
    # Make peptide features
    peptide_feature = generate_mutant_features(index_list,
                                               peptide_list,
                                               aa_matrix)
    
    # Peptide activation categories
    peptide_activity = peptide_data['peptide_activity'].to_numpy()
    
    # Infer weights    
    inferred_weights = pooled_inference_weights_only(peptide_feature, # features for each peptide
                                      peptide_activity,
                                      seed=seed, steps=steps)
    
    return inferred_weights
    ##############################################################################
