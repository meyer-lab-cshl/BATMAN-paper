import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
import sys
import os

# Join path to BATMAN functions file
sys.path.append('../../../batman_functions')
from batman_functions import *


#%% '''load data'''
tcr_pmhc_data = pd.read_excel(''.join(['../../../tcr_epitope_datasets/',
          'mutational_scan_datasets/train_test_data_folds.xlsx']),index_col=0)
tcr_pmhc_data.rename(columns={'index_peptide': 'index'}, inplace=True)

# Subset selected TCRs
TCR_list_9mer = ['a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11'] 

TCR_list_10mer = ['TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va']

TCR_list_mhcii = ['TCR-F5','TCR-3598-2','MBP-TCR','B3K508']

TCR_list = TCR_list_9mer + TCR_list_10mer + TCR_list_mhcii

tcr_pmhc_data = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list)]

#%%'''Run within TCR predictions: pooled on 9mer-binding TCRs'''
# Pool within 9-mer-binding TCRs
tcr_pmhc_data_9mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_9mer)].copy()
aucs_9mer = np.zeros((len(TCR_list_9mer),5)) #Array to store AUCs: #TCRs X #folds

# Loop over 5 folds
for fold in np.arange(5):
    train_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.training_fold!=fold].copy()
    test_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.training_fold==fold].copy()
    
    # Add index peptide in training data for all fold
    train_data = pd.concat([train_data,
                            tcr_pmhc_data_9mer[pd.isna(tcr_pmhc_data_9mer.training_fold)]])

    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_mhc_effect = train(train_data,
                                                'weights_only', # no AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Get test set mutant2index distance and MHC effect as TCR activation predictor
    test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              'blosum100', 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                          test_data['BindLevel']
                                                          ].to_numpy())        
    
    # Calculate TCR-specific AUC, averaged over 3 class-pairs
    for tcr_index in np.arange(len(TCR_list_9mer)):
        tcr = TCR_list_9mer[tcr_index]
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        auc_tcr = np.zeros((1,3))
        
        # Loop over class pairs
        for i in np.arange(0,3):
            test_data_tcr = test_data[test_data.tcr==tcr].copy()
            fpr[i], tpr[i], _ = roc_curve(test_data_tcr.activation[test_data_tcr.activation !=i],
                                          -test_data_tcr.mutant_distance[test_data_tcr.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            auc_tcr[0,i] = auc(fpr[i], tpr[i])
        
        # Record average AUC
        aucs_9mer[tcr_index,fold] = np.nanmean(auc_tcr)#if 2 classes only, avoids nan


# Average over folds and save
aucs_9mer = np.mean(aucs_9mer,axis=1)    

'''###############################################################'''

#%%'''Run within TCR predictions: pooled on 10-mer-binding TCRs'''
# Pool within 10-mer-binding TCRs
tcr_pmhc_data_10mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_10mer)].copy()
aucs_10mer = np.zeros((len(TCR_list_10mer),5)) #Array to store AUCs: #TCRs X #folds

# Loop over 5 folds
for fold in np.arange(5):
    train_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.training_fold!=fold].copy()
    test_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.training_fold==fold].copy()
    
    # Add index peptide in training data for all fold
    train_data = pd.concat([train_data,
                            tcr_pmhc_data_10mer[pd.isna(tcr_pmhc_data_10mer.training_fold)]])

    # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
    inferred_weights,inferred_mhc_effect = train(train_data,
                                                'weights_only', # No AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Get mutant2index distance from test set as TCR activation predictor
    test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              'blosum100', 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                          test_data['BindLevel']
                                                          ].to_numpy()) 
    
    # Calculate TCR-specific AUC, averaged over 3 class-pairs
    for tcr_index in np.arange(len(TCR_list_10mer)):
        tcr = TCR_list_10mer[tcr_index]
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        auc_tcr = np.zeros((1,3))
        
        # Loop over class pairs
        for i in np.arange(0,3):
            test_data_tcr = test_data[test_data.tcr==tcr].copy()
            fpr[i], tpr[i], _ = roc_curve(test_data_tcr.activation[test_data_tcr.activation !=i],
                                          -test_data_tcr.mutant_distance[test_data_tcr.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            auc_tcr[0,i] = auc(fpr[i], tpr[i])
        
        # Record average AUC
        aucs_10mer[tcr_index,fold] = np.nanmean(auc_tcr) #if 2 classes only, avoids nan


# Average over folds and save
aucs_10mer = np.mean(aucs_10mer,axis=1)    

'''##############################################################'''
#%% '''Run within TCR predictions: MHCII TCRs (unpooled only)'''
# Array to store results
aucs_mhcii = np.zeros((len(TCR_list_mhcii),5)) #Array to store AUCs: #TCRs X #folds

# Calculate TCR-specific AUC, individually for each MHCII TCR
for tcr_index in np.arange(len(TCR_list_mhcii)):
    tcr = TCR_list_mhcii[tcr_index]
    tcr_pmhc_data_mhcii = tcr_pmhc_data[tcr_pmhc_data.tcr==tcr].copy()
    
    # Loop over 5 folds
    for fold in np.arange(5):
        train_data = tcr_pmhc_data_mhcii[tcr_pmhc_data_mhcii.training_fold!=fold].copy()
        test_data = tcr_pmhc_data_mhcii[tcr_pmhc_data_mhcii.training_fold==fold].copy()
        
        # Add index peptide in training data for all fold
        train_data = pd.concat([train_data,
                                tcr_pmhc_data_mhcii[pd.isna(
                                    tcr_pmhc_data_mhcii.training_fold)]])
    
        # Run BATMAN to infer AA matrix and TCR-specific weights and MHC effect
        inferred_weights,inferred_mhc_effect = train(train_data,
                                                    'weights_only', # No AA matrix
                                                    'blosum100', #AA matrix prior
                                                    consider_mhc_binding=True,
                                                    steps = 80000,
                                                    seed = 111)
        
        # Get mutant2index distance from test set as TCR activation predictor
        test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                                  test_data['peptide'].tolist(),
                                  'blosum100', 
                                  inferred_weights.loc[test_data['tcr'],:].to_numpy()
                                  ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                              test_data['BindLevel']
                                                              ].to_numpy())  
    
        # Calculate pairwise AUC for mutant2index distances
        fpr = dict()
        tpr = dict()
        auc_tcr = np.zeros((1,3))
        
        # Loop over class pairs
        for i in np.arange(0,3):
            fpr[i], tpr[i], _ = roc_curve(test_data.activation[test_data.activation !=i],
                                          -test_data.mutant_distance[test_data.activation !=i],
                                          pos_label=max(set(np.arange(0,3)).difference({i})))
            auc_tcr[0,i] = auc(fpr[i], tpr[i])
        
        # Record average AUC
        aucs_mhcii[tcr_index,fold] = np.nanmean(auc_tcr) #if 2 classes only, avoids nan

# Average over folds and save
aucs_mhcii = np.mean(aucs_mhcii,axis=1)    

# All AUCs: Within-TCR
aucs_all_within=pd.DataFrame(np.concatenate((aucs_9mer,aucs_10mer,aucs_mhcii),axis=0),
                      index=TCR_list,columns=['AUC_within_TCR'])


'''##############################################################'''
#%%'''Run LOO-TCR predictions: 9mer-binding TCRs'''
# Subset to 9-mer-binding TCRs
tcr_pmhc_data_9mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_9mer)].copy()
# Array to store results
aucs_9mer = np.zeros((len(TCR_list_9mer),1)) #Array to store AUCs: #TCRs

# Calculate TCR-specific AUC, individually for each MHCII TCR
for tcr_index in np.arange(len(TCR_list_9mer)):
    tcr = TCR_list_9mer[tcr_index]    
    
    # Loo TCR data
    train_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.tcr!=tcr].copy()
    test_data = tcr_pmhc_data_9mer[tcr_pmhc_data_9mer.tcr==tcr].copy()
    
    # since it is pan-TCR, discard all TCR names
    train_data.tcr = 'pan-TCR'
    test_data.tcr = 'pan-TCR'
    
    # Run BATMAN to infer pan-TCR weights and MHC effect
    inferred_weights,inferred_mhc_effect = train(train_data,
                                                'weights_only', # No AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Get mutant2index distance from test set as TCR activation predictor
    test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              'blosum100', 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                          test_data['BindLevel']
                                                          ].to_numpy())  

    # Calculate pairwise AUC for mutant2index distances
    fpr = dict()
    tpr = dict()
    auc_tcr = np.zeros((1,3))
    
    # Loop over class pairs
    for i in np.arange(0,3):
        fpr[i], tpr[i], _ = roc_curve(test_data.activation[test_data.activation !=i],
                                      -test_data.mutant_distance[test_data.activation !=i],
                                      pos_label=max(set(np.arange(0,3)).difference({i})))
        auc_tcr[0,i] = auc(fpr[i], tpr[i])
    
    # Record average AUC
    aucs_9mer[tcr_index,0] = np.nanmean(auc_tcr) #if 2 classes only, avoids nan


'''##############################################################'''
#%%'''Run LOO-TCR predictions: 10mer-binding TCRs'''
# Subset to 9-mer-binding TCRs
tcr_pmhc_data_10mer = tcr_pmhc_data[tcr_pmhc_data.tcr.isin(TCR_list_10mer)].copy()
# Array to store results
aucs_10mer = np.zeros((len(TCR_list_10mer),1)) #Array to store AUCs: #TCRs

# Calculate TCR-specific AUC, individually for each MHCII TCR
for tcr_index in np.arange(len(TCR_list_10mer)):
    tcr = TCR_list_10mer[tcr_index]    
    
    # Loo TCR data
    train_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.tcr!=tcr].copy()
    test_data = tcr_pmhc_data_10mer[tcr_pmhc_data_10mer.tcr==tcr].copy()
    
    # since it is pan-TCR, discard all TCR names
    train_data.tcr = 'pan-TCR'
    test_data.tcr = 'pan-TCR'
    
    # Run BATMAN to infer pan-TCR weights and MHC effect
    inferred_weights,inferred_mhc_effect = train(train_data,
                                                'weights_only', # No AA matrix
                                                'blosum100', #AA matrix prior
                                                consider_mhc_binding=True,
                                                steps = 80000,
                                                seed = 111)
    
    # Get mutant2index distance from test set as TCR activation predictor
    test_data['mutant_distance'] = peptide2index(test_data['index'].tolist(), 
                              test_data['peptide'].tolist(),
                              'blosum100', 
                              inferred_weights.loc[test_data['tcr'],:].to_numpy()
                              ) + np.diag(inferred_mhc_effect.loc[test_data['tcr'],
                                                          test_data['BindLevel']
                                                          ].to_numpy()) 

    # Calculate pairwise AUC for mutant2index distances
    fpr = dict()
    tpr = dict()
    auc_tcr = np.zeros((1,3))
    
    # Loop over class pairs
    for i in np.arange(0,3):
        fpr[i], tpr[i], _ = roc_curve(test_data.activation[test_data.activation !=i],
                                      -test_data.mutant_distance[test_data.activation !=i],
                                      pos_label=max(set(np.arange(0,3)).difference({i})))
        auc_tcr[0,i] = auc(fpr[i], tpr[i])
    
    # Record average AUC
    aucs_10mer[tcr_index,0] = np.nanmean(auc_tcr) #if 2 classes only, avoids nan


# All AUCs: LOO-TCR
aucs_all_loo=pd.DataFrame(np.concatenate((aucs_9mer,aucs_10mer),axis=0),
                      index=TCR_list_9mer+TCR_list_10mer,columns=['AUC_LOO_TCR'])


#%%# All AUCs
aucs_all = pd.concat([aucs_all_within,aucs_all_loo],axis=1)
aucs_all.index.name = 'TCR'
# Save data
aucs_all.to_csv('BATMAN_BLOSUM100_AUCs.csv')








