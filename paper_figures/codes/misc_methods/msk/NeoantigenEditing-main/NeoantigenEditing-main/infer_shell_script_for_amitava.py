#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 21:32:05 2023

@author: zacharysethna
"""

import os
import numpy as np
import time
import sys
sys.path.insert(0, 'C:/Users/amita/Downloads/epitope_analysis/antigen-availability-master/antigen-availability-master/TCR-crossreactivity/final_codes/local_outputs/NeoantigenEditing-main/NeoantigenEditing-main')
from InferEpitopeDist import TCRpMHCAvidity, InferEpitopeDist

#%
import json
def save_inferred_model(av_m, m, outfile_name):
    outmodel_dict = {}
    for prop, val in m.items():
        if type(val) is float:
            outmodel_dict[prop] = val
        elif prop == 'euclid_coords':
            outmodel_dict[prop] = {aa: list(m['euclid_coords'][i]) for i, aa in enumerate(av_m.amino_acids)}
        elif prop == 'M_ab':
            outmodel_dict[prop] = {aaA + '->' + aaB: m['M_ab'][i, j] for i, aaA in enumerate(av_m.amino_acids) for j, aaB in enumerate(av_m.amino_acids)}
        elif type(val) is np.ndarray:
            outmodel_dict[prop] = list(val)
                    
#    outmodel_dict = {#'blosum62_reg': m['blosum62_reg'],
#                     'log_units': m['log_units'],
#                     'd_i': list(m['d_i']),
#                     #'euclid_coords': {aa: list(m['euclid_coords'][i]) for i, aa in enumerate(av_m.amino_acids)},
#                     'M_ab': {aaA + '->' + aaB: m['M_ab'][i, j] for i, aaA in enumerate(av_m.amino_acids) for j, aaB in enumerate(av_m.amino_acids)}}
    with open(outfile_name, 'w') as outfile:
        json.dump(outmodel_dict, outfile)        
data_folder = 'C:/Users/amita/Downloads/epitope_analysis/antigen-availability-master/antigen-availability-master/TCR-crossreactivity/final_codes/local_outputs/NeoantigenEditing-main/NeoantigenEditing-main'

cmv_peptide = 'NLVPMVATV'
tcr_pMHC_pairs = {'NLV_tcr1': cmv_peptide,
                  'NLV_tcr2': cmv_peptide,
                  'NLV_tcr3': cmv_peptide
                  }

start_time = time.time()
K_a_reg = 0.001
log_Ka_range = [-4, 4]
all_tcr_info = {}

#no wt datafiles loaded
for tcr in tcr_pMHC_pairs.keys():
    all_tcr_info[tcr] = TCRpMHCAvidity(mut_pep_avidity_infile = os.path.join(data_folder, tcr + '_reactivity_curves.tsv'), K_a_reg = K_a_reg, log_Ka_range= log_Ka_range)

print('Loaded all tcrs in %.2f seconds'%(time.time() - start_time))
infer_ep_dist = InferEpitopeDist()
#%
nlv_tcrs = {tcr_name: c_tcr for tcr_name, c_tcr in all_tcr_info.items() if tcr_name.startswith('NLV')}

#%%
infer_ep_dist.n_dim = 2
infer_ep_dist.reactivity_cutoff = 0.1
infer_ep_dist.blosum62_reg = 0

print('Starting NLV model')
nlv_model_w_reg = infer_ep_dist.infer_euclid_coord_model_from_tcrs(list(nlv_tcrs.values()), seed = 17, use_all_pep_combos = True)
#%%
outdir = 'C:/Users/amita/Downloads/epitope_analysis/antigen-availability-master/antigen-availability-master/TCR-crossreactivity/final_codes/local_outputs/NeoantigenEditing-main/NeoantigenEditing-main/outdir/'
save_inferred_model(infer_ep_dist, nlv_model_w_reg, os.path.join(outdir, 'nlv_model_w_reg_all_combos_model_bl62reg_%.0e_rt_%.0e.json'%(infer_ep_dist.blosum62_reg,infer_ep_dist.reactivity_cutoff)))

