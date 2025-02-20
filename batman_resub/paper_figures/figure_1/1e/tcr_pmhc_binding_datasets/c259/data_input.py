import numpy as np
import pandas as pd
from PIL import Image
import itertools
from matplotlib import pyplot as plt


# Load image
im = Image.open('F6.large.jpg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

# Set up a grid to extract RGB values for each box in heatmap
xmin = 448
xmax = 553
ymin = 98
ymax = 385
xpts = 9
ypts = 19

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
TCR_rgb = np.zeros((xpts*ypts,3))

for pixel in np.arange(xpts*ypts):
    TCR_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
TCR_rgb = np.reshape(TCR_rgb, (xpts,ypts,3))


'''Colorbar'''
# sample 100 equally spaced collinear points between upper & lower ticks

xmin = 870
xmax = 870
ymin = 392
ymax = 94
xpts= 1
ypts = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
cb_rgb = np.zeros((ypts,3))

for pixel in np.arange(ypts):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Interpolate colorbar scale
rgb2log_FC = np.linspace(0,2*1E6,ypts)

'''Assign killing values based on RGB'''
TCR_killing = np.zeros((9,19))

for position in np.arange(9):
    for aa in np.arange(19):
        
        # APN-TCR
        rgb = TCR_rgb[position,aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        TCR_killing[position,aa] = rgb2log_FC[cb_loc] #assign logFC
        
'''Record data'''

# APN-TCR
tcr_data = pd.DataFrame({
'tcr_name':['c259'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["CAVRPLYGGSYIPTF"],
'cdr3b':["CASSYVGNTGELFF"],
'trav':["21*01"],
'traj':["6*01"],
'trbv':["6-5*01"],
'trbd':["NA"],
'trbj':["2-2*01"],
'assay':["killing"],
'tcr_source_organism':["human"],
'index_peptide':["SLLMWITQV"],
'mhc':["HLA-A*02:01"],
'pmid	':["TBA"],
'peptide_type':["cancer antigen"]
})

index_peptide = "SLLMWITQV"
aa_order = ['A','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_killing = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(9):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = TCR_killing[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_killing = pd.concat([mutant_killing,data_add])

# the index peptide has slightly different values due to RGB conversion, average them
mutant_killing.peptide_activity[
    mutant_killing.peptide==index_peptide
    ]=np.mean(mutant_killing.peptide_activity[mutant_killing.peptide==index_peptide])
mutant_killing = mutant_killing.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data= tcr_data.loc[tcr_data.index.repeat(len(mutant_killing))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_killing],axis=1)


# Export data
tcr_data.to_excel('epitope_data_c259.xlsx', index=False)
        
        





















