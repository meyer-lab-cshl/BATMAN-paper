import numpy as np
import pandas as pd
from PIL import Image
import itertools
from matplotlib import pyplot as plt


# Load image
im = Image.open('heatmaps_1_hamming.jpg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

'''APN-responsive TCR'''
# Set up a grid to extract RGB values for each box in heatmap
xmin = 107
xmax = 450
ymin = 60
ymax = 557
xpts = 10
ypts = 20

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
APN_rgb = np.zeros((200,3))

for pixel in np.arange(200):
    APN_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
APN_rgb = np.reshape(APN_rgb, (10,20,3))

'''EWW-responsive TCR'''
# Set up a grid to extract RGB values for each box in heatmap
xmin = 640
xmax = 983
ymin = 60
ymax = 557
xpts = 10
ypts = 20

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
EWW_rgb = np.zeros((200,3))

for pixel in np.arange(200):
    EWW_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
EWW_rgb = np.reshape(EWW_rgb, (10,20,3))

'''Colorbar'''
# sample 100 equally spaced collinear points between upper (4) & lower (-6) ticks

xmin = 1050
xmax = 1051
ymin = 436.5
ymax = 565.5
xpts= 1
ypts = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
cb_rgb = np.zeros((ypts,3))

for pixel in np.arange(ypts):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Interpolate logFC scale
rgb2log_FC = np.linspace(4,-6,ypts)

'''Assign logFC values based on RGB'''
APN_logFC = np.zeros((10,20))
EWW_logFC = np.zeros((10,20))

for position in np.arange(10):
    for aa in np.arange(20):
        
        # APN-TCR
        rgb = APN_rgb[position,aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        APN_logFC[position,aa] = rgb2log_FC[cb_loc] #assign logFC
        
        # EWW-TCR
        rgb = EWW_rgb[position,aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        EWW_logFC[position,aa] = rgb2log_FC[cb_loc] #assign logFC

'''Normalize to [0,1]'''

#First set negative logFC to zero
APN_logFC[APN_logFC<=0]=0
EWW_logFC[EWW_logFC<=0]=0


APN_logFC = (APN_logFC - np.min(APN_logFC))/(np.max(APN_logFC) - np.min(APN_logFC))        
EWW_logFC = (EWW_logFC - np.min(EWW_logFC))/(np.max(EWW_logFC) - np.min(EWW_logFC))        

'''Record data'''

# APN-TCR
tcr_data_APN = pd.DataFrame({
'tcr_name':['APN-TCR'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["DNA barcode enrichment"],
'tcr_source_organism':["human"],
'index_peptide':["APNCYGNIPL"],
'index_peptide_activity':[APN_logFC[0,0]],
'mhc':["HLA-B*07:02"],
'pmid	':["30451992"],
'peptide_type':["cancer antigen"]
})

index_peptide = "APNCYGNIPL"
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_logFC = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(10):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = APN_logFC[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_logFC = pd.concat([mutant_logFC,data_add])

mutant_logFC = mutant_logFC.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data_APN = tcr_data_APN.loc[tcr_data_APN.index.repeat(len(mutant_logFC))
                        ].reset_index(drop=True)
tcr_data_APN = pd.concat([tcr_data_APN,mutant_logFC],axis=1)

# EWW-TCR
tcr_data_EWW = pd.DataFrame({
'tcr_name':['EWW-TCR'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["DNA barcode enrichment"],
'tcr_source_organism':["human"],
'index_peptide':["EWWRSGGFSF"],
'index_peptide_activity':[EWW_logFC[0,3]],
'mhc':["HLA-A*24:02"],
'pmid	':["30451992"],
'peptide_type':["cancer antigen"]
})

index_peptide = "EWWRSGGFSF"
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_logFC = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(10):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = EWW_logFC[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_logFC = pd.concat([mutant_logFC,data_add])

mutant_logFC = mutant_logFC.drop_duplicates().reset_index(drop=True)


# Merge all data columns
tcr_data_EWW = tcr_data_EWW.loc[tcr_data_EWW.index.repeat(len(mutant_logFC))
                        ].reset_index(drop=True)
tcr_data_EWW = pd.concat([tcr_data_EWW,mutant_logFC],axis=1)

tcr_data = pd.concat([tcr_data_APN,tcr_data_EWW])


# Export data
#tcr_data.to_excel('epitope_data_fingerprint_TCRs.xlsx', index=False)
        
        





















