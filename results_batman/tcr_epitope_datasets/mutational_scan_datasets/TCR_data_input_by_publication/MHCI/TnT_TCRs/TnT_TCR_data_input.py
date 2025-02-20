import numpy as np
import pandas as pd
from PIL import Image
import itertools
from matplotlib import pyplot as plt


# Load image
im = Image.open('heatmaps_lrg.jpg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

'''A3-05 TCR'''
# Set up a grid to extract RGB values for each box in heatmap
xmin = 60
xmax = 458
ymin = 200
ymax = 829
xpts = 9
ypts = 20

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
a305_rgb = np.zeros((180,3))

for pixel in np.arange(180):
    a305_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
a305_rgb = np.reshape(a305_rgb, (9,20,3))


#Colorbar

# The length of the colorbar is measured to be (ymax-ymin)=843-186=657 px
# of which length till 100 mark is (ymax-y_100) = 843-303 =540 px, 
# so total colorbar is (100/540)*657=121.667 in %CD69 unit
# sample 100 equally spaced collinear points across the colorbar

xmin = 509
xmax = 510
ymin = 843
ymax = 186
xpts= 1
ypts = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
cb_rgb = np.zeros((ypts,3))

for pixel in np.arange(ypts):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Interpolate CD69 scale
rgb2CD69 = np.linspace(0,121.667,ypts)

#Assign CD69 values based on RGB
CD69_values = np.zeros((9,20))

for position in np.arange(9):
    for aa in np.arange(20):        
        # RGB value in heatmap
        rgb = a305_rgb[position,aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        CD69_values[position,aa] = rgb2CD69[cb_loc] #assign CD69
        
#Record data
tcr_data_a305 = pd.DataFrame({
'tcr_name':['A3-05'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["CD69-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["EVDPIGHLY"],
'mhc':["HLA-A*01:01"],
'pmid':["36174557"],
'peptide_type':["cancer antigen"]
})

index_peptide = "EVDPIGHLY"
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_CD69 = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(9):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = CD69_values[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_CD69 = pd.concat([mutant_CD69,data_add])

# Remove multiple instances of the index
mutant_CD69.peptide_activity[mutant_CD69.peptide==index_peptide] = 100
mutant_CD69 = mutant_CD69.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data_a305 = tcr_data_a305.loc[tcr_data_a305.index.repeat(len(mutant_CD69))
                        ].reset_index(drop=True)
tcr_data_a305 = pd.concat([tcr_data_a305,mutant_CD69],axis=1)


'''A3-10 TCR'''
# Set up a grid to extract RGB values for each box in heatmap
xmin = 739
xmax = 1137
ymin = 200
ymax = 829
xpts = 9
ypts = 20

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
a310_rgb = np.zeros((180,3))

for pixel in np.arange(180):
    a310_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
a310_rgb = np.reshape(a310_rgb, (9,20,3))


#Colorbar

# The length of the colorbar is measured to be (ymax-ymin)=843-186=657 px
# of which length till 100 mark is (ymax-y_100) = 843-245.5 =597.5 px, 
# so total colorbar is (100/597.5)*657=109.958 in %CD69 unit
# sample 100 equally spaced collinear points across the colorbar

xmin = 509
xmax = 510
ymin = 843
ymax = 186
xpts= 1
ypts = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
cb_rgb = np.zeros((ypts,3))

for pixel in np.arange(ypts):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Interpolate CD69 scale
rgb2CD69 = np.linspace(0,109.958,ypts)

#Assign CD69 values based on RGB
CD69_values = np.zeros((9,20))

for position in np.arange(9):
    for aa in np.arange(20):        
        # RGB value in heatmap
        rgb = a310_rgb[position,aa,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        CD69_values[position,aa] = rgb2CD69[cb_loc] #assign CD69
        
#Record data
tcr_data_a310 = pd.DataFrame({
'tcr_name':['A3-10'],
'va':["NA"],
'vb':["NA"],
'cdr3a':["NA"],
'cdr3b':["NA"],
'trav':["NA"],
'traj':["NA"],
'trbv':["NA"],
'trbd':["NA"],
'trbj':["NA"],
'assay':["CD69-GFP"],
'tcr_source_organism':["human"],
'index_peptide':["EVDPIGHLY"],
'mhc':["HLA-A*01:01"],
'pmid':["36174557"],
'peptide_type':["cancer antigen"]
})

index_peptide = "EVDPIGHLY"
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_CD69 = pd.DataFrame(columns=['peptide','peptide_activity'])
for position in np.arange(9):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[position] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = CD69_values[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'peptide_activity':[peptide_activity]})
        mutant_CD69 = pd.concat([mutant_CD69,data_add])

# Remove multiple instances of the index
mutant_CD69.peptide_activity[mutant_CD69.peptide==index_peptide] = 100
mutant_CD69 = mutant_CD69.drop_duplicates().reset_index(drop=True)

# Merge all data columns
tcr_data_a310 = tcr_data_a310.loc[tcr_data_a310.index.repeat(len(mutant_CD69))
                        ].reset_index(drop=True)
tcr_data_a310 = pd.concat([tcr_data_a310,mutant_CD69],axis=1)

tcr_data = pd.concat([tcr_data_a305,tcr_data_a310]).reset_index(drop=True)

# Export data
tcr_data.to_excel('epitope_data_fingerprint_TCRs.xlsx', index=False)
        
        





















