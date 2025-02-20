import numpy as np
import pandas as pd
from PIL import Image
import itertools


# Load image
im = Image.open('fig5_ddg_heatmap.jpeg','r') 
pix = list(im.getdata())
width, height = im.size
pixel_values = np.array(pix).reshape((width, height, 3))

# Set up a grid to extract RGB values for each box in heatmap
xmin = 94
xmax = 735
ymin = 566
ymax = 832
xpts = 20
ypts = 7

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
heatmap_rgb = np.zeros((140,3))

for pixel in np.arange(140):
    heatmap_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Reshape to heatmap size y channel
heatmap_rgb = np.reshape(heatmap_rgb, (20,7,3))


#Colorbar
# sample 100 equally spaced collinear points across the colorbar

xmin = 825
xmax = 826
ymin = 859
ymax = 548.6
xpts= 1
ypts = 100

grid = list(itertools.product(np.linspace(xmin,xmax,xpts),
                              np.linspace(ymin,ymax,ypts)))

# Extract RGB values
cb_rgb = np.zeros((ypts,3))

for pixel in np.arange(ypts):
    cb_rgb[pixel,:] = im.getpixel(grid[pixel])[0:3]

# Interpolate DDG scale
rgb2DDG = np.linspace(4.5,-0.2,ypts)

#Assign CD69 values based on RGB
DDG_values = np.zeros((7,20))

for position in np.arange(7):
    for aa in np.arange(20):        
        # RGB value in heatmap
        rgb = heatmap_rgb[aa,position,:] #Read recorded rgb value
        # Find closest RGB combination in the colorbar
        cb_loc=np.argmin(np.sum(abs(cb_rgb-np.tile(rgb,(len(cb_rgb),1))),axis=1))
        DDG_values[position,aa] = rgb2DDG[cb_loc] #assign CD69

#Record data
tcr_data = pd.DataFrame({
'tcr':['TIL1383I '],
'index_peptide':["YMDGTMSQV"],
})

index_peptide = "YMDGTMSQV"
aa_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y']
mutant_DDG = pd.DataFrame(columns=['peptide','Kd'])
mutation_positions = [0,2,3,4,5,6,7]
for position in np.arange(7):
    for aa in np.arange(len(aa_order)):
        peptide = list(index_peptide)
        peptide[mutation_positions[position]] = aa_order[aa]
        peptide = ''.join(peptide)
        peptide_activity = DDG_values[position,aa]
        data_add = pd.DataFrame({'peptide':[peptide],
                                 'Kd':[peptide_activity]})
        mutant_DDG = pd.concat([mutant_DDG,data_add])

# Remove multiple instances of the index
mutant_DDG.Kd[mutant_DDG.peptide==index_peptide] = 0 #Since DDG is wrt index
mutant_DDG = mutant_DDG.drop_duplicates().reset_index(drop=True)

# Change DDG data to Kd by using Kd = exp(DG/RT), with DG = DG_index + DDG
# DG of index peptide is -6.25 kcal/mol from Fig. 6 table
# T is assumed to be 298K, so RT=0.592 (this also matches Fig. 6 DG and Kd conversions)
mutant_DDG['Kd'] = np.exp((mutant_DDG['Kd'] - 6.25)/0.592)


# Merge all data columns
tcr_data = tcr_data.loc[tcr_data.index.repeat(len(mutant_DDG))
                        ].reset_index(drop=True)
tcr_data = pd.concat([tcr_data,mutant_DDG],axis=1)

# Export data
tcr_data.to_excel('Kd_TIL1383I.xlsx', index=False)
        
        





















