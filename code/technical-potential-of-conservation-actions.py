# !/bin/env python3
#
# technical-potential-of-conservation-actions.py 
#
# This script scopes the utility of "technical abatement/mitigation potential" analysis to quantify the potential 
# of different conservation action to meet conservation and sustainability targets.
# 
# By Carla Arcihbald (c.archibaldc@deakin.edu.au)
# Created: 2022-02-02
# Last modified: 2022-02-02
#
#
################################################

################ Set up ########################
# Load python modules
import numpy as np
import pandas as pd
import rasterio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Input files
# Cell_df administrative zones
cell_zones_df_file = "N://Planet-A//Data-Master//LUTO_2.0_input_data//Input_data//2D_Spatial_Snapshot//cell_zones_df.h5"
cell_zones_df = pd.read_hdf(cell_zones_df_file)

# NLUM
infile = "N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif"

# Annual Native title determinations database, this is used allong with Indigenou Tenure to identify Indigenous lands
#in_cell_df_path_3 = "\\\school-les-m.shares.deakin.edu.au\\school-les-m\\Planet-A\\Data-Master\\Indigenous_lands\\cell_indigenous_df.pkl"

# Set graph colour pallette
cmap = ListedColormap(['#E0F5EB', '#34a169','white'])

################ Create some helper data and functions ########################
 
# View columns and unique values
for col in cell_zones_df.columns:
    print(col)

dfUnique = pd.DataFrame(cell_zones_df.PRIMARY_V7.unique())  
pd.set_option('display.max_rows', dfUnique.shape[0]+1)
dfUnique

pd.set_option('display.max_rows', 200)
cell_df.groupby(['NNTT_EXCLUSIVE_TITLE_2020','SPREAD_DESC'], observed=True)['SPREAD_DESC'].count().head(100)

# Open NLUM_ID as mask raster and get metadata
with rasterio.open(infile) as rst:
    
    # Load a 2D masked array with nodata masked out
    NLUM_ID_raster = rst.read(1, masked=True) 
    
    # Create a 0/1 binary mask
    NLUM_mask = NLUM_ID_raster.mask == False
    
    # Get metadata and update parameters
    NLUM_transform = rst.transform
    NLUM_crs = rst.crs
    meta = rst.meta.copy()
    meta.update(compress='lzw', driver='GTiff') # , dtype='int32', nodata='0')
    [meta.pop(key) for key in ['dtype', 'nodata']] # Need to add dtype and nodata manually when exporting GeoTiffs
    
    # Set up some data structures to enable conversion on 1D arrays to 2D
    array_2D = np.zeros(NLUM_ID_raster.shape)
    xy = np.nonzero(NLUM_mask)
    
# Convert 1D column to 2D spatial array
def conv_1D_to_2D(in_1D_array):
    array_2D[xy] = np.array(in_1D_array)
    array_2D[array_2D == 0] = np.nan
    return array_2D

# Convert 1D column to 2D spatial array and plot map
def map_in_2D(col, data): # data = 'continuous' or 'categorical'

    a2D = conv_1D_to_2D(col)
    if data == 'categorical':
        n = col.nunique()
        cmap = matplotlib.colors.ListedColormap(np.random.rand(n,3))
        plt.imshow(a2D, cmap=cmap, resample=False)
    elif data == 'continuous':
        plt.imshow(a2D, cmap='pink', resample=False)
    plt.show()

# Convert object columns to categories and downcast int64 columns to save memory and space
def downcast(dframe):
    obj_cols = dframe.select_dtypes(include = ["object"]).columns
    dframe[obj_cols] = dframe[obj_cols].astype('category')
    int64_cols = dframe.select_dtypes(include = ["int64"]).columns
    dframe[int64_cols] = dframe[int64_cols].apply(pd.to_numeric, downcast = 'integer')
    
# Process land use sieve mask                   
def create_lu_sieve_mask(name,index,dframe):  
#def create_lu_sieve_mask(name,index,outfile,dframe):  
                    
    # Add a new field to cell_df and calculate values for selected rows using pandas
    dframe[name] = 0
    dframe[name] = dframe[name].astype(np.uint8)
    dframe.loc[index,name] = 1
    
    # Open NLUM mask raster and get metadata
    with rasterio.open(infile) as src:
        
        # Read geotiff to numpy array, loads a 2D masked array with nodata masked out
        NLUM_mask = src.read(1) 
        
        # Update the metadata
        meta = src.meta.copy()
        meta.update(compress = 'lzw', driver = 'GTiff', dtype = 'uint8', nodata = 255)
        
        # Convert 1D column created above from the cell_df to a 2D array
        np2dArray = np.zeros(NLUM_mask.shape) + 255
        xy = np.nonzero(NLUM_mask)
        np2dArray[xy] = dframe[name]
 
        # Plot GeoTiff
        plt.imshow(np2dArray, cmap=cmap, vmin=0, vmax=2)
        
        # Save the output file as GeoTiff
        # r'N:\Planet-A\LUF-Modelling\Land-Use-Sieve\Data\Masks\SOL_Conservation_Government-protected-area_Potential-area-mask_1km.tif'
        #with rasterio.open(outfile, "w+", **meta) as dst:
        #    dst.write(np2dArray.astype(np.uint8), indexes = 1)
                      
################ Determine spatial footprint of consevration action ########################
#
# Heuristic rues are developed for each conservation action based upon land tenure, land cover and land use. 
# 
# --- [PR-GPA] Government protected areas
#              Heuristic: Government protected areas can occur on all land exept freehold land, and must 
PrGpaIndex = cell_zones_df.query("(TENURE_DESC in ['Public land - other','Nature conservation areas']) and \
                                  (PRIMARY_V7 in ['1 Conservation and natural environments'])").index

reservePublicMask = create_lu_sieve_mask('PR-GPA_MASK', PrGpaIndex, cell_zones_df)

# --- [PR-PPA] Private land conservation (Ecosystem regeneration and protection on private land)
PrPpaIndex  = cell_zones_df.query("(TENURE_DESC in ['Private land']) and \
                                   (VEG_COND_DESC == ['No modification of native vegetation'])").index

reservePrivateMask = create_lu_sieve_mask('PR-PPA_MASK', PrPpaIndex, cell_zones_df)

# --- [PR-IPA] Indigenous protected areas
PrIpaIndex  = cell_zones_df.query("(TENURE_DESC in ['Aboriginal land - traditional indigenous uses','Aboriginal land - other non-agricultural']) and \
                                   (VEG_COND_DESC == ['No modification of native vegetation'])").index

reserveIndigenousMask = create_lu_sieve_mask('PR-IPA_MASK', PrIpaIndex, cell_zones_df)

# --- [PR-AHL] Avoided habitat loss
PrAhlIndex = cell_zones_df.query("(TENURE_DESC in ['Private land']) and \
                                  (PRIMARY_V7 in ['2 Production from relatively natural environments'])").index

reservePrivateMask = create_lu_sieve_mask('PR-PPA_MASK', PrAhlIndex, cell_zones_df)

# --- [PR-NRV] Native revegetation
PrNrvIndex  = cell_zones_df.query("(TENURE_DESC in ['Private land']) and \
                                   (VEG_COND_DESC == ['Modification of native vegetation'])").index

natiRevegMask = create_lu_sieve_mask('PR-NRV_MASK', PrNrvIndex, cell_zones_df)

# --- [PR-RIP] Riparian restoration
#PrRipIndex = cell_zones_df.query("(RIPARIAN_WATER_BODY >= 100) or \
#                                  (RIPARIAN_RIVERS_STRA >= 1) and \
#                                  (STUDY_AREA_MASK == 1)").index

#ripManMask = create_lu_sieve_mask('PR-RIP_MASK', PrRipIndex, cell_df)

################ Determine technical potential of consevration action to contribute towards conservation and sustainability ########################


# End :) 