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
import numpy as npimport 
import pandas as pd 


# Input files
# Cell_df administrative zones
cell_zones_df_file = "N://Planet-A//Data-Master//LUTO_2.0_input_data//Input_data//2D_Spatial_Snapshot//cell_zones_df.h5"
cell_zones_df = pd.read_hdf(cell_zones_df_file)[['CELL_ID','CELL_HA','X','Y','NLUM_ID','STE_NAME11','LGA_CODE','LGA_NAME10','SA2_ID','TENURE','TENURE_DESC','PRIMARY_V7','SECONDARY_V7','TERTIARY_V7','NRM_CODE','NRM_NAME','IBRA_ID','IBRA_SUB_CODE_7','IBRA_SUB_NAME_7','IBRA_REG_CODE_7','IBRA_REG_NAME_7','NVIS_EXTANT_MVG_ID','NVIS_EXTANT_MVG_NAME','NVIS_EXTANT_MVS_ID','NVIS_EXTANT_MVS_NAME','NVIS_PRE-EURO_MVG_ID','NVIS_PRE-EURO_MVG_NAME','NVIS_PRE-EURO_MVS_ID','NVIS_PRE-EURO_MVS_NAME',]]
cell_zones_df_rn = cell_zones_df.rename(columns={"STE_NAME11": "STATE_NAME", "LGA_NAME10": "LGA_NAME","SA2_ID":"SA2_CODE","TENURE":"TENURE_CODE","TENURE_DESC":"TENURE_NAME","IBRA_SUB_CODE_7":"IBRA_SUB_CODE","IBRA_SUB_NAME_7":"IBRA_SUB_NAME","IBRA_REG_CODE_7":"IBRA_REG_CODE","IBRA_REG_NAME_7":"IBRA_REG_NAME"})

cell_biophysical_df_file = "N://Planet-A//Data-Master//LUTO_2.0_input_data//Input_data//2D_Spatial_Snapshot//cell_biophysical_df.h5"
cell_biophysical_df = pd.read_hdf(cell_biophysical_df_file)[['CELL_ID',"WATER_YIELD_DR_ML_HA","WATER_YIELD_SR_ML_HA","REMNANT_VEG_T_CO2_HA","EP_BLOCK_TREES_AVG_T_CO2_HA_YR","EP_BLOCK_DEBRIS_AVG_T_CO2_HA_YR","EP_BLOCK_SOIL_AVG_T_CO2_HA_YR","EP_RIP_TREES_AVG_T_CO2_HA_YR","EP_RIP_DEBRIS_AVG_T_CO2_HA_YR","EP_RIP_SOIL_AVG_T_CO2_HA_YR","EP_BELT_TREES_AVG_T_CO2_HA_YR","EP_BELT_DEBRIS_AVG_T_CO2_HA_YR","EP_BELT_SOIL_AVG_T_CO2_HA_YR","CP_BLOCK_TREES_AVG_T_CO2_HA_YR","CP_BLOCK_DEBRIS_AVG_T_CO2_HA_YR","CP_BLOCK_SOIL_AVG_T_CO2_HA_YR","CP_BELT_TREES_AVG_T_CO2_HA_YR","CP_BELT_DEBRIS_AVG_T_CO2_HA_YR","CP_BELT_SOIL_AVG_T_CO2_HA_YR","HIR_BLOCK_TREES_AVG_T_CO2_HA_YR","HIR_BLOCK_DEBRIS_AVG_T_CO2_HA_YR","HIR_BLOCK_SOIL_AVG_T_CO2_HA_YR","HIR_RIP_TREES_AVG_T_CO2_HA_YR","HIR_RIP_DEBRIS_AVG_T_CO2_HA_YR","HIR_RIP_SOIL_AVG_T_CO2_HA_YR","RIP_LENGTH_M_CELL","SOC_T_HA_TOP_30CM","NATURAL_AREA_CONNECTIVITY","VAST_CLASS","VAST_LANDSCAPE","FD_RISK_PERC_5TH","FD_RISK_MEDIAN","FD_RISK_PERC_95TH","INVEST_RKLS","INVEST_SDR","C_FACTOR_VEG","P_FACTOR_AG"]]

cell_other_df_file = "N://Planet-A//Data-Master//LUTO_2.0_input_data//Input_data//2D_Spatial_Snapshot//cell_lu_sieve_df.pkl"
cell_other_df = pd.read_pickle(cell_other_df_file)

cell_ag_df_file = "N://Planet-A//Data-Master//Profit_map//cell_ag_data.h5"
cell_ag_df = pd.read_hdf(cell_ag_df_file)

# Join cell_df files imported above
cell_df_temp1 = cell_zones_df_rn.join(cell_biophysical_df, on='CELL_ID', how='left', lsuffix='_zones', rsuffix='_biophysical', sort=False)
cell_df_rn1 = cell_df_temp1.rename(columns={"CELL_ID_zones": "CELL_ID"})
cell_df_temp2 = cell_df_rn1.join(cell_other_df, on='CELL_ID', how='left', lsuffix='_zones', rsuffix='_other', sort=False)
cell_df_rn2 = cell_df_temp2.rename(columns={"CELL_ID_zones": "CELL_ID"})
cell_df_temp3 = cell_df_rn2.join(cell_ag_df, on='CELL_ID', how='left', lsuffix='_zones', rsuffix='_ag', sort=False)
cell_df_rn3  = cell_df_temp3.rename(columns={"CELL_ID_zones": "CELL_ID","CELL_HA_zones": "CELL_HA"})
cell_df  = cell_df_rn3.drop(columns={"CELL_ID_biophysical","CELL_ID_ag","CELL_HA_ag","PRIMARY_V7_ag"})

# Remove unesseasry files
del cell_zones_df
del cell_zones_df_rn
del cell_biophysical_df
del cell_other_df
del cell_ag_df
del cell_df_temp1
del cell_df_rn1
del cell_df_temp2
del cell_df_rn2
del cell_df_temp3

# View columns and unique values
for col in cell_df.columns:
    print(col)

dfUnique = pd.DataFrame(cell_df_sub.TENURE_NAME.unique())  
pd.set_option('display.max_rows', dfUnique.shape[0]+1)
dfUnique

# NLUM
infile = "N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif"

# Annual Native title determinations database, this is used allong with Indigenou Tenure to identify Indigenous lands
#in_cell_df_path_3 = "\\\school-les-m.shares.deakin.edu.au\\school-les-m\\Planet-A\\Data-Master\\Indigenous_lands\\cell_indigenous_df.pkl"

# Set graph colour pallette
cmap = ListedColormap(['#E0F5EB', '#34a169','white'])

################ Create some helper data and functions ########################

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
def footprint(name,index,dframe, infile):  
#def footprint(name,index,outfile,dframe):  
                    
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
 
# --- Private land mask
privateIndex  = cell_df.query("(TENURE_NAME in ['Private land'])").index

privateMask = footprint('PRIVATELAND_MASK', privateIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [NATURAL_MASK] Current natural lands
protectIndex  = cell_df.query("(VAST_CLASS in ['Bare','Residual'])").index

naturalMask = footprint('NATURAL_MASK', protectIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-ALL_MASK] Current available area for protected areas as an action
protectRemainIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Bare','Residual'])").index

naturalRemainMask = footprint('PR-ALL_MASK', protectRemainIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RESTORE_MASK] Current modified lands
restoreIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Modified','Transformed','Removed','Replaced'])").index

modifiedMask = footprint('RESTORE_MASK', restoreIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-ALL_MASK] Current available area for restoration as an action
restoreRemainIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Modified','Transformed','Removed','Replaced'])").index

modifiedRemainMask = footprint('RG-ALL_MASK', restoreRemainIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

################ PROTECTION ACTIONS 

# --- [PR-NRS] Current protected areas
nrsCapadIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 1)").index

paMask = footprint('PR-NRS_MASK', nrsCapadIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-GPA] Government protected areas
#              Heuristic: Government protected areas can occur on all land exept freehold land, and must 
PrGpaIndex = cell_df.query("(TENURE_NAME in ['Public land - other','Nature conservation areas']) and \
                            (VAST_CLASS in ['Bare','Residual']) and \
                            (PROTECTED_AREAS_2020 == 0)").index

reservePublicMask = footprint('PR-GPA_MASK', PrGpaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-PPA] Private protected areas (Protect only high quality remnant veg)
PrPpaIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Bare','Residual']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

reservePrivateMask = footprint('PR-PPA_MASK', PrPpaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-PLC] Private land conservation (Protect and more actively manage native veg)
PrPlcIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Bare','Residual','Modified']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

managePrivateMask = footprint('PR-PLC_MASK', PrPlcIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-IPA] Indigenous protected areas
#PrIpaIndex  = cell_df.query("(TENURE_NAME in ['Aboriginal land - traditional indigenous uses','Aboriginal land - other non-agricultural']) and \
#                             (VAST_CLASS in ['Bare','Residual']) and \
#                             (PROTECTED_AREAS_2020 == 0)").index

#reserveIndigenousMask = footprint('PR-IPA_MASK', PrIpaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-ICA] Indigenous collaborations
#PrIcaIndex  = cell_df.query("(NNTT_EXCLUSIVE_TITLE_2020 ==1) and \
#                             (VAST_CLASS in ['Bare','Residual','Modified']) and \
#                             (PROTECTED_AREAS_2020 == 0)").index
#
#reserveIndigCollabMask = footprint('PR-ICA_MASK', PrIcaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

################ RESTORATION ACTIONS 

# --- [RG-HIR] Human induced regeneration
RgHirIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Modified','Transformed']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

assistRegenMask = footprint('RG-HIR_MASK', RgHirIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-NRV] Environmental revegetation
RgNrvIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Removed','Replaced']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

natiRevegMask = footprint('RG-NRV_MASK', RgNrvIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-CMP] Environmental revegetation
RgCarbonMonoIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Removed','Replaced']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

monoRevegMask = footprint('RG-CMP_MASK', RgCarbonMonoIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")


## ADD RIPARIAN AS AN ACTION

############# Determine technical potential of consevration action to contribute towards conservation and sustainability ########################
#
# 

# Create lists of actions
actions = ['NATURAL_MASK','RESTORE_MASK','PR-ALL_MASK','PR-NRS_MASK','RG-ALL_MASK','PR-GPA_MASK','PR-PPA_MASK','RG-HIR_MASK','RG-CMP_MASK','RG-NRV_MASK']
conservation_actions = ['NATURAL_MASK','PR-NRS_MASK','PR-ALL_MASK','PR-GPA_MASK','PR-PPA_MASK']#,'PR-PLC_MASK','PR-IPA_MASK','PR-ICA_MASK']
restoration_actions = ['RESTORE_MASK','RG-ALL_MASK','RG-HIR_MASK','RG-NRV_MASK','RG-CMP_MASK']

# Build data frame
actions_df = pd.DataFrame({'Action type':[], 
'Indicator type':[],
'Source':[],
'Indicator':[],
'Metric':[],
'Value':[]},
  dtype="float32")

actions_df = actions_df.astype(dtype= {"Action type":"object","Indicator":"string","Indicator type":"string","Metric":"string","Source":"string"})

#### Calculate Area 
for action in actions:
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(cell_df.pivot_table(index= action, values=['CELL_HA'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Indicator',0:'Value'})
    # Filter out 0's
    pivot_temp_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_temp_rn.drop(columns={"Action status"})
    # Populate Indicator Type columns
    pivot_drop['Indicator type'] = "Area"
    # Populate Indicator Type columns
    pivot_drop['Indicator'] = "True area"
    # Populate Source data column
    pivot_drop['Source'] = "CELL_HA"
    # Populate Metric columns
    pivot_drop['Metric'] = "km2"
    # Convert hectares to km2, and round
    pivot_drop['Value'] = round((pivot_drop['Value'] / 100),0)
    # Populate Action Type columns
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)


#### Calculate overlap with suitable climate space (biodiversity indicator) ,   ~~~ Only do this on private land for avoided loss
for action in actions:
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(cell_df.pivot_table(index= action, values=['BIODIV_HIST_1990','BIODIV_SSP126_2030','BIODIV_SSP245_2030','BIODIV_SSP370_2030','BIODIV_SSP585_2030',
                                                                    'MAMMALS_HIST_1990','MAMMALS_SSP126_2030','MAMMALS_SSP245_2030','MAMMALS_SSP370_2030','MAMMALS_SSP585_2030',
                                                                    'BIRDS_HIST_1990','BIRDS_SSP126_2030','BIRDS_SSP245_2030','BIRDS_SSP370_2030','BIRDS_SSP585_2030',
                                                                    'REPTILES_HIST_1990','REPTILES_SSP126_2030','REPTILES_SSP245_2030','REPTILES_SSP370_2030','REPTILES_SSP585_2030',
                                                                    'FROGS_HIST_1990','FROGS_SSP126_2030','FROGS_SSP245_2030','FROGS_SSP370_2030','FROGS_SSP585_2030',
                                                                    'PLANTS_HIST_1990','PLANTS_SSP126_2030','PLANTS_SSP245_2030','PLANTS_SSP370_2030','PLANTS_SSP585_2030'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename colums
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
    # Filter out 0's
    pivot_temp_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_temp_rn.drop(columns={"Action status"})
    # Populate Source data column
    pivot_drop['Indicator'] = "Species climate niche"
    # Populate Metric columns
    pivot_drop['Metric'] = "%"
    # Populate Indicator Type colum
    pivot_drop['Indicator type'] = "Biodiversity"
    # Populate Action Type colum
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)

#### Calculate overlap with biome and vegetation community (ecosystem indicator)

# Create list of ecosystem data I want to report on
ecosystems = ['IBRA_SUB_NAME','IBRA_REG_NAME','NVIS_EXTANT_MVG_NAME','NVIS_EXTANT_MVS_NAME','NVIS_PRE-EURO_MVG_NAME','NVIS_PRE-EURO_MVS_NAME']

for action in actions:
    # iterate over ecosystem data
    for ecosystem in ecosystems:
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(cell_df.pivot_table(index= [action,ecosystem], values='CELL_HA', aggfunc=np.sum).stack())
        # Index to columns
        pivot.reset_index(inplace=True) 
        # Rename columns
        pivot_temp = pivot.rename(columns={action:'Action status',ecosystem:'Indicator',0:'Value'})
        # Filter out 0's
        pivot_temp_rn = pivot_temp[pivot_temp['Action status']==1]
        # Drop column
        pivot_drop  = pivot_temp_rn.drop(columns={"level_2","Action status"})
        # Convert hectares to km2, and round
        pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
        # Populate Source data column
        pivot_drop['Source'] = ecosystem
        # Populate Metric columns
        pivot_drop['Metric'] = "km2"
        # Populate Indicator Type column
        pivot_drop['Indicator type'] = "Ecosystems"
        # Populate Action Type column 
        pivot_drop['Action type'] = action
        # Reorder columns
        pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
        # Append data to actions_df
        actions_df = actions_df.append(pivot_drop)

outpath = "N://Planet-A//LUF-Modelling//Technical-potential-analysis//Conservation-actions//Data//actions_df_area-bio-eco.csv"
actions_df.to_csv(outpath,index=False, line_terminator='\n')

#### Calculate carbon emmissions abatement for conservation actions (carbon indicator, avoided loss)

# Run procedure for soil carbon
for action in conservation_actions:
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(cell_df.pivot_table(index= action, values=['SOC_T_HA_TOP_30CM'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
    # Filter out 0's
    pivot_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_rn.drop(columns={"Action status"})
    # Convert Tonnes per hectare to Tonnes per km2, and round
    pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
    # Populate Source data column
    pivot_drop['Indicator'] = "Soil organic carbon (top 30cm)"
    # Populate Metric columns
    pivot_drop['Metric'] = "Tonnes per km2"
    # Populate Indicator Type column
    pivot_drop['Indicator type'] = "Carbon"
    # Populate Action Type column
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)

# Run procedure for remnant vegetation (avoided loss)
for action in conservation_actions:
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(cell_df.pivot_table(index= action, values=['REMNANT_VEG_T_CO2_HA'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
    # Filter out 0's
    pivot_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_rn.drop(columns={"Action status"})
    # Convert Tonnes per hectare to Tonnes per km2, and round
    pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
    # Populate Source data column
    pivot_drop['Indicator'] = "Carbon abatement"
    # Populate Metric columns
    pivot_drop['Metric'] = "CO2 per km2" 
    # Populate Indicator Type column
    pivot_drop['Indicator type'] = "Carbon"
    # Populate Action Type column
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)


#### Calculate carbon emmissions sequestration for restoration actions (carbon indicator, revegetation)

# Pivot to calculate carbon/emmissions sequestration of land use types
pivot = pd.DataFrame(cell_df.pivot_table(index= 'RG-HIR_MASK', values=['HIR_BLOCK_TREES_AVG_T_CO2_HA_YR','HIR_BLOCK_DEBRIS_AVG_T_CO2_HA_YR','HIR_BLOCK_SOIL_AVG_T_CO2_HA_YR','HIR_RIP_TREES_AVG_T_CO2_HA_YR','HIR_RIP_DEBRIS_AVG_T_CO2_HA_YR','HIR_RIP_SOIL_AVG_T_CO2_HA_YR'], aggfunc=np.sum).stack()) # Need to consider the proportion of a cell under Belt planting
# Index to columns
pivot.reset_index(inplace=True) 
# Rename columns
pivot_temp = pivot.rename(columns={'RG-HIR_MASK':'Action status','level_1':'Source',0:'Value'})
# Filter out 0's
pivot_rn = pivot_temp[pivot_temp['Action status']==1]
# Drop column
pivot_drop  = pivot_rn.drop(columns={"Action status"})
# Populate Source data column
pivot_drop['Indicator'] = "Carbon sequestration"
# Convert Tonnes per hectare to Tonnes per km2, and round
pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
# Populate Metric columns
pivot_drop['Metric'] = "CO2 per km2" 
# Populate Indicator Type column
pivot_drop['Indicator type'] = "Carbon"
# Populate Action Type column
pivot_drop['Action type'] = "RG-HIR_MASK"
# Reorder columns
pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
# Append data to actions_df
actions_df = actions_df.append(pivot_drop)

# Pivot to calculate carbon/emmissions sequestration of land use types
pivot = pd.DataFrame(cell_df.pivot_table(index= 'RG-NRV_MASK', values=['EP_BLOCK_TREES_AVG_T_CO2_HA_YR','EP_BLOCK_DEBRIS_AVG_T_CO2_HA_YR','EP_BLOCK_SOIL_AVG_T_CO2_HA_YR','EP_RIP_TREES_AVG_T_CO2_HA_YR','EP_RIP_DEBRIS_AVG_T_CO2_HA_YR','EP_RIP_SOIL_AVG_T_CO2_HA_YR'], aggfunc=np.sum).stack())
# Index to columns
pivot.reset_index(inplace=True) 
# Rename columns
pivot_temp = pivot.rename(columns={'RG-NRV_MASK':'Action status','level_1':'Source',0:'Value'})
# Filter out 0's
pivot_rn = pivot_temp[pivot_temp['Action status']==1]
# Drop column
pivot_drop  = pivot_rn.drop(columns={"Action status"})
# Convert Tonnes per hectare to Tonnes per km2, and round
pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
# Populate Source data column
pivot_drop['Indicator'] = "Carbon sequestration"
# Populate Metric columns
pivot_drop['Metric'] = "CO2 per km2" 
# Populate Indicator Type column
pivot_drop['Indicator type'] = "Carbon"
# Populate Action Type column
pivot_drop['Action type'] = "RG-NRV_MASK"
# Reorder columns
pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
# Append data to actions_df
actions_df = actions_df.append(pivot_drop)

# Pivot to calculate carbon/emmissions sequestration of land use types
pivot = pd.DataFrame(cell_df.pivot_table(index= 'RG-CMP_MASK', values=['CP_BLOCK_TREES_AVG_T_CO2_HA_YR','CP_BLOCK_DEBRIS_AVG_T_CO2_HA_YR','CP_BLOCK_SOIL_AVG_T_CO2_HA_YR','CP_BELT_TREES_AVG_T_CO2_HA_YR','CP_BELT_DEBRIS_AVG_T_CO2_HA_YR','CP_BELT_SOIL_AVG_T_CO2_HA_YR'], aggfunc=np.sum).stack())
# Index to columns
pivot.reset_index(inplace=True) 
# Rename columns
pivot_temp = pivot.rename(columns={'RG-CMP_MASK':'Action status','level_1':'Source',0:'Value'})
# Filter out 0's
pivot_rn = pivot_temp[pivot_temp['Action status']==1]
# Drop column
pivot_drop  = pivot_rn.drop(columns={"Action status"})
# Convert Tonnes per hectare to Tonnes per km2, and round
pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
# Populate Source data column
pivot_drop['Indicator'] = "Carbon sequestration"
# Populate Metric columns
pivot_drop['Metric'] = "CO2 per km2" 
# Populate Indicator Type column
pivot_drop['Indicator type'] = "Carbon"
# Populate Action Type column
pivot_drop['Action type'] = "RG-CMP_MASK"
# Reorder columns
pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
# Append data to actions_df
actions_df = actions_df.append(pivot_drop)


# Add soil indicators 
#  INVEST_RKLS * C_FACTOR_VEG * P_FACTOR_AG (for veg p = 1) .... 

# Conservation
for action in conservation_actions:
    cell_df_sub = cell_df[[action,"CELL_HA", "TENURE_NAME", "INVEST_RKLS","C_FACTOR_VEG", "P_FACTOR_AG"]].copy()
    cell_df_sub['Value']  = cell_df_sub['CELL_HA'] * (cell_df_sub['INVEST_RKLS'] * cell_df_sub['C_FACTOR_VEG']* cell_df_sub['P_FACTOR_AG'])
    # Filter for only private land
    cell_df_sub = cell_df_sub[cell_df_sub['TENURE_NAME']=='Private land']
    cell_df_sub_2 = cell_df_sub[[action,"Value"]]
    # Pivot to calculate biodiversity consition of land use types, index= action,LAND TENURE... risk on land clearing on [private] free or leasehold land
    pivot = pd.DataFrame(cell_df_sub.pivot_table(index= action, values=['Value'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
    #if pivot_temp['Action status']==1:
    # Filter out 0's
    pivot_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_rn.drop(columns={"Action status"})
    # Convert Tonnes per hectare to Tonnes per km2, and round
    pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
    # Populate Source data column
    pivot_drop['Indicator'] = "Avoided soil loss"
    # Populate Source data column
    pivot_drop['Source'] = "INVEST_RKLS*C_FACTOR_VEG*P_FACTOR_AG"
    # Populate Source columns
    pivot_drop['Metric'] = "[Add metric]" 
    # Populate Indicator Type column
    pivot_drop['Indicator type'] = "Soil"
    # Populate Action Type column
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)

# Restoration

# - INVEST_RKLS = See Invest materials
# - INVEST_SDR = Sediment delivery ratio, see Invest materials
# - C_FACTOR_VEG = Teng paper.... Do in the same way for water
# - P_FACTOR_AG = Teng paper

############################################################################################################################################
# Applying P and C factors in the LUTO model
############################################################################################################################################

# C_FACTOR_VEG_PCT should be applied in cells supporting native vegetation, whereas for other land-uses C factors need to be applied 
# in the LUTO model based on land-use (from Teng et al. 2016)...  Cropping (dryland) = 7%, Cropping (irrigated) = 10%, Pasture = 8%, All forest = 3%   
# P_FACTOR_AG should be applied in areas of agricultural production and modified by regenerative agriculture methods, and be set to 1 elsewhere

# Need to calculate this in the same way as water MASK * CELL_HA * (

# Add Water indicators 
#
# - WATER_YIELD_DR_ML_HA
# - WATER_YIELD_SR_ML_HA
# Look for the difference between the 2, for conservation actions (avoided loss)  MASK * CELL_HA * (WATER_YIELD_DR_ML_HA - WATER_YIELD_SR_ML_HA)
# Look for the difference between the 2, for restoration actions MASK * CELL_HA * (WATER_YIELD_SR_ML_HA - WATER_YIELD_DR_ML_HA), reduced water yield from restoration

# Run procedure for conservation actions,   ~~~ Only do this on private land
for action in conservation_actions:
    cell_df_sub = cell_df[[action,"CELL_HA", "TENURE_NAME", "WATER_YIELD_SR_ML_HA","WATER_YIELD_DR_ML_HA"]].copy()
    cell_df_sub['Value']  = cell_df_sub['CELL_HA'] * (cell_df_sub['WATER_YIELD_DR_ML_HA'] - cell_df_sub['WATER_YIELD_SR_ML_HA'])
    # Filter for only private land
    cell_df_sub = cell_df_sub[cell_df_sub['TENURE_NAME']=='Private land']
    cell_df_sub_2 = cell_df_sub[[action,"Value"]]
    # Pivot to calculate biodiversity consition of land use types, index= action,LAND TENURE... risk on land clearing on [private] free or leasehold land
    pivot = pd.DataFrame(cell_df_sub.pivot_table(index= action, values=['Value'], aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
    #if pivot_temp['Action status']==1:
    # Filter out 0's
    pivot_rn = pivot_temp[pivot_temp['Action status']==1]
    # Drop column
    pivot_drop  = pivot_rn.drop(columns={"Action status"})
    # Convert Tonnes per hectare to Tonnes per km2, and round
    pivot_drop['Value'] = round(pivot_drop['Value'] / 100,0)
    # Populate Source data column
    pivot_drop['Indicator'] = "Avoided water loss"
    # Populate Source data column
    pivot_drop['Source'] = "WATER_YIELD_SR/DR_ML_HA/CELL_HA"
    # Populate Source columns
    pivot_drop['Metric'] = "GL" 
    # Populate Indicator Type column
    pivot_drop['Indicator type'] = "Water"
    # Populate Action Type column
    pivot_drop['Action type'] = action
    # Reorder columns
    pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop)


#### Calculate overlap with land use types (production indicator)

# Create list of ecosystem data I want to report on
landUses = ['SECONDARY_V7','LU_DESC']

for action in actions:
    # iterate over ecosystem data
    for landUse in landUses:
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(cell_df.pivot_table(index= [action,landUse], values='CELL_HA', aggfunc=np.sum).stack())
        # Index to columns
        pivot.reset_index(inplace=True) 
        # Rename columns
        pivot_rn = pivot.rename(columns={action:'Action status',landUse:'Indicator',0:'Value'})
        # Drop column
        pivot_drop_rn = pivot_rn.drop(columns={"level_2","Action status"})
        # Convert hectares to km2, and round
        pivot_drop_rn['Value'] = round(pivot_drop_rn['Value'] / 100,0)
        # Populate Source data column
        pivot_drop_rn['Source'] = landUse
        # Populate Metric columns
        pivot_drop_rn['Metric'] = "km2"
        # Populate Indicator Type column
        pivot_drop_rn['Indicator type'] = "Land use"
        # Populate Action Type column 
        pivot_drop_rn['Action type'] = action
        # Reorder columns
        pivot_drop_rn = pivot_drop_rn[['Action type','Indicator type','Source','Indicator','Metric','Value']]
        # Append data to actions_df
        actions_df = actions_df.append(pivot_drop_rn)

for action in actions:
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(cell_df.pivot_table(index= [action,'LU_DESC'], values='Yield_cell', aggfunc=np.sum).stack())
    # Index to columns
    pivot.reset_index(inplace=True) 
    # Rename columns
    pivot_rn = pivot.rename(columns={action:'Action status',landUse:'Source',0:'Value'})
    # Drop column
    pivot_drop_rn = pivot_rn.drop(columns={"level_2","Action status"})
    # Populate Source data column
    pivot_drop_rn['Indicator'] = "Yield"
    # Populate Metric columns
    pivot_drop_rn['Metric'] = "Varies depending on commodity"
    # Populate Indicator Type column
    pivot_drop_rn['Indicator type'] = "Production"
    # Populate Action Type column 
    pivot_drop_rn['Action type'] = action
    # Reorder columns
    pivot_drop_rn = pivot_drop_rn[['Action type','Indicator type','Source','Indicator','Metric','Value']]
    # Append data to actions_df
    actions_df = actions_df.append(pivot_drop_rn)





#regions = ['STATE_NAME','LGA_NAME','NRM_NAME','TENURE_NAME']

#for action in actions:
#    for region in regions:
        # Pivot to calculate biodiversity consition of land use types
#        pivot = pd.DataFrame(cell_df.pivot_table(index= [action,region], values='CELL_HA', aggfunc=np.sum).stack())
        # Index to columns
#        pivot.reset_index(inplace=True) 
#        pivot_drop  = pivot.drop(columns={"level_2"})
#        pivot_drop_rn = pivot_drop.rename(columns={action:'Action status',region:'Indicator',0:'Value'})
#        pivot_drop_rn['Indicator type'] = "Location"
#        pivot_drop_rn['Action type'] = action
#        pivot_drop_rn['Action type'] = action
#        pivot_drop_rn = pivot_drop_rn[['Action type','Action status','Indicator type','Indicator','Metric','Value']]
#        actions_df = actions_df.append(pivot_drop_rn)




# End :) 

# Supply Curve Analysis

for col in cell_df.columns:
    print(col)

cell_df_biodiv_sub = cell_df[["CELL_HA", "BIODIV_HIST_1990"]].copy()

cell_df_biodiv_sort = cell_df_biodiv_sub.sort_values("BIODIV_HIST_1990", ascending=False)
cell_df_biodiv_sort.reset_index(inplace=True) 


# Plot cumulative % species habitat against area

cell_df_biodiv_sort["BIODIV_HIST_1990_CUMSUM"] = cell_df_biodiv_sort["BIODIV_HIST_1990"].sort_index().cumsum()

#cell_df_biodiv_sort["BIODIV_HIST_1990_CUMSUM"].sort_index().plot()

plt.figure()
plt.axhline(y=30, linestyle="--", color="black", linewidth=1)
plt.axvline(x=2100000, linestyle="--", color="black", linewidth=1)
plt.plot(cell_df_biodiv_sort["BIODIV_HIST_1990_CUMSUM"].sort_index(), color="#0D8572")
plt.xlabel("Area (km2)")
plt.ylabel("Species climate space (%)")
plt.title("Species climate space supply curve (Historical)")
plt.show()





# No of data points usedN= 500  
# normal distributiondata= np.random.randn(N)  
# sort the data in ascending orderx= np.sort(data)  
# get the cdf values of yy= np.arange(N)/ float(N)  
# plottingplt.xlabel('x-axis')plt.ylabel('y-axis')  

plt.title('CDF using sorting the data')  
plt.plot(cell_df_biodiv_sort["CELL_HA"], cell_df_biodiv_sort["BIODIV_HIST_1990_CUMSUM"], marker='o')


# - Minimize food tradeoffs.
# - Minimize area tradeoffs.


