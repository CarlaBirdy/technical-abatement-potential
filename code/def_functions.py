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


################ Create some helper data and functions ########################

# Convert 1D column to 2D spatial array
def conv_1D_to_2D(in_1D_array,infile):
    """Convert 1D column to 2D spatial array using the NLUM mask."""
    # Open NLUM_ID as mask raster and get metadata
    with rasterio.open(infile) as rst:
        
        # Load a 2D masked array with nodata masked out
        NLUM_ID_raster = rst.read(1, masked=True) 
        
        # Create a 0/1 binary mask
        NLUM_mask = NLUM_ID_raster.mask == False

        # Set up some data structures to enable conversion on 1D arrays to 2D
        array_2D = np.zeros(NLUM_ID_raster.shape)
        xy = np.nonzero(NLUM_mask)
        
    array_2D[xy] = np.array(in_1D_array)
    array_2D[array_2D == 0] = np.nan
    return array_2D

# Convert 1D column to 2D spatial array and plot map
def map_in_2D(col, in_data): # data = 'continuous' or 'categorical'
    """Convert 1D column to 2D spatial array and plot map."""
    a2D = conv_1D_to_2D(col)
    if in_data == 'categorical':
        n = col.nunique()
        cmap = matplotlib.colors.ListedColormap(np.random.rand(n,3))
        plt.imshow(a2D, cmap=cmap, resample=False)
    elif in_data == 'continuous':
        plt.imshow(a2D, cmap='pink', resample=False)
    plt.show()

# Convert object columns to categories and downcast int64 columns to save memory and space
def downcast(dframe):
    """Convert object columns to categories and downcast int64 columns to save memory and space."""
    obj_cols = dframe.select_dtypes(include = ["object"]).columns
    dframe[obj_cols] = dframe[obj_cols].astype('category')
    int64_cols = dframe.select_dtypes(include = ["int64"]).columns
    dframe[int64_cols] = dframe[int64_cols].apply(pd.to_numeric, downcast = 'integer')
    
# Process land use sieve mask                   
def footprint(name,index,dframe,infile):  
    """Appends mask to `cell_df` and returns map of the spatial footprint of an action."""                    
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
 
        # Set graph colour pallette
        cmap = ListedColormap(['#E0F5EB', '#34a169','white'])

        # Plot GeoTiff
        plt.imshow(np2dArray, cmap=cmap, vmin=0, vmax=2)
        
        # Save the output file as GeoTiff
        # r'N:\Planet-A\LUF-Modelling\Land-Use-Sieve\Data\Masks\SOL_Conservation_Government-protected-area_Potential-area-mask_1km.tif'
        #with rasterio.open(outfile, "w+", **meta) as dst:
        #    dst.write(np2dArray.astype(np.uint8), indexes = 1)
            

# Returns a tabular report `out_df` of the full technical potential impact of an `action`                
def technical_potential(action,in_df): 
    """Returns a tabular report `out_df` of the full technical potential impact of an `action` to conserve biodiversity and ecosystems, contribute towards sustainability outcomes and impact agricultural production using `in_data`."""
  
    # Build data frame
    out_df = pd.DataFrame({'Action type':[], 
    'Indicator type':[],
    'Source':[],
    'Indicator':[],
    'Metric':[],
    'Value':[]},
      dtype="float32")
    
    out_df = out_df.astype(dtype= {"Action type":"object","Indicator":"string","Indicator type":"string","Metric":"string","Source":"string"})

    # Extract action type, RG for restoration and PR for protection
    actionType = action[:2]

# --- AREA
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(in_df.pivot_table(index= action, values=['CELL_HA'], aggfunc=np.sum).stack())
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
    out_df = out_df.append(pivot_drop)

# --- BIODIVERISTY
    #### Calculate overlap with suitable climate space (biodiversity indicator) ,   ~~~ Only do this on private land for avoided loss

    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(in_df.pivot_table(index= action, values=['BIODIV_HIST_1990','BIODIV_SSP126_2030','BIODIV_SSP245_2030','BIODIV_SSP370_2030','BIODIV_SSP585_2030',
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
    out_df = out_df.append(pivot_drop)

# --- ECOSYSTEMS    
    #### Calculate overlap with biome and vegetation community (ecosystem indicator)

    # Create list of ecosystem data I want to report on
    ecosystems = ['IBRA_SUB_NAME','IBRA_REG_NAME','NVIS_EXTANT_MVG_NAME','NVIS_EXTANT_MVS_NAME','NVIS_PRE-EURO_MVG_NAME','NVIS_PRE-EURO_MVS_NAME']

    # iterate over ecosystem data
    for ecosystem in ecosystems:
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(in_df.pivot_table(index= [action,ecosystem], values='CELL_HA', aggfunc=np.sum).stack())
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
        out_df = out_df.append(pivot_drop)

# --- CARBON
#### Calculate carbon emmissions abatement for conservation actions (carbon indicator, avoided loss)
    if actionType == "PR":
    
        #--- SOIL CARBON
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(in_df.pivot_table(index= action, values=['SOC_T_HA_TOP_30CM'], aggfunc=np.sum).stack())
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
        out_df = out_df.append(pivot_drop)
     
        #--- AVOIDED LOSS (of remnant vegetation)
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(in_df.pivot_table(index= action, values=['REMNANT_VEG_T_CO2_HA'], aggfunc=np.sum).stack())
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
        out_df = out_df.append(pivot_drop)

    if actionType == "RG":
        
        restorationType = action[3:-5]
        
        if restorationType == 'NRV':
            fulcamType = 'EP'
            
        if restorationType == 'HIR':
            fulcamType = 'HIR'
        
        if restorationType == 'CMP':
            fulcamType = 'CP'
            
        #### CARBON SQUESTRATION (Through revegetation)
        
        if restorationType == 'CMP':
            columnNames = [fulcamType+'_BLOCK_TREES_AVG_T_CO2_HA_YR',fulcamType+'_BLOCK_DEBRIS_AVG_T_CO2_HA_YR',fulcamType+'_BLOCK_SOIL_AVG_T_CO2_HA_YR',fulcamType+'_BELT_TREES_AVG_T_CO2_HA_YR',fulcamType+'_BELT_DEBRIS_AVG_T_CO2_HA_YR',fulcamType+'_BELT_SOIL_AVG_T_CO2_HA_YR']
        else:
            columnNames = [fulcamType+'_BLOCK_TREES_AVG_T_CO2_HA_YR',fulcamType+'_BLOCK_DEBRIS_AVG_T_CO2_HA_YR',fulcamType+'_BLOCK_SOIL_AVG_T_CO2_HA_YR',fulcamType+'_RIP_TREES_AVG_T_CO2_HA_YR',fulcamType+'_RIP_DEBRIS_AVG_T_CO2_HA_YR',fulcamType+'_RIP_SOIL_AVG_T_CO2_HA_YR']
        
        # Pivot to calculate carbon/emmissions sequestration of land use types
        pivot = pd.DataFrame(in_df.pivot_table(index= action, values=columnNames, aggfunc=np.sum).stack())
        # Index to columns
        pivot.reset_index(inplace=True) 
        # Rename columns
        pivot_temp = pivot.rename(columns={action:'Action status','level_1':'Source',0:'Value'})
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
        pivot_drop['Action type'] = action
        # Reorder columns
        pivot_drop = pivot_drop[['Action type','Indicator type','Source','Indicator','Metric','Value']]
        # Append data to actions_df
        out_df = out_df.append(pivot_drop)


# --- WATER 

    #### AVOIDED WATER LOSS
    # Run procedure for conservation actions
    if actionType == "PR":
        cell_df_sub = in_df[[action,"CELL_HA", "TENURE_NAME", "WATER_YIELD_SR_ML_HA","WATER_YIELD_DR_ML_HA"]].copy()
        cell_df_sub['Value']  = cell_df_sub['CELL_HA'] * (cell_df_sub['WATER_YIELD_DR_ML_HA'] - cell_df_sub['WATER_YIELD_SR_ML_HA'])
        # Filter for only private land
        cell_df_sub = cell_df_sub[cell_df_sub['TENURE_NAME']=='Private land'] # Only do this on private land
        cell_df_sub_2 = cell_df_sub[[action,"Value"]]
        # Pivot to calculate biodiversity consition of land use types, index= action,LAND TENURE... risk on land clearing on [private] free or leasehold land
        pivot = pd.DataFrame(cell_df_sub_2.pivot_table(index= action, values=['Value'], aggfunc=np.sum).stack())
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
        out_df = out_df.append(pivot_drop)

    #### POTENTIAL WATER LOSS
    # Run procedure for restoration actions
    if actionType == "RG":
        cell_df_sub = in_df[[action,"CELL_HA", "TENURE_NAME", "WATER_YIELD_SR_ML_HA","WATER_YIELD_DR_ML_HA"]].copy()
        cell_df_sub['Value']  = cell_df_sub['CELL_HA'] * (cell_df_sub['WATER_YIELD_SR_ML_HA'] - cell_df_sub['WATER_YIELD_DR_ML_HA'])
        # Filter for only private land
        cell_df_sub = cell_df_sub[cell_df_sub['TENURE_NAME']=='Private land']
        cell_df_sub_2 = cell_df_sub[[action,"Value"]]
        # Pivot to calculate biodiversity consition of land use types, index= action,LAND TENURE... risk on land clearing on [private] free or leasehold land
        pivot = pd.DataFrame(cell_df_sub_2.pivot_table(index= action, values=['Value'], aggfunc=np.sum).stack())
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
        pivot_drop['Indicator'] = "Potential water loss"
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
        out_df = out_df.append(pivot_drop)


# --- LAND USE
     
    #### LAND USE TYPES
    #### Calculate overlap with land use types (production indicator)
    
    # Create list of ecosystem data I want to report on
    landUses = ['SECONDARY_V7','LU_DESC']

    # iterate over land use data
    for landUse in landUses:
        # Pivot to calculate biodiversity consition of land use types
        pivot = pd.DataFrame(in_df.pivot_table(index= [action,landUse], values='CELL_HA', aggfunc=np.sum).stack())
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
        out_df = out_df.append(pivot_drop_rn)
   
    #### AGRICULTURAL PRODUCTION
    # Pivot to calculate biodiversity consition of land use types
    pivot = pd.DataFrame(in_df.pivot_table(index= [action,'LU_DESC'], values='Yield_cell', aggfunc=np.sum).stack())
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
    out_df = out_df.append(pivot_drop_rn)


# --- SOILS
# Add soil indicators 
#  INVEST_RKLS * C_FACTOR_VEG * P_FACTOR_AG (for veg p = 1) .... 

# Conservation
    if actionType == "PR":
        cell_df_sub = in_df[[action,"CELL_HA", "TENURE_NAME", "INVEST_RKLS","C_FACTOR_VEG", "P_FACTOR_AG"]].copy()
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
        out_df = out_df.append(pivot_drop)

# Restoration
    if actionType == "RG":
        cell_df_sub = in_df[[action,"CELL_HA", "TENURE_NAME", "INVEST_RKLS","C_FACTOR_VEG", "P_FACTOR_AG"]].copy()
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
        out_df = out_df.append(pivot_drop)


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

# --- LOCATION
#regions = ['STATE_NAME','LGA_NAME','NRM_NAME','TENURE_NAME']

#for action in actions:
#    for region in regions:
        # Pivot to calculate biodiversity consition of land use types
#        pivot = pd.DataFrame(in_df.pivot_table(index= [action,region], values='CELL_HA', aggfunc=np.sum).stack())
        # Index to columns
#        pivot.reset_index(inplace=True) 
#        pivot_drop  = pivot.drop(columns={"level_2"})
#        pivot_drop_rn = pivot_drop.rename(columns={action:'Action status',region:'Indicator',0:'Value'})
#        pivot_drop_rn['Indicator type'] = "Location"
#        pivot_drop_rn['Action type'] = action
#        pivot_drop_rn['Action type'] = action
#        pivot_drop_rn = pivot_drop_rn[['Action type','Action status','Indicator type','Indicator','Metric','Value']]
#        actions_df = actions_df.append(pivot_drop_rn)


    return out_df

# End :) 
