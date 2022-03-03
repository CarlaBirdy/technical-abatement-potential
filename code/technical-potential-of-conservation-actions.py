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

# --- Load python modules
import numpy as np
import pandas as pd
import rasterio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as npimport 
import pandas as pd 
import os
os.chdir("N://Planet-A//LUF-Modelling//Technical-potential-analysis//Conservation-actions//technical-abatement-potential//code//")
import def_functions

# --- Input files

# NLUM
infile = "N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif"

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

# --- Remove unesseasry files
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

# --- Inspect data
for col in cell_df.columns:
    print(col)

dfUnique = pd.DataFrame(cell_df_sub.TENURE_NAME.unique())  
pd.set_option('display.max_rows', dfUnique.shape[0]+1)
dfUnique
        
################ Determine spatial footprint of consevration action ########################
 
# --- Private land mask
privateIndex  = cell_df.query("(TENURE_NAME in ['Private land'])").index

privateMask = def_functions.footprint('PRIVATELAND_MASK', privateIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [NATURAL_MASK] Current natural lands
protectIndex  = cell_df.query("(VAST_CLASS in ['Bare','Residual'])").index

naturalMask = def_functions.footprint('NATURAL_MASK', protectIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-ALL_MASK] Current available area for protected areas as an action
protectRemainIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Bare','Residual'])").index

naturalRemainMask = def_functions.footprint('PR-ALL_MASK', protectRemainIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RESTORE_MASK] Current modified lands
restoreIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Modified','Transformed','Removed','Replaced'])").index

modifiedMask = footprint('RESTORE_MASK', restoreIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-ALL_MASK] Current available area for restoration as an action
restoreRemainIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 0) and \
                              (VAST_CLASS in ['Modified','Transformed','Removed','Replaced'])").index

modifiedRemainMask = def_functions.footprint('RG-ALL_MASK', restoreRemainIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

### PROTECTION ACTIONS 

# --- [PR-NRS] Current protected areas
nrsCapadIndex  = cell_df.query("(PROTECTED_AREAS_2020 == 1)").index

paMask = def_functions.footprint('PR-NRS_MASK', nrsCapadIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-GPA] Government protected areas
#              Heuristic: Government protected areas can occur on all land exept freehold land, and must 
PrGpaIndex = cell_df.query("(TENURE_NAME in ['Public land - other','Nature conservation areas']) and \
                            (VAST_CLASS in ['Bare','Residual']) and \
                            (PROTECTED_AREAS_2020 == 0)").index

reservePublicMask = def_functions.footprint('PR-GPA_MASK', PrGpaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-PPA] Private protected areas (Protect only high quality remnant veg)
PrPpaIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Bare','Residual']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

reservePrivateMask = def_functions.footprint('PR-PPA_MASK', PrPpaIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [PR-PLC] Private land conservation (Protect and more actively manage native veg)
PrPlcIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Bare','Residual','Modified']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

managePrivateMask = def_functions.footprint('PR-PLC_MASK', PrPlcIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

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

### RESTORATION ACTIONS 

# --- [RG-HIR] Human induced regeneration
RgHirIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Modified','Transformed']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

assistRegenMask = def_functions.footprint('RG-HIR_MASK', RgHirIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-NRV] Environmental revegetation
RgNrvIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Removed','Replaced']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

natiRevegMask = def_functions.footprint('RG-NRV_MASK', RgNrvIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")

# --- [RG-CMP] Environmental revegetation
RgCarbonMonoIndex  = cell_df.query("(TENURE_NAME in ['Private land']) and \
                             (VAST_CLASS in ['Removed','Replaced']) and \
                             (PROTECTED_AREAS_2020 == 0)").index

monoRevegMask = def_functions.footprint('RG-CMP_MASK', RgCarbonMonoIndex, cell_df,"N://Planet-A//Data-Master//National_Landuse_Map//NLUM_2010-11_mask.tif")


## ADD RIPARIAN AS AN ACTION

############# Determine technical potential of consevration action to contribute towards conservation and sustainability ########################

# --- List of actions
actions = ['PR-ALL_MASK','PR-NRS_MASK','PR-GPA_MASK','PR-PPA_MASK','RG-HIR_MASK','RG-CMP_MASK','RG-NRV_MASK']

# --- Run procedure 

for action in actions:
    print(action)
    out_df = def_functions.technical_potential(action,cell_df)
    # Append data to actions_df
    out_df = out_df.append(out_df)
    len(out_df)

# --- Save out file
outpath = "N://Planet-A//LUF-Modelling//Technical-potential-analysis//Conservation-actions//Data//actions_df_03032022.csv"
out_df.to_csv(outpath,index=False, line_terminator='\n')


############# Supply Curve Analysis #############

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

# End :) 

