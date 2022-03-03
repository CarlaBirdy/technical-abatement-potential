#
# summaise-technical-potentials.R
#
# This script summarises the "technical abatement/mitigation potential"  
# of different conservation actions to meet conservation and sustainability targets.
# 
# By Carla Arcihbald (c.archibaldc@deakin.edu.au)
# Created: 2022-02-02
# Last modified: 2022-02-02
#
#
################################################

################ Set up ########################
# Load packages
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)

# Set wd
wd <- "N:/Planet-A/LUF-Modelling/Technical-potential-analysis/Conservation-actions"

# Read in files
actions_df <- read_csv(paste0(wd, "/Data/actions_df_area-bio-eco.csv"))


View(actions_df %>% 
       dplyr::filter(`Indicator type` %in% "Biodiversity"))

# Plotting Area #####

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK"))  %>%
  dplyr::filter(`Indicator type` %in% "Area")  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  ggplot(., aes(`Action type`, CELL_HA,fill = `Action type`)) + 
  geom_col() +
  labs(title = "Area (km2) of each action type")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+#,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plotting Biodiversity #####

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
            dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
            dplyr::filter(grepl("BIODIV", `Source`, fixed = TRUE))  %>%
            #dplyr::filter(`Action type` %in% c("PR-GPA_MASK","PR-PPA_MASK","RG-NRV_MASK","RG-HIR_MASK"))  %>%
            dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
            tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
            dplyr::mutate(BIODIV_SSP126_2020 = BIODIV_HIST_1990)%>%
            dplyr::mutate(BIODIV_SSP245_2020 = BIODIV_HIST_1990)%>%
            dplyr::mutate(BIODIV_SSP370_2020 = BIODIV_HIST_1990)%>%
            dplyr::mutate(BIODIV_SSP585_2020 = BIODIV_HIST_1990)%>%
            #dplyr::select(-BIODIV_HIST_1990)%>%
            tidyr::pivot_longer(!c(`Action type`,BIODIV_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
            dplyr::mutate(SSP = Source) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"BIODIV_")) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
            dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
            dplyr::filter(SSP %in% "SSP245")  %>%
            dplyr::mutate(Year = Source) %>%
            dplyr::mutate(Year = str_remove_all(Year,"BIODIV_")) %>%
            dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
            dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
            dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
            dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
            dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
            dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
            dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
            dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
            dplyr::select(-c(Source))%>%
            ggplot(., aes(Year, Value, fill = `Action type`)) + 
            geom_col(position="dodge") +
            labs(title = "Percentage of all species climate space represented")+
            xlab("Year")+
            ylab("% of species climate space")+
            scale_fill_viridis_d()+
            facet_wrap(~Domain,scale="free_x")+
            theme_bw()
                  
# MAMMALS 

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
  dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
  dplyr::filter(grepl("MAMMALS", `Source`, fixed = TRUE))  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(MAMMALS_SSP126_2020 = MAMMALS_HIST_1990)%>%
  dplyr::mutate(MAMMALS_SSP245_2020 = MAMMALS_HIST_1990)%>%
  dplyr::mutate(MAMMALS_SSP370_2020 = MAMMALS_HIST_1990)%>%
  dplyr::mutate(MAMMALS_SSP585_2020 = MAMMALS_HIST_1990)%>%
  tidyr::pivot_longer(!c(`Action type`,MAMMALS_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
  dplyr::mutate(SSP = Source) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"MAMMALS_")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
  dplyr::filter(SSP %in% "SSP245")  %>%
  dplyr::mutate(Year = Source) %>%
  dplyr::mutate(Year = str_remove_all(Year,"MAMMALS_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source))%>%
  ggplot(., aes(Year, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Percentage of mammal species climate space represented")+
  xlab("Year")+
  ylab("% of species climate space")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+
  theme_bw()

# REPTILES

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
  dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
  dplyr::filter(grepl("REPTILES", `Source`, fixed = TRUE))  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(REPTILES_SSP126_2020 = REPTILES_HIST_1990)%>%
  dplyr::mutate(REPTILES_SSP245_2020 = REPTILES_HIST_1990)%>%
  dplyr::mutate(REPTILES_SSP370_2020 = REPTILES_HIST_1990)%>%
  dplyr::mutate(REPTILES_SSP585_2020 = REPTILES_HIST_1990)%>%
  tidyr::pivot_longer(!c(`Action type`,REPTILES_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
  dplyr::mutate(SSP = Source) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"REPTILES_")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
  dplyr::filter(SSP %in% "SSP245")  %>%
  dplyr::mutate(Year = Source) %>%
  dplyr::mutate(Year = str_remove_all(Year,"REPTILES_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source))%>%
  ggplot(., aes(Year, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Percentage of reptile species climate space represented")+
  xlab("Year")+
  ylab("% of species climate space")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+
  theme_bw()

# BIRDS

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
  dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
  dplyr::filter(grepl("BIRDS", `Source`, fixed = TRUE))  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(BIRDS_SSP126_2020 = BIRDS_HIST_1990)%>%
  dplyr::mutate(BIRDS_SSP245_2020 = BIRDS_HIST_1990)%>%
  dplyr::mutate(BIRDS_SSP370_2020 = BIRDS_HIST_1990)%>%
  dplyr::mutate(BIRDS_SSP585_2020 = BIRDS_HIST_1990)%>%
  tidyr::pivot_longer(!c(`Action type`,BIRDS_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
  dplyr::mutate(SSP = Source) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"BIRDS_")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
  dplyr::filter(SSP %in% "SSP245")  %>%
  dplyr::mutate(Year = Source) %>%
  dplyr::mutate(Year = str_remove_all(Year,"BIRDS_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source))%>%
  ggplot(., aes(Year, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Percentage of bird species climate space represented")+
  xlab("Year")+
  ylab("% of species climate space")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+
  theme_bw()

# FROGS

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
  dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
  dplyr::filter(grepl("FROGS", `Source`, fixed = TRUE))  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(FROGS_SSP126_2020 = FROGS_HIST_1990)%>%
  dplyr::mutate(FROGS_SSP245_2020 = FROGS_HIST_1990)%>%
  dplyr::mutate(FROGS_SSP370_2020 = FROGS_HIST_1990)%>%
  dplyr::mutate(FROGS_SSP585_2020 = FROGS_HIST_1990)%>%
  tidyr::pivot_longer(!c(`Action type`,FROGS_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
  dplyr::mutate(SSP = Source) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"FROGS_")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
  dplyr::filter(SSP %in% "SSP245")  %>%
  dplyr::mutate(Year = Source) %>%
  dplyr::mutate(Year = str_remove_all(Year,"FROGS_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source))%>%
  ggplot(., aes(Year, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Percentage of amphibian species climate space represented")+
  xlab("Year")+
  ylab("% of species climate space")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+
  theme_bw()


# PLANTS

actions_df %>% 
  dplyr::filter(!`Action type` %in% c("NATURAL_MASK","RESTORE_MASK")) %>%
  dplyr::filter(`Indicator type` %in% "Biodiversity")  %>%
  dplyr::filter(grepl("PLANTS", `Source`, fixed = TRUE))  %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`))%>%
  tidyr::pivot_wider(names_from = Source, values_from = Value)%>%
  dplyr::mutate(PLANTS_SSP126_2020 = PLANTS_HIST_1990)%>%
  dplyr::mutate(PLANTS_SSP245_2020 = PLANTS_HIST_1990)%>%
  dplyr::mutate(PLANTS_SSP370_2020 = PLANTS_HIST_1990)%>%
  dplyr::mutate(PLANTS_SSP585_2020 = PLANTS_HIST_1990)%>%
  tidyr::pivot_longer(!c(`Action type`,PLANTS_HIST_1990), names_to = 'Source', values_to = 'Value')%>%
  dplyr::mutate(SSP = Source) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"PLANTS_")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_1990")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2030")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2050")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2070")) %>%
  dplyr::mutate(SSP = str_remove_all(SSP,"_2090")) %>%
  dplyr::filter(SSP %in% "SSP245")  %>%
  dplyr::mutate(Year = Source) %>%
  dplyr::mutate(Year = str_remove_all(Year,"PLANTS_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"HIST_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP126_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP245_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP370_")) %>%
  dplyr::mutate(Year = str_remove_all(Year,"SSP585_")) %>%
  dplyr::mutate(Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "PR", "Protection")) %>%
  dplyr::mutate(Domain = str_replace_all(Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source))%>%
  ggplot(., aes(Year, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Percentage of plant species climate space represented")+
  xlab("Year")+
  ylab("% of species climate space")+
  scale_fill_viridis_d()+
  facet_wrap(~Domain,scale="free_x")+
  theme_bw()

# Plotting Ecosystems #####

ecosystem_groups <- read_csv("N:/Planet-A/LUF-Modelling/Technical-potential-analysis/Conservation-actions/Data/MVG_grouping.csv") %>%
                    rename(Indicator = `Major Vegetation Group`)

action_eco_t1 <- actions_df %>% 
  dplyr::filter(`Source` %in% "NVIS_PRE-EURO_MVG_NAME") %>%
  dplyr::filter(`Action type` %in% c("RG-ALL_MASK","PR-ALL_MASK","PR-GPA_MASK","PR-PPA_MASK","RG-NRV_MASK","RG-HIR_MASK","PR-NRS_MASK")) 
  
action_eco_t2 <-  merge(action_eco_t1,ecosystem_groups,by="Indicator")%>%
                  dplyr::select(-"MVG_ID")
# Rainforest
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Rainforests") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Rainforests represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Eucalyptus Forests and Woodlands
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Eucalyptus Forests and Woodlands") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Eucalyptus Forests and Woodlands represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Other Forests and Woodlands
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Other Forests and Woodlands") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Other Forests and Woodlands represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Open Woodlands
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Open Woodlands") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Open Woodlands represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Shrublands and Grasslands
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Shrublands and Grasslands") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Shrublands and Grasslands represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Other
action_eco_t2 %>%
  dplyr::select(-c(Indicator, Metric,`Indicator type`)) %>%
  dplyr::mutate(Action_Domain = str_sub(`Action type`, start = 1, end = 2)) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "PR", "Protection")) %>%
  dplyr::mutate(Action_Domain = str_replace_all(Action_Domain, "RG", "Restoration")) %>%
  dplyr::select(-c(Source)) %>% 
  dplyr::filter(`Domain` %in% "Other") %>%
  ggplot(., aes(`Action type`, Value, fill = `Action type`)) + 
  geom_col(position="dodge") +
  labs(title = "Other habitats represented by actions")+
  ylab("Area (km2)")+
  scale_fill_viridis_d()+
  facet_wrap(~Action_Domain,scale="free_x")+ #,scale="free"
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())





















