##  This is a subfile of Data_analysis.R
##  Data matching for Taiwan
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)
library(terra)
library(furrr)
library(lubridate)

##=======Read in Observations and covariates========
year = 2013
list.dirs(paste0("D:/B620_ZW/OMI_Completed/3days-rf/OMI_Completed_Results_fc"))
## Landscan population
Pop <- fread("data/Urban_Pop_0.1g_2000-2020.csv")%>%
  rename(X_Lon = x, Y_Lat = y, 
         Pop = Pop_all, UPop = Pop_U, RPop = Pop_O)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop,UPop,RPop)

target_name <- c('GridID','X_Lon','Y_Lat','Time','Year','Label','Site_NO2','ColumnNO2',
                 'u10', 'v10', 't2m', 'blh', 'tcc', 'tco3',
                 'sp', 'e', 'tcw', 'lsp', 'lai_hv', 'lai_lv','Pop',
                 'crop_ratio', 'forest_ratio','grass_ratio', 'water_ratio', 'constr_ratio')

Grid <- fread("data/Grid_Information.csv")%>%
  rast(type = "xyz",crs = "WGS84")

NO2_taiwan <- fread("data/AQI_observation_taiwan.csv", encoding = "UTF-8")%>%
  filter(ItemEngName == "NO2")

NO2_GridID <- terra::extract(Grid,NO2_taiwan%>%dplyr::select(Longitude,Latitude),ID= FALSE)

## match on 0.1 degree cells (calculated mean NO2)
Grid_NO2_obs <- NO2_taiwan%>%cbind(NO2_GridID)%>%
  group_by(GridID,Year,DOY)%>%
  summarise(value = mean(value)*1.88)%>%ungroup()%>%
  rename(Site_NO2 = value)%>%
  as.data.frame()%>%
  mutate(Date = strptime(paste(Year, DOY), format="%Y %j")%>%as.Date(),
         Month = month(Date))

## create label
Grid_NO2_label <- Grid_NO2_obs%>%
  group_by(GridID, Month)%>%
  summarise(threshold = mean(Site_NO2) + 2*sd(Site_NO2))%>%
  ungroup()%>%left_join(Grid_NO2_obs,.)%>%
  mutate(Label = ifelse(Site_NO2>threshold, 1, 0))%>%
  dplyr::select(-threshold)%>%
  na.omit()

NO2_GridID_distinct <- distinct(Grid_NO2_label%>%dplyr::select(GridID))

ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"
OMI_path <- "D:/B620_ZW/OMI_Completed/2-ground_est/result/OMI_prediction_csv/"

for(year in 2007:2020){
  cat("=========Processing Year: ",year,"========\n")
  
  ## New OMI data
  plan(multisession, workers = 2)
  OMI <- list.files(paste0(OMI_path,year),full.names = TRUE)%>%
    future_map_dfr(function(file){
      doy <- substr(file,96,98)%>%as.numeric()

      OMI_daily <- fread(file)%>%
        rename(ColumnNO2_Pre = ColumnNO2)%>%
        mutate(ColumnNO2 = ifelse(!is.na(ColumnNO2_True),
                                  ColumnNO2_True,ColumnNO2_Pre))%>%
        mutate(DOY = doy)%>%
        dplyr::select(GridID,ColumnNO2,DOY)
      
      NO2_GridID_distinct%>%
        left_join(OMI_daily)
    })
  future:::ClusterRegistry("stop")
  
  ## ERA5 data
  plan(multisession, workers = 2)
  ERA5 <- list.files(paste0(ERA5_path,year),full.names = TRUE)%>%
    future_map_dfr(function(file){
      doy <- yday(substr(file,70,79))
      ERA5_daily <- fread(file)%>%mutate(DOY = doy)%>%
        dplyr::select(-X_Lon,-Y_Lat)
      
      NO2_GridID_distinct%>%
        left_join(ERA5_daily)
    })
  future:::ClusterRegistry("stop")
  
  ## Landuse classification
  lucc <- fread(paste0("D:/B620_ZW/Landuse/full/lucc_full_",year,".csv"))%>%
    dplyr::select(c('GridID',"X_Lon","Y_Lat","crop_ratio",
                    "forest_ratio","grass_ratio",
                    "water_ratio","constr_ratio"))

  ## Data matching
  Taiwan_data <- Grid_NO2_label%>%
    filter(Year == year)%>%
    list(OMI, ERA5, lucc, Pop%>%filter(Year == year))%>%
    reduce(left_join)%>%
    mutate(Date = strptime(paste(Year, DOY), format="%Y %j")%>%as.Date(),
           Time = difftime(Date,as.Date("2013-01-01"), units = "days")%>%as.numeric())%>%
    dplyr::select(all_of(target_name))

  fwrite(Taiwan_data%>%na.omit, paste0("result/Data_analysis/Training_dataset/Taiwan_NO2_0.1g_",year,".csv"))
}
