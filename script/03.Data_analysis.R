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
library(geosphere)

##=======Read in Observations and covariates========
year = 2013
list.dirs(paste0("D:/B620_ZW/OMI_Completed/3days-rf/OMI_Completed_Results_fc"))
## Landscan population
Pop <- fread("data/Urban_Pop_0.1g_2000-2020.csv")%>%
  rename(X_Lon = x, Y_Lat = y, 
         Pop = Pop_all, UPop = Pop_U, RPop = Pop_O)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop,UPop,RPop)%>%
  left_join(fread("data/Grid_Information.csv")%>%
              dplyr::select(X_Lon,Y_Lat,GridID))%>%
  na.omit()

target_name <- c('GridID','X_Lon','Y_Lat','Time','Year','Label','Site_NO2','ColumnNO2',
                 'u10', 'v10', 't2m', 'blh', 'tcc', 'tco3',
                 'sp', 'e', 'tcw', 'lsp', 'lai_hv', 'lai_lv','Pop',
                 'crop_ratio', 'forest_ratio','grass_ratio', 'water_ratio', 'constr_ratio')

# ## Grid inoformation
# fread(paste0("D:/B620_ZW/Landuse/full/lucc_full_2013.csv"))%>%
#   dplyr::select(X_Lon, Y_Lat, GridID)%>%
#   rowwise()%>%
#   mutate(DC = distHaversine(c(X_Lon,Y_Lat),c(104.15,35.85)),
#          D1 = distHaversine(c(X_Lon,Y_Lat),c(73.5,18.1)),
#          D2 = distHaversine(c(X_Lon,Y_Lat),c(73.5,53.6)),
#          D3 = distHaversine(c(X_Lon,Y_Lat),c(134.8,18.1)),
#          D4 = distHaversine(c(X_Lon,Y_Lat),c(134.8,53.6)))%>%
#   fwrite("data/Grid_Information.csv")

Grid <- fread("data/Grid_Information.csv")%>%
  rast(type = "xyz",crs = "WGS84")

ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"
OMI_path <- "D:/B620_ZW/OMI_Completed/2-ground_est/result/OMI_prediction_csv/"

for(year in 2013:2020){
  cat("=========Processing Year: ",year,"========\n")
  # summary(fread(paste0("D:/B620_ZW/China_Site_based/NO2/NO2_Average/china_sites_NO2_2014.csv")))
  NO2_obs <- paste0("data/AQI_obs/china_sites_daily_",year,".csv")%>%
    fread()%>%filter(type == "NO2")%>%
    dplyr::select(Lon, Lat, value, date, count)%>%
    filter(count >= 18, value >0)%>% ## filter hourly observations > 18 in one day
    mutate(value = ifelse(date<20180900, value*0.92,value))

  NO2_GridID <- terra::extract(Grid,NO2_obs%>%dplyr::select(Lon,Lat),ID= FALSE)
  
  ## match on 0.1 degree cells (calculated mean NO2)
  Grid_NO2_obs <- NO2_obs%>%cbind(NO2_GridID)%>%
    group_by(GridID,date)%>%
    summarise(value = mean(value))%>%ungroup()%>%
    mutate(date = date%>%str_c()%>%as_date(),
           DOY = yday(date), Month = month(date))%>%
    rename(Site_NO2 = value)%>%
    as.data.frame()
  
  ## create label
  Grid_NO2_label <- Grid_NO2_obs%>%
    group_by(GridID, Month)%>%
    summarise(threshold = mean(Site_NO2) + 2*sd(Site_NO2))%>%
    ungroup()%>%left_join(Grid_NO2_obs,.)%>%
    mutate(Label = ifelse(Site_NO2>threshold, 1, 0))%>%
    dplyr::select(-threshold)%>%
    na.omit()
  
  NO2_GridID_distinct <- distinct(Grid_NO2_label%>%dplyr::select(GridID))
  
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

  ## Modify time variable to indicate day distance between 2013-01-01
  ML_data <- list(Grid_NO2_label, OMI, ERA5, lucc, Pop%>%filter(Year == year))%>%
    reduce(left_join)%>%
    mutate(Date = strptime(paste(Year, DOY), format="%Y %j")%>%as.Date(),
           Time = difftime(Date,as.Date("2013-01-01"), units = "days")%>%as.numeric())%>%
    dplyr::select(all_of(target_name))

  ## read in Taiwan data
  TW_data <- fread(paste0("result/Data_analysis/Training_dataset/Taiwan_NO2_0.1g_",year,".csv"))
  
  Train_data <- full_join(ML_data, TW_data)

  fwrite(Train_data%>%na.omit, paste0("result/Data_analysis/Training_dataset/Training_NO2_0.1g_",year,".csv"))
}

## Please run python smote script first to simulate data
##========Data after resampling========
NO2_smote <- 2013:2020%>%
  map_dfr(~ fread(paste0("Train_Raw/Resampled_Train_Raw_",.x,".csv")))

NO2_raw <- 2013:2020%>%
  map_dfr(~ fread(paste0("result/Data_analysis/Training_dataset/Training_NO2_0.1g_",.x,".csv"))%>%
            dplyr::select(GridID,Time,Year,Site_NO2))%>%
  mutate(Raw_label = 1)

NO2_resample <- NO2_smote%>%na.omit()%>%
  left_join(NO2_raw)%>%
  mutate(Raw_label = replace_na(Raw_label, 0))
NO2_resample%>%colnames

fwrite(NO2_resample,paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))


##========Observations marked in specific region========
source("function/Spatial_analysis.R")

polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

SR1 <- st_read(dsn = "data/specific_region/jingjinji", layer = "jingjinji")
SR2 <- st_read(dsn = "data/specific_region/changsanjiao", layer = "changsanjiao")
SR3 <- st_read(dsn = "data/specific_region/zhusanjiao", layer = "zhusanjiao")


X_train <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time))%>%
  filter(Label == 0)%>%
  dplyr::select(all_of(target_name))

## In specific region
X_train_SR <- list(SR1,SR2,SR3)%>%
  map2_dfr(.x = ., .y = c("BTH","YRD","PRD"),
           ~ X_train%>%
             points_mask(x = "X_Lon", y = "Y_Lat", polygon = .x)%>%
             mutate(region = .y))
fwrite(X_train_SR,"Train_Raw/Resampled_Train_Raw_For_Specific_Region.csv")


##========Sites observations during 2007-2020========
X_train <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time))%>%
  filter(Label == 0)%>%
  dplyr::select(all_of(target_name))

Site_info <- distinct(X_train, GridID)

Grid <- fread("data/Grid_Information.csv")%>%
  rast(type = "xyz",crs = "WGS84")

ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"
OMI_path <- "D:/B620_ZW/OMI_Completed/2-ground_est/result/OMI_prediction_csv/"
target_name
for(year in 2007:2012){
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
      
      Site_info%>%
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
  
      Site_info%>%
        left_join(ERA5_daily)
    })
  future:::ClusterRegistry("stop")
  
  ## Landuse classification
  lucc <- fread(paste0("D:/B620_ZW/Landuse/full/lucc_full_",year,".csv"))%>%
    dplyr::select(c('GridID',"X_Lon","Y_Lat","crop_ratio",
                    "forest_ratio","grass_ratio",
                    "water_ratio","constr_ratio"))
  
  ## Modify time variable to indicate day distance between 2013-01-01
  ## Historical Site Grid data
  HS_data <- list(OMI, ERA5, lucc, Pop%>%filter(Year == year))%>%
    reduce(left_join)%>%
    mutate(Date = strptime(paste(Year, DOY), format="%Y %j")%>%as.Date(),
           Time = difftime(Date,as.Date("2013-01-01"), units = "days")%>%as.numeric())%>%
    dplyr::select(all_of(setdiff(target_name,c("Site_NO2","Label"))))
  fwrite(HS_data%>%na.omit, paste0("Train_Raw/Historical_Sites_Dataset_",year,".csv"))
}


##========Predictor distribution========
## remove the annotation to calculate the annual mean predictor
plan(multisession, workers = 1)
Predictor_annual_mean <- 2007:2020%>%
  future_map_dfr(function(year){
    cat(paste0("====",year,"====\n"))
    
    ## New OMI data
    OMI <- list.files(paste0(OMI_path,year),full.names = TRUE)%>%
      map_dfr(function(file){
        doy <- substr(file,96,98)%>%as.numeric()
        OMI_daily <- fread(file)%>%
          rename(ColumnNO2_Pre = ColumnNO2)%>%
          mutate(ColumnNO2 = ifelse(!is.na(ColumnNO2_True),
                                    ColumnNO2_True,ColumnNO2_Pre))%>%
          mutate(DOY = doy)%>%
          dplyr::select(GridID,ColumnNO2,DOY)
        OMI_daily
      })
    
    ## ERA5 data
    ERA5 <- list.files(paste0(ERA5_path,year),full.names = TRUE)%>%
      map_dfr(function(file){
        doy <- yday(substr(file,70,79))
        ERA5_daily <- fread(file)%>%mutate(DOY = doy)%>%
          dplyr::select(-X_Lon,-Y_Lat)
      })

    lucc <- fread(paste0("D:/B620_ZW/Landuse/full/lucc_full_",year,".csv"))%>%
      dplyr::select(c('GridID',"X_Lon","Y_Lat","crop_ratio",
                      "forest_ratio","grass_ratio",
                      "water_ratio","constr_ratio"))
    
    left_join(OMI, ERA5)%>%
      group_by(GridID)%>%
      summarise_all(mean, na.rm = TRUE)%>%
      ungroup()%>%mutate(Year = year)%>%
      left_join(lucc)%>%left_join(Pop%>%filter(Year == year))

  })
fwrite(Predictor_annual_mean, "result/Data_analysis/Predictor_annual_mean.csv")
