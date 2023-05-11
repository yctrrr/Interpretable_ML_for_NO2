setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)
library(furrr)
library(Cairo)
source("function/Sf_visualization.R")
source("function/Grid_visualization.R")
polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

## read in tuning model
lrn_xgboost_tune <- read_rds("result/Model/Final_model.Rds")
lrn_xgboost_tune$model$feature_names
fnames <- lrn_xgboost_tune$model$feature_names

## read in predictors
Grid <- fread("data/Grid_Information.csv")%>%
  dplyr::select(-X_Lon,-Y_Lat)
Pop <- fread("data/Urban_Pop_0.1g_2000-2020.csv")%>%
  rename(X_Lon = x, Y_Lat = y,
         Pop = Pop_all, UPop = Pop_U, RPop = Pop_O)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop,UPop,RPop)

#========Daily prediction========
lrn_xgboost_tune$model$feature_names
ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"
OMI_path <- "D:/B620_ZW/OMI_Completed/2-ground_est/result/OMI_prediction_csv/"
year = 2007
for (year in c(2007:2020)) {
  doy_list <- 1:yday(as.Date(paste0(year,"-12-31")))
  cat("Processing NO2 prediction in year: ",year,"\n")
  ## Landuse classification
  lucc <- fread(paste0("D:/B620_ZW/Landuse/full/lucc_full_",year,".csv"))%>%
    dplyr::select(c('GridID',"X_Lon","Y_Lat","crop_ratio",
                    "forest_ratio","grass_ratio",
                    "water_ratio","constr_ratio"))

  future::plan(multisession, workers = 10)
  doy_list%>%
    future_map(function(i){
      doy <- formatC(i,width = 3,flag = "0")
      date <- strptime(paste(year, i), format="%Y %j")%>%as.character()
      
      OMI <- fread(paste0(OMI_path,year,"/OMI_final_prediction_0.1g_",doy,".csv"))%>%
        rename(ColumnNO2_Pre = ColumnNO2)%>%
        mutate(ColumnNO2 = ifelse(!is.na(ColumnNO2_True),
                                  ColumnNO2_True,ColumnNO2_Pre))%>%
        dplyr::select(GridID,ColumnNO2)
      
      ERA5 <- fread(paste0(ERA5_path,year,"/era5-reanalysis-single-an-g0.25-",date,".csv"))
      
      NO2_predictorc <- list(OMI, ERA5, lucc, Pop%>%filter(Year == year))%>%
        reduce(left_join)%>%
        mutate(Date = as.Date(date),
               Time = difftime(Date,as.Date("2013-01-01"), units = "days")%>%as.numeric(),
               Month = month(Date), DOY = yday(Date), Year = year,
               Utime = as.POSIXct(Date)%>%
                 as.numeric())%>%
        na.omit()

      ## prediction
      prediction <- lrn_xgboost_tune$predict_newdata(NO2_predictorc)
      
      ## remove value less than 0
      ## save prediction
      NO2_predictorc%>%
        dplyr::select(X_Lon,Y_Lat)%>%ungroup()%>%
        mutate(NO2_response = prediction$response)%>%
        mutate(NO2_response = ifelse(NO2_response<0, 0, NO2_response))%>%
        write.csv(file = paste0("result/NO2_prediction_csv/",year,
                                "/NO2_prediction_0.1g_",year,"_",doy,".csv"),row.names = FALSE)
    })
  future:::ClusterRegistry("stop")
}

## Aggregate to annual prediction
NO2_annual <- 2007:2020%>%
  furrr::future_map_dfr(function(year){
    dir_path <- paste0("result/NO2_prediction_csv/",year)
    doy_list <- list.files(dir_path)%>%
      substr(26,28)

  doy_list%>%
    map_dfr(~fread(paste0(dir_path,"/NO2_prediction_0.1g_",year,"_",.x,".csv")))%>%
    dplyr::group_by(X_Lon,Y_Lat)%>%
    dplyr::summarise(NO2_response = mean(NO2_response), .groups = "drop")%>%ungroup()%>%
    as.data.frame()%>%mutate(Year = year)
  },.options = furrr_options(seed = TRUE))
write.csv(NO2_annual, paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_2007-2020.csv"), row.names = FALSE)

## Aggregate to monthly prediction
NO2_monthly <- 2007:2020%>%
  furrr::future_map_dfr(function(year){
    dir_path <- paste0("result/NO2_prediction_csv/",year)
    doy_list <- list.files(dir_path)%>%
      substr(26,28)
    
    doy_list%>%
      map_dfr(~fread(paste0(dir_path,"/NO2_prediction_0.1g_",year,"_",.x,".csv"))%>%
                mutate(month = strptime(paste(year,.x), format="%Y %j")%>%
                       lubridate::month()))%>%
      group_by(X_Lon,Y_Lat,month)%>%
      summarise(NO2_response = mean(NO2_response))%>%ungroup()%>%
      as.data.frame()%>%mutate(Year = year)
  })
write.csv(NO2_monthly, paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_monthly_2007-2020.csv"), row.names = FALSE)

## Aggregate to seasonal prediction
NO2_seasonal <- fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_monthly_2007-2020.csv"))%>%
  mutate(Season = ifelse(month%in%c(1,2,12),"Winter",
                  ifelse(month%in%c(3,4,5),"Spring",
                  ifelse(month%in%c(6,7,8),"Summer","Autumn"))))%>%
  mutate(Season = factor(Season, levels = c("Spring","Summer","Autumn","Winter")))%>%
  group_by(X_Lon,Y_Lat,Season)%>%
  summarise(NO2_response = mean(NO2_response))%>%ungroup()%>%
  as.data.frame()
write.csv(NO2_seasonal, paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_seasonal_2007-2020.csv"), row.names = FALSE)

#========Spatial distribution========
library(terra)
spplot_set <- function(base){
  base+
    scale_fill_manual(name = expression(paste("NO"[2]," concentration (¦Ìg/m"^"3",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.hjust = 0.5,
                                           byrow = TRUE, ticks = F))+
    geom_sf(data = polygon_sf, fill = NA, linewidth = 0.2)+
    geom_sf(data = Nline_sf, fill = NA, linewidth = 0.5)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    scale_x_continuous(breaks = c(90,110))+
    theme(axis.text = element_text(size = 25),
          strip.background = element_blank(),
          legend.key.width=unit(1.2, "cm"),
          legend.key.height=unit(1.2, "cm"),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 28),
          plot.margin=unit(c(5,0,5,0),"mm"),
          panel.grid.major = element_line(colour = alpha("black", 0.2), linewidth = .2))
}

breaks <- c(-Inf,seq(15,45,5),+Inf)
labels<- label_gen(breaks)
colors <- RColorBrewer::brewer.pal(length(labels),"Spectral")%>%rev()


## result visualization
coords <- 2013:2020%>%
  map_dfr(~fread(paste0("result/Data_analysis/Training_dataset/Training_NO2_0.1g_",.x,".csv")))%>%
  distinct(GridID, X_Lon, Y_Lat)

NO2_station <- 2013:2020%>%
  map_dfr(~ fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
      filter(Year == .x, Raw_label == 1)%>%
      group_by(GridID)%>%
      summarise(Site_NO2 = mean(Site_NO2))%>%ungroup()%>%
      left_join(coords)%>%dplyr::select(X_Lon,Y_Lat,Site_NO2)%>%
      dfproj(method = "bilinear")%>%mutate(Year = .x))

NO2_annual <- 2007:2020%>%
  map_dfr(function(year){
    fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_2007-2020.csv"))%>%
      filter(Year == year)%>%dplyr::select(-Year)%>%
      dfproj(method = "bilinear")%>%mutate(Year = year)
    })

NO2_annual%>%
  mutate(Year_Label = ifelse(Year < 2013, "2007-2012","2013-2020"))%>%
  group_by(Year)%>%
  summarise(NO2 = mean(NO2_response))

NO2_station_mod <- NO2_station%>%
  mutate(Site_NO2 = cut(Site_NO2, breaks, labels))

summary(NO2_annual)
NO2_annual_mod <- NO2_annual%>%
  left_join(NO2_station)%>%
  mutate_at(c("NO2_response","Site_NO2"), ~ cut(.,breaks,labels))

## annual NO2 during 2007-2020
NO2_annual_plot <- NO2_annual_mod%>%
  grid_facet_plot(x = "x", y = "y",
                  value = "NO2_response", facet = "Year",
                  facet_row = 4, facet_text = 30)%>%
  spplot_set()

CairoPNG(paste0("result/NO2_prediction/NO2_annual_plot.png"), 
         width = 3500,height = 2800,res = 160)
print(NO2_annual_plot)
dev.off()

## Grid estimate + Sites observation
for (year in c(2007:2020)) {
  NO2_year_mod <- NO2_annual_mod%>%
    filter(Year == year)
  ## Map purpose only. In case values are all less than 50
  NO2_year_mod[1,"NO2_response"] <- ">50"

  base_map <- NO2_year_mod%>%
    grid_plot(x = "x", y = "y", value = "NO2_response")+
    scale_fill_manual(name = expression(paste("NO"[2]," concentration (¦Ìg/m"^"3",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.vjust = 0.5,
                                           byrow = TRUE, ticks = F))+
    geom_sf(data = polygon_sf, fill = NA, linewidth = 0.3)+
    geom_sf(data = Nline_sf, fill = NA, linewidth = 0.5)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(axis.text = element_text(size = 30),
          legend.title = element_text(colour = "black", size = 35, angle = -90),
          legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),##legend height
          plot.margin=unit(c(5,0,5,0),"mm"),
          panel.border = element_rect(size = 0.5,colour = "black", fill = NA),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = alpha("black", 0.2), linewidth = .2))
  if(year < 2013){
    NO2_annual_plot <- base_map
  }else{
    NO2_annual_plot <- base_map+
      geom_point(data = NO2_station_mod%>%filter(Year == year),
                 aes(x, y, fill = Site_NO2),
                 colour = "black", shape = 21, size = 2.5, stroke = .2)
  }
  
  CairoPNG(paste0("result/NO2_prediction/NO2_annual_plot_",year,".png"), 
           width = 2500,height = 1500,res = 160)
  print(NO2_annual_plot)
  dev.off()
}

##========Multiyear average plot========
NO2_annual <- 2007:2020%>%
  map_dfr(function(year){
    fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_2007-2020.csv"))%>%
      filter(Year == year)%>%dplyr::select(-Year)%>%
      dfproj(method = "bilinear")%>%mutate(Year = year)
  })

breaks <- c(-Inf,seq(15,45,5),+Inf)
labels<- label_gen(breaks)
colors <- RColorBrewer::brewer.pal(length(labels),"Spectral")%>%rev()

## multiyear average
NO2_mya <- NO2_annual%>%
  group_by(x,y)%>%
  summarise(NO2_response = mean(NO2_response))%>%
  ungroup()

NO2_mya_mod <- NO2_mya%>%
  mutate_at(c("NO2_response"), ~ cut(.,breaks,labels))

NO2_mya_plot <- NO2_mya_mod%>%
  grid_plot(x = "x", y = "y", value = "NO2_response")%>%
  spplot_set()+
  scale_fill_manual(name = expression(paste("NO"[2]," concentration (¦Ìg/m"^"3",")")),
                    values = colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                         title.position="left",title.hjust = 0.5,
                                         byrow = TRUE, ticks = F))+
  scale_x_continuous(breaks = c(80,90,100,110,120))+
  theme(axis.text = element_text(size = 28),
        legend.title = element_text(colour = "black", size = 38, angle = 90),
        legend.text = element_text(size = 28),
        legend.key.width=unit(1.2, "cm"),
        legend.key.height=unit(1.2, "cm"),##legend height
        legend.position = "left",
        panel.border = element_rect(size = 1,colour = "black", fill = NA),
        panel.grid.major = element_line(colour = alpha("black", 0.2), linewidth = .2))

CairoPNG(paste0("result/NO2_prediction/NO2_multiyear_average_plot.png"), 
         width = 2500,height = 1500,res = 160)
print(NO2_mya_plot)
dev.off()


## Three hotspots
SRU1 <- st_read(dsn = "data/specific_region/jingjinji", layer = "jingjinji")%>%
  st_union()
SRU2 <- st_read(dsn = "data/specific_region/changsanjiao", layer = "changsanjiao")%>%
  st_union()
SRU3 <- st_read(dsn = "data/specific_region/zhusanjiao", layer = "zhusanjiao")%>%
  st_union()

## coordinate range for BTH,YRD and PRD
SR_range <- data.table(
  xlim = list(c(500000,1500000),c(1000000,1800000),c(500000,1200000)),
  ylim = list(c(3850000,4700000),c(3000000,3700000),c(2150000,2650000)),
  SR = c("BTH","YRD","PRD"),linewidth = c(1.2,1.2,1.2))

list(SRU1,SRU2,SRU3)%>%
  map2_dfr(.x = ., .y = 1:3,function(x,y){
    SR = SR_range$SR[y]
    NO2_mya_SR <- NO2_mya_plot+
      scale_x_continuous(limits = SR_range$xlim[[y]],expand = c(0, 0))+
      scale_y_continuous(limits = SR_range$ylim[[y]],expand = c(0, 0))+
      geom_sf(data = polygon_sf, fill = NA, linewidth = 0.5)+
      geom_sf(data = x, fill = NA, linewidth = SR_range$linewidth[y])+
      theme(legend.position = "none",
            panel.grid.major = element_blank(),
            panel.border = element_rect(size = 2.5,colour = "black", fill = NA),
            axis.text = element_blank())
    
    CairoPNG(paste0("result/NO2_prediction/NO2_multiyear_average_",SR,"_plot.png"), 
             width = 2800,height = 2800,res = 160)
    print(NO2_mya_SR)
    dev.off()
  })


##========Seasonal spatial plot========
season_label <- c("Spring","Summer","Autumn","Winter")
NO2_seasonal <- season_label%>%
  map_dfr(function(season){
    fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_seasonal_2007-2020.csv"))%>%
      filter(Season == season)%>%dplyr::select(-Season)%>%
      dfproj(method = "bilinear")%>%mutate(Season = season)
  })%>%mutate(Season = factor(Season, levels = season_label))

NO2_seasonal%>%
  group_by(Season)%>%
  summarise(mean = mean(NO2_response), sd = sd(NO2_response))

breaks <- c(-Inf,seq(15,50,5),+Inf)
labels<- label_gen(breaks)
colors <- RColorBrewer::brewer.pal(length(labels),"Spectral")%>%rev()

NO2_seasonal_mod <- NO2_seasonal%>%
  mutate_at(c("NO2_response"), ~ cut(.,breaks,labels))


NO2_seasonal_plot <- NO2_seasonal_mod%>%
  grid_facet_plot(x = "x", y = "y",
                  value = "NO2_response", facet = "Season", 
                  facet_row = 2, facet_text = 25)%>%
  spplot_set()+
  theme(axis.text = element_text(size = 20),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 25),
        legend.key.width=unit(1, "cm"),
        legend.key.height=unit(1, "cm"))

CairoPNG(paste0("result/NO2_prediction/NO2_seasonal_plot.png"), 
         width = 2400,height = 1500,res = 160)
print(NO2_seasonal_plot)
dev.off()

##=======Specific region calculation======
source("function/Spatial_analysis.R")
NO2_annual <- fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_2007-2020.csv"))

SR1 <- st_read(dsn = "data/specific_region/jingjinji", layer = "jingjinji")
SR2 <- st_read(dsn = "data/specific_region/changsanjiao", layer = "changsanjiao")
SR3 <- st_read(dsn = "data/specific_region/zhusanjiao", layer = "zhusanjiao")

NO2_SR_annual <- list(SR1,SR2,SR3)%>%
  map2_dfr(.x = ., .y = c("BTH","YRD","PRD"),
           ~ NO2_annual%>%
             dplyr::select(X_Lon,Y_Lat,NO2_response,Year)%>%
             points_mask(x = "X_Lon", y = "Y_Lat", polygon = .x)%>%
             mutate(region = .y))
NO2_annual_join <- NO2_annual%>%
  left_join(NO2_SR_annual)%>%
  mutate(region = ifelse(is.na(region),"Other",region))

NO2_annual_join%>%
  count(region)
NO2_SR_annual%>%
  count(region)

NO2_annual_join%>%
  group_by(region)%>%
  summarise(NO2 = mean(NO2_response),sd = sd(NO2_response))%>%ungroup()
  summarise(NO2 = mean(NO2_response))%>%ungroup()%>%
  pivot_wider(names_from = region, values_from = NO2)
