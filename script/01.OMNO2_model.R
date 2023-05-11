## Data analysis for OMI NO2
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)
library(geosphere)
library(terra)
library(furrr)

Grid <- fread("data/Grid_Information.csv")
##========Match training dataset========
year = 2007
GOME_path <- "D:/B620_ZW/OMI_Completed/GOME-2_and_Merra/"
OMI_path <- "D:/B620_ZW/OMINO2_to_point_one_grid/"
Overlap_path <- "D:/B620_ZW/OMI_and_GOME-2/"
ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"

for(year in 2007:2020){
  plan("multisession",workers = 6)
  ## read in GOME and OMI data
  GOME <- list.files(paste0(GOME_path,year),full.names = TRUE)%>%
    map_dfr(~fread(.x))%>%dplyr::select(GridID,DOY,NO2_GOME)

  OMI <- list.files(paste0(OMI_path,year),full.names = TRUE)%>%
    future_map_dfr(~fread(.x)%>%na.omit())%>%
    dplyr::select(GridID,DOY,ColumnNO2)%>%
    distinct(GridID,DOY,.keep_all = TRUE)

  Overlap <- list.files(paste0(Overlap_path,year),full.names = TRUE)%>%
    map_dfr(~fread(.x))%>%
    dplyr::select(Grid_ID,DOY,NO2_GOME,NO2_OMI)%>%
    rename(GridID = Grid_ID)

  TROPNO2 <- left_join(OMI,Overlap)%>%
    dplyr::select(-NO2_OMI)%>%
    full_join(GOME)

  ## read in ERA5 data and join with observations
  Train_data <- list.files(paste0(ERA5_path,year),full.names = TRUE)%>%
    future_map_dfr(function(file){
      doy <- yday(substr(file,70,79))

      ERA5_daily <- fread(file)%>%mutate(DOY = doy)
      TROPNO2%>%
        filter(DOY == doy)%>%
        left_join(ERA5_daily)
    })

  fwrite(Train_data, paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))
}

##========Step 1 Data fusion=========
## OMI and GOME2 matching
lrn_rf <- lrn("regr.ranger",num.trees = 500,
              min.node.size = 1, num.threads = 10, 
              importance = 'impurity')

Fusion_18_20 <- 2018:2020%>%
  map_dfr(~fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",.x,".csv"))%>%
            na.omit())%>%dplyr::select(-X_Lon,-Y_Lat,-GridID)

Fusion_task_18_20 <- TaskRegr$new("Fusion_18_20", Fusion_18_20, target = "ColumnNO2")

lrn_rf_fusion_18_20 <- lrn_rf$train(Fusion_task_18_20)

## Model performance
Fusion_mp <- data.frame(OOB_Rsquare = numeric(),MSE = numeric(),Year = numeric())
## Data availability (available days)
Fusion_ad <- data.frame(DOY = numeric(),Year = numeric())

for(year in 2007:2020){
  if(year >= 2018){
    lrn_rf_fusion <- lrn_rf_fusion_18_20
  }else{
    Fusion_train <- fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))%>%
      na.omit()%>%dplyr::select(-X_Lon,-Y_Lat,-GridID)
    
    Fusion_task <- TaskRegr$new("OMI_singleyear", Fusion_train, target = "ColumnNO2")

    lrn_rf_fusion <- lrn_rf$train(Fusion_task)
  }

  mp <- data.frame(OOB_Rsquare = lrn_rf_fusion$model$r.squared,
                   MSE = lrn_rf_fusion$oob_error(),Year = year)
  ## Model performance
  Fusion_mp <- full_join(Fusion_mp, mp)

  ## Prediction
  Fusion_fill <- fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))%>%
    filter(is.na(ColumnNO2))%>%dplyr::select(-ColumnNO2)%>%na.omit()
  New_OMI_pred <- lrn_rf_fusion$predict_newdata(Fusion_fill)

  Fusion_result <- Fusion_fill%>%
    dplyr::select(GridID, DOY)%>%
    mutate(ColumnNO2 = New_OMI_pred$response)
  
  ## Available days
  ad <- count(Fusion_result, DOY)%>%
    mutate(Year = year, Date = as.Date(strptime(paste(year, DOY), format="%Y %j")))

  Fusion_ad <- full_join(Fusion_ad, ad)
  
  fwrite(Fusion_result, paste0("result/Data_analysis/OMI_training_dataset/New_OMI_0.1g_",year,".csv"))
}
## Model performance
fwrite(Fusion_mp, paste0("result/OMI_prediction/New_OMI_model_performance.csv"))
## Available days
fwrite(Fusion_ad, paste0("result/OMI_prediction/New_OMI_available_days.csv"))

## Merge step1 result
## Raw OMI data
Fusion_ad <- fread(paste0("result/OMI_prediction/New_OMI_available_days.csv"))

Raw_ad <- 2007:2020%>%
  map_dfr(function(year){
    fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))%>%
      dplyr::select(GridID, DOY, ColumnNO2)%>%na.omit()%>%
      count(DOY, name = "n_OMI")%>%mutate(Year = year)
  })

Fusion_ad%>%
  filter(Year == 2007, DOY==365)

Fusion_ad_comp <- full_join(Fusion_ad,Raw_ad)%>%
  mutate(Date = as.Date(strptime(paste(Year, DOY), format="%Y %j")))%>%
  mutate_at(c("n","n_OMI"),~replace_na(.x, 0))%>%
  fwrite(paste0("result/OMI_prediction/New_OMI_available_days.csv"))


##========Summary of model fusion in Step 1=========
Model_OOB <- fread(paste0("result/OMI_prediction/New_OMI_model_performance.csv"))

Data_enhance <- fread(paste0("result/OMI_prediction/New_OMI_available_days.csv"))%>%
  group_by(Year)%>%summarise(n = sum(n),n_OMI = sum(n_OMI))%>%
  ungroup()%>%mutate(Inc = n/n_OMI*100)%>%
  left_join(Model_OOB)%>%
  mutate(RMSE = MSE^(1/2))

fwrite(Data_enhance,paste0("result/OMI_prediction/New_OMI_performance.csv"))
Data_enhance
##========Step 2 Gap filling=========
## Require a lot of time and computational memory
ERA5_path <- "D:/B620_SYC/data/ERA5/daily_csv/"
year = 2007
month = 3
for(year in 2007:2020){
  ## Model performance
  Final_mp <- data.frame(OOB_Rsquare = numeric(),MSE = numeric(),
                         Year = numeric(),Month = numeric())
  ## Training dataset
  Raw_OMI <- fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))%>%
    dplyr::select(GridID, DOY, ColumnNO2)%>%na.omit()%>%mutate(Fusion = 0)
  
  Cov_OMI <- fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_",year,".csv"))%>%
    dplyr::select(-NO2_GOME, -ColumnNO2)%>%na.omit()
  
  New_OMI <- fread(paste0("result/Data_analysis/OMI_training_dataset/New_OMI_0.1g_",year,".csv"))%>%
    mutate(Fusion = 1)
  ## Training model
  Final_train <- full_join(Raw_OMI,New_OMI)%>%
    left_join(Cov_OMI)%>%dplyr::select(-X_Lon,-Y_Lat)%>%
    mutate(Month = month(as.Date(strptime(paste(year, DOY), format="%Y %j"))))
  
  ## Fusion result(stage 1)
  Fuse_OMI <- full_join(Raw_OMI,New_OMI)%>%
    rename(ColumnNO2_True = ColumnNO2)
  
  ## Prediction dir
  if(!dir.exists(paste0("result/OMI_prediction_csv/",year))){
    dir.create(paste0("result/OMI_prediction_csv/",year))
  }
  
  ## Prediction dataset (ERA 5)
  month_list <- list.files(paste0(ERA5_path,year),full.names = TRUE)%>%
    substr(70,79)%>%month()
  for (month in 1:12) {
    cat("====Processing ",year," Month ",month,"====\n")
    Sub_train <- Final_train%>%filter(Month == month)
    Final_task <- TaskRegr$new("OMI_singleyear", Sub_train, 
                               target = "ColumnNO2")
    
    Final_task$col_roles$feature <- setdiff(Final_task$feature_names, 
                                            c("Fusion","GridID","Month"))

    lrn_rf <- lrn("regr.ranger",num.trees = 500,
                  min.node.size = 1, num.threads = 8)
    lrn_rf_final <- lrn_rf$train(Final_task)

    mp <- data.frame(OOB_Rsquare = lrn_rf_final$model$r.squared,
                     MSE = lrn_rf_final$oob_error(),
                     Year = year, Month = month)
    
    ## Model performance
    Final_mp <- full_join(Final_mp, mp)
    
    ## Monthly prediction
    ERA5_list <- list.files(paste0(ERA5_path,year),full.names = TRUE)[month_list==month]
    ERA5_list%>%map(function(file){
      
        doy <- yday(substr(file,70,79))
        
        ERA5_daily <- fread(file)%>%mutate(DOY = doy)%>%
          mutate(Month = month(as.Date(strptime(paste(year, DOY), format="%Y %j"))))
        
        ## May take a long time
        Final_pred <- lrn_rf_final$predict_newdata(ERA5_daily)
        
        Gap_OMI <- ERA5_daily%>%
          dplyr::select(GridID)%>%
          mutate(ColumnNO2 = Final_pred$response)
        
        ## combine raw data
        Fuse_OMI%>%filter(DOY == doy)%>%
          dplyr::select(-DOY)%>%
          left_join(Gap_OMI,.)%>%
          fwrite(paste0("result/OMI_prediction_csv/",
                        year,"/OMI_final_prediction_0.1g_",
                        formatC(doy,width = 3,flag = "0"),".csv"))
      })
  }
  rm(lrn_rf_final)
  ## Model performance
  fwrite(Final_mp, paste0("result/OMI_prediction_csv/OMI_final_model_performance_",year,".csv"))
}

## Monthly mean
OMI_path <- "D:/B620_ZW/OMI_Completed/2-ground_est/result/OMI_prediction_csv/"
OMI_monthly <- 2007:2020%>%
  map_dfr(function(year){
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
    
    OMI%>%
      mutate(Month = month(strptime(paste(year,DOY), format="%Y %j")))%>%
      group_by(GridID, Month)%>%
      summarise(ColumnNO2 = mean(ColumnNO2))%>%ungroup()%>%
      mutate(Year = year)
  })
write.csv(OMI_monthly, paste0("result/OMI_prediction_csv/OMNO2_prediction_0.1g_monthly_2007-2020.csv"), row.names = FALSE)


## Typical days
Grid_info <- fread("data/Grid_Information.csv")%>%
  dplyr::select(X_Lon,Y_Lat,GridID)

GOME_td <- fread(paste0("result/Data_analysis/OMI_training_dataset/Training_OMI_0.1g_2015.csv"))%>%
  dplyr::select(GridID, DOY, NO2_GOME)%>%na.omit()%>%filter(DOY%in%c(1:3))

1:3%>%
  map_dfr(~fread(paste0("result/OMI_prediction_csv/2015/OMI_final_prediction_0.1g_",
                        formatC(.x,width = 3,flag = "0"),".csv"))%>%
            mutate(DOY = .x))%>%
  left_join(GOME_td)%>%left_join(Grid_info)%>%
  fwrite(paste0("result/OMI_prediction_csv/OMI_final_prediction_20150101-03_0.1g.csv"))



##========Model performance summary========
Final_mp <- 2007:2020%>%
  map_dfr(~fread(paste0("result/OMI_prediction_csv/OMI_final_model_performance_",.x,".csv"))%>%
            mutate(RMSE = MSE^(1/2)))

Final_mp_plot <- ggplot(Final_mp)+
  geom_boxplot(aes(x = Month,y = OOB_Rsquare,group = Month,color = "R2"),
               lwd = .8)+
  geom_boxplot(aes(x = Month,y = RMSE/scale,group = Month,color = "RMSE"),
               lwd = .8)+
  scale_x_continuous(name = "Month", breaks = 1:12)+
  scale_y_continuous(name = expression(paste("OOB  R "^"2")),
                     sec.axis = sec_axis(~.*scale, 
                                         name = expression(paste("RMSE (10"^"15"," mol/cm"^"2",")"))))+
  scale_color_manual(values = c("#feb24c", "#31a354"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth = 0.5, fill = NA, colour = 'black'),
        panel.border = element_rect(linewidth = 1.5, colour = "black", fill = NA)
  )+
  theme(axis.text = element_text(size = 25),
        axis.text.y = element_text(color = "#feb24c"),
        axis.text.y.right = element_text(color = "#31a354"),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28, color = "#feb24c"),
        axis.title.y.right = element_text(color = "#31a354"),
        legend.text = element_text(size = 25),
        legend.position = "none")
ggsave("result/OMI_prediction/OMI_Final_Model_Performance.png",Final_mp_plot,
       width = 12, height = 7)

##========OMI on typical days plot=========
library(Cairo)
source("function/Sf_visualization.R")
source("function/Grid_visualization.R")
spplot_set <- function(base){
  base+
    scale_fill_manual(name = expression(paste("NO"[2]," VCD (10"^"15"," mol/cm"^"2",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.vjust = 0.5,
                                           byrow = TRUE, ticks = F,
                                           override.aes = list(size=10)))+
    geom_sf(data = polygon_sf, fill = NA, linewidth = 0.2)+
    geom_sf(data = Nline_sf, fill = NA, linewidth = 0.5)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(strip.text.x = element_text(
      size = 35, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
      strip.background = element_rect(color=NA, fill = NA))+
    theme(axis.text = element_text(size = 25),
          legend.title = element_text(colour = "black", size = 30, angle = -90),
          legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),##legend height
          legend.text = element_text(colour = "black", size = 25),
          plot.margin=unit(c(5,0,5,0),"mm"),#
          panel.border = element_rect(size = 1.5,
                                      colour = "black", fill = NA),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = alpha("black", 0.2), linewidth = .2))
}

polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

## GOME and OMI fusion
OMI_final <- fread(paste0("result/OMI_prediction_csv/OMI_final_prediction_20150101-03_0.1g.csv"))%>%
  filter(DOY == 1)

Raw_OMI <- OMI_final%>%filter(Fusion == 0)%>%
  dplyr::select(X_Lon,Y_Lat,ColumnNO2_True)%>%na.omit()%>%
  rename(ColumnNO2 = ColumnNO2_True)%>%mutate(Data = "OMNO2")

GOME <- dplyr::select(OMI_final,X_Lon,Y_Lat,NO2_GOME)%>%na.omit()%>%
  rename(ColumnNO2 = NO2_GOME)%>%mutate(Data = "GOME-2")

Fusion_OMI <- dplyr::select(OMI_final,X_Lon,Y_Lat,ColumnNO2_True)%>%na.omit()%>%
  rename(ColumnNO2 = ColumnNO2_True)%>%mutate(Data = "Data Fusion")

Stage1 <- list(Raw_OMI,GOME,Fusion_OMI)%>%
  reduce(full_join)%>%
  mutate(Data = factor(Data, levels = c("OMNO2","GOME-2","Data Fusion")))%>%
  st_as_sf(coords = c("X_Lon","Y_Lat"),crs = "+proj=longlat +datum=WGS84")%>%
  st_transform(crs = "+proj=aea +lat_0=0 +lon_0=105 +lat_1=25 +lat_2=47 +x_0=0 +y_0=0 +datum=WGS84")

breaks <- c(-Inf,c(3,4,5,10,15,20,25,30),+Inf)
labels<- label_gen(breaks)
colors <- RColorBrewer::brewer.pal(length(labels),"OrRd")

Stage1_mod <- Stage1%>%
  mutate(value = cut(ColumnNO2, breaks, labels),.before = "geometry")

Stage1_plot <- ggplot(data = Stage1_mod)+
  ggplot2::geom_sf(mapping = aes(fill = value), shape = 21,
                   size = 2,stroke = 0)+
  facet_wrap("Data",nrow = 1)

CairoPNG(paste0("result/OMI_prediction/OMI_Fusion_201501_plot.png"),
         width = 3500,height = 1500,res = 160)
print(spplot_set(Stage1_plot))
dev.off()

##=======Spatial plot for OMI NO2========
## Please run Data_analysis to calculate covariaite annual mean value first
library(Cairo)
library(terra)
source("function/Sf_visualization.R")
source("function/Grid_visualization.R")
polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

spplot_set <- function(base){
  base+
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

Predictor_am <- fread("result/Data_analysis/Predictor_annual_mean.csv")

Pop <- fread("data/Urban_Pop_0.1g_2000-2020.csv")%>%
  rename(X_Lon = x, Y_Lat = y, 
         Pop = Pop_all, UPop = Pop_U, RPop = Pop_O)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop,UPop,RPop)

Predictor_am%>%
  left_join(Pop)%>%
  na.omit%>%
  group_by(Year)%>%
  summarise(OMI_Pop = weighted.mean(ColumnNO2,Pop),
            OMI = mean(ColumnNO2), OMI_sd = sd(ColumnNO2))

OMINO2 <- 2007:2020%>%
  map_dfr(~ Predictor_am%>%
            dplyr::select(X_Lon, Y_Lat, Year, ColumnNO2)%>%
            filter(Year == .x)%>%dplyr::select(-Year)%>%
            dfproj(method = "bilinear")%>%mutate(Year = .x))


summary(OMINO2)
breaks <- c(-Inf,c(3,4,5,10,15,20,25,30),+Inf)
labels<- label_gen(breaks)

colors <- RColorBrewer::brewer.pal(length(labels),"OrRd")

OMINO2_annual_mod <- OMINO2%>%
  mutate(value = cut(ColumnNO2, breaks, labels))
OMINO2_annual_mod


OMI_annual_plot <- OMINO2_annual_mod%>%
  grid_facet_plot(x = "x", y = "y",
                  value = "value", facet = "Year", facet_row = 4, facet_text = 30)+
  scale_fill_manual(name = expression(paste("NO"[2]," VCD (10"^"15"," mol/cm"^"2",")")),
                    values = colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                         title.position="right",title.hjust = 0.5,
                                         byrow = TRUE, ticks = F))

CairoPNG(paste0("result/OMI_prediction/OMI_NO2_annual_plot.png"), 
         width = 3500,height = 2800,res = 160)
print(spplot_set(OMI_annual_plot))
dev.off()

##=======Specific region calculation======
Predictor_am <- fread("result/Data_analysis/Predictor_annual_mean.csv")

SR1 <- st_read(dsn = "data/specific_region/jingjinji", layer = "jingjinji")
SR2 <- st_read(dsn = "data/specific_region/changsanjiao", layer = "changsanjiao")
SR3 <- st_read(dsn = "data/specific_region/zhusanjiao", layer = "zhusanjiao")

OMNO2_SR_annual <- list(SR1,SR2,SR3)%>%
  map2_dfr(.x = ., .y = c("BTH","YRD","PRD"),
           ~ Predictor_am%>%
             dplyr::select(X_Lon,Y_Lat,ColumnNO2,Year)%>%
             points_mask(x = "X_Lon", y = "Y_Lat", polygon = .x)%>%
             mutate(region = .y))

OMNO2_SR_annual%>%
  group_by(region,Year)%>%
  summarise(NO2 = mean(ColumnNO2))%>%
  pivot_wider(names_from = region, values_from = NO2)
