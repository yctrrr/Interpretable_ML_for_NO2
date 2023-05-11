setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)
library(geosphere)

##======Read in dataset======
target_name <- c('GridID', 'DOY', "Month",'Year','Raw_label','Site_NO2','ColumnNO2',
                 'u10', 'v10', 't2m', 'blh', 'tcc', 'tco3',
                 'sp', 'e', 'tcw', 'lsp', 'lai_hv', 'lai_lv',
                 'crop_ratio', 'forest_ratio','grass_ratio', 'water_ratio', 'constr_ratio','Pop')
Grid <- fread("data/Grid_Information.csv")%>%
  dplyr::select(-X_Lon,-Y_Lat)

NO2_resample <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
            mutate(Month = month(as.Date("2013-01-01")+Time),
                   DOY = yday(as.Date("2013-01-01")+Time),
                   Utime = as.POSIXct(as.Date("2013-01-01")+Time)%>%
                     as.numeric())%>%
            dplyr::select(all_of(target_name))
NO2_resample$Site_NO2%>%mean
NO2_raw <- NO2_resample%>%filter(Raw_label == 1)

NO2_task <- TaskRegr$new("NO2_multiyear", NO2_resample, target = "Site_NO2")
NO2_raw_task <- TaskRegr$new("NO2_multiyear", NO2_raw, target = "Site_NO2")

## After arbitrary selection
rmv <- c("TQI","TROPPB","QV2M","TO3","U2M","V2M","unused_ratio")
feature_clean <- setdiff(NO2_task$feature_names,rmv)
NO2_task$select(feature_clean)
NO2_task$filter(which(complete.cases(NO2_task$data())))
NO2_raw_task$select(feature_clean)

##  resample methods
cv10 <- rsmp("cv", folds = 10)
cv3 <- rsmp("cv", folds = 3)
by_year <- rsmp("loo")
by_site <- rsmp("cv", folds = 10)
measure <- msrs(c("regr.rsq","regr.rmse"))

##  assign different tasks
NO2_overall_task <- NO2_task$clone()
NO2_overall_task$col_roles$feature <- setdiff(NO2_overall_task$feature_names, 
                                           c("GridID","Raw_label"))
NO2_overall_task

NO2_back_task <- NO2_task$clone()
NO2_back_task$col_roles$group <- "Year"
NO2_back_task$col_roles$feature <- setdiff(NO2_back_task$feature_names, 
                                           c("GridID","Raw_label"))
NO2_back_task

NO2_site_task <- NO2_task$clone()
NO2_site_task$col_roles$group = "GridID"
NO2_site_task$col_roles$feature <- setdiff(NO2_site_task$feature_names, 
                                           c("GridID","Raw_label"))

##  Task for raw dataset
NO2_overallraw_task <- NO2_raw_task$clone()
NO2_overallraw_task$col_roles$feature <- setdiff(NO2_overallraw_task$feature_names,
                                                 c("GridID","Raw_label"))
NO2_overallraw_task
NO2_backraw_task <- NO2_raw_task$clone()
NO2_backraw_task$col_roles$group <- "Year"
NO2_backraw_task$col_roles$feature <- setdiff(NO2_backraw_task$feature_names, 
                                              c("GridID","Raw_label"))

##======model tuning for xgboost======
lrn_xgboost_tune <- lrn("regr.xgboost", objective="reg:squarederror", 
                        eval_metric = "rmse", booster = "gbtree",
                        nrounds=500L, eta=0.1, gamma=0.1, nthread = 15)

search_space <- ps(
  max_depth = p_int(lower = 4, upper = 10),
  min_child_weight = p_dbl(lower = 3, upper = 8),
  subsample = p_dbl(lower = .5, upper = .9),
  colsample_bytree = p_dbl( lower = .5, upper = .9)
)

mbo_tuner <- tnr("mbo")
grid_tuner <- tnr("grid_search")

at_xgboost = mlr3verse::auto_tuner(
  tuner = tnr("mbo"),
  learner = lrn_xgboost_tune,
  resampling = by_year,
  measure = msr("regr.rsq"),
  search_space = search_space,
  term_evals = 10,
  store_models = TRUE
)

tic("Model tuning...")
NO2_tune_xgboost <- at_xgboost$train(task = NO2_back_task)
toc()

## This is only the model fitting result
NO2_tune_xgboost$predict(NO2_back_task)$score(measure)
tuning_result <- as.data.table(NO2_tune_xgboost$archive)%>%
  dplyr::select(-resample_result)
fwrite(tuning_result,"result/Model/Tuning_lookup_table.csv")
saveRDS(NO2_tune_xgboost$learner,"result/Model/Tuning_model.Rds")

## apply tunning result
# lrn_xgboost_tune <- read_rds("result/Model/Tuning_model.Rds")
# lrn_xgboost_done <- lrn("regr.xgboost")
# lrn_xgboost_done$param_set$values <- lrn_xgboost_tune$param_set$values

lrn_xgboost_done <- lrn("regr.xgboost", objective="reg:squarederror",
                        eval_metric = "rmse", booster = "gbtree",
                        nrounds=500L, eta=0.1, gamma=0.1, nthread = 15,
                        max_depth = 9, min_child_weight = 5.485718,
                        subsample = 0.5173684, colsample_bytree = 0.8428992)

tic("Sample-based CV")
xgboost_sample <- mlr3::resample(NO2_overall_task, lrn_xgboost_done, cv10)
xgboost_sample$prediction()$score(measure)
toc()
saveRDS(xgboost_sample, "result/Model/Xgboost_sample_CV.Rds")


tic("Site-based CV")
xgboost_site <- mlr3::resample(NO2_site_task, lrn_xgboost_done, by_site)
xgboost_site$prediction()$score(measure)
toc()
saveRDS(xgboost_site, "result/Model/Xgboost_site_CV.Rds")

tic("By-year CV")
xgboost_back <- mlr3::resample(NO2_back_task, lrn_xgboost_done, by_year)
xgboost_back$prediction()$score(measure)
toc()
saveRDS(xgboost_back, "result/Model/Xgboost_by_year_CV.Rds")

tic("Final model")
xgboost_final <- lrn_xgboost_done
learner <- xgboost_final$train(NO2_overall_task)
toc()
saveRDS(learner,"result/Model/Final_model.Rds")

tic("Sample-based CV for raw data")
xgboost_sampleraw <- mlr3::resample(NO2_overallraw_task, lrn_xgboost_done, cv10)
xgboost_sampleraw$prediction()$score(measure)
toc()
saveRDS(xgboost_sampleraw, "result/Model/Xgboost_raw_sample_CV.Rds")

# tic("By-year CV for raw data")
# xgboost_backraw <- mlr3::resample(NO2_backraw_task, lrn_xgboost_done, by_year)
# xgboost_backraw$prediction()$score(measure)
# toc()
# saveRDS(xgboost_backraw, "result/Model/Xgboost_raw_by_year_CV.Rds")

##========CV validation========
cv_measure <- function(df){
  data.table(R2 = cor(df$truth,df$response)^2,
             RMSE = sqrt(mean((df$truth - df$response)^2)))
}
cv_list <- c("sample_CV","site_CV","by_year_CV")

cv_combo <- 1:length(cv_list)%>%
  map_dfr(function(i){
    xgboost_cv <- read_rds(paste0("result/Model/Xgboost_",cv_list[i],".Rds"))
    xgboost_cv_df <- xgboost_cv$prediction()%>%as.data.table()
    
    xgboost_cv_dfd <- NO2_resample%>%
      as.data.table()%>%
      mutate(row_ids = 1:n(),
             Month = month(as.Date("2013-01-01")+Time),
             DOY = yday(as.Date("2013-01-01")+Time))%>%
      dplyr::select(row_ids,GridID,Time,DOY,Month,Year,Raw_label)%>%
      left_join(xgboost_cv_df)

    
    cor(xgboost_cv_dfd%>%filter(Raw_label == 1))^2
    xgboost_cv_dfm <- xgboost_cv_dfd%>%
      group_by(GridID,Month,Year)%>%
      summarise(truth = mean(truth), response = mean(response))%>%
      ungroup()
    
    xgboost_cv_dfa <- xgboost_cv_dfd%>%
      group_by(GridID,Year)%>%
      summarise(truth = mean(truth), response = mean(response))%>%
      ungroup()
    
    list(cv_measure(xgboost_cv_dfd)%>%mutate(level = "daily"),
         cv_measure(xgboost_cv_dfm)%>%mutate(level = "monthly"),
         cv_measure(xgboost_cv_dfa)%>%mutate(level = "annual"))%>%
      reduce(full_join)%>%
      mutate(method = cv_list[i])
  })

##========CV visualization========
library(Cairo)
source("function/Scatter_plot.R")
## Derive daily, monthly and annual value from mlr3 prediction
mlr_pred <- xgboost_cv_df

cv_derive <- function(data, mlr_pred, raw = TRUE){
  cv_df <- mlr_pred%>%as.data.table()%>%
    mutate(response = ifelse(response<0,0,response),
           truth = ifelse(truth<0,0,truth))

  if(raw == TRUE){
    cv_dfd <- data%>%
      as.data.table()%>%
      mutate(row_ids = 1:n())%>%
      dplyr::select(row_ids,GridID,DOY,Month,Year,Raw_label)%>%
      left_join(cv_df)%>%filter(Raw_label == 1)
  }else{
    cv_dfd <- data%>%
      as.data.table()%>%
      mutate(row_ids = 1:n())%>%
      dplyr::select(row_ids,GridID,DOY,Month,Year)%>%
      left_join(cv_df)
  }
  
  cv_dfm <- cv_dfd%>%
    group_by(GridID,Month,Year)%>%
    summarise(truth = mean(truth), response = mean(response))%>%
    ungroup()

  cv_dfa <- cv_dfd%>%
    group_by(GridID,Year)%>%
    summarise(truth = mean(truth), response = mean(response))%>%
    ungroup()
  
  list(cv_dfd%>%mutate(level = "Daily"),
       cv_dfm%>%mutate(level = "Monthly"),
       cv_dfa%>%mutate(level = "Annual"))%>%
    reduce(full_join)
}

## Setting for NO2 scatter plot
plot_set <- function(base){
  base+
    stat_smooth(method = "lm",
                se = FALSE,colour="black",
                size = .3, fullrange = TRUE)+
    geom_abline(slope = 1, linetype = "dashed",
                intercept=0, colour="blue", size = .3)+
    stat_binhex(bins = 300)+
    scale_fill_gradientn(name = expression("Count"),
                         trans = "log10",
                         colours = c("#114bf9","#0ffff4","#60f818",
                                     "#f2f607","#f2d01e","#f20505"),
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth=.8,
                                                title.position="top",title.hjust=0.5,title.vjust = 1,ticks= F))+
    scale_x_continuous(
      name = expression(paste("Observed NO"[2]," (¦Ìg/m"^"3",")")),
      expand = c(0.01, 0))+
    scale_y_continuous(
      name = expression(paste("Predicted NO"[2]," (¦Ìg/m"^"3",")")),
      expand = c(0.01, 0))+
    theme(strip.text.x = element_text(
      size = 28, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
      strip.text.y = element_text(size = 28, color = "black"),
      strip.background = element_rect(color="black", fill=NA, linewidth = .8))+
    theme(legend.position = "right",
          legend.key.width = unit(.6, "cm"),
          legend.key.height = unit(4, "cm"),
          legend.text = element_text(colour = "black", size = rel(2)),
          legend.title = element_text(colour = "black", size = rel(2.5)),
          legend.margin = margin(r = 10))+
    theme(panel.spacing.x = unit(2, "lines"))+
    theme(axis.text = element_text(size = 25),
          axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"), size = 30),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"), size = 30))
}

## Main program
cv_list <- c("sample_CV","site_CV","by_year_CV")
cv_list_mod <- c("Sample-based CV","Grid-based CV","By-year CV")
NO2_resample <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))

cv_point_combo <- 1:length(cv_list)%>%
  map_dfr(function(i){
    xgboost_cv <- read_rds(paste0("result/Model/Xgboost_",cv_list[i],".Rds"))

    cv_derive(NO2_resample, xgboost_cv$prediction())%>%
      mutate(method = cv_list_mod[i])
  })%>%
  rename(y = truth, ypred = response)%>%
  mutate(method = factor(method, levels = cv_list_mod),
         level = factor(level, levels = c("Daily","Monthly","Annual")))%>%
  as.data.frame()

cv_grid <- cv.plot.facet.grid(cv_point_combo, nrow = 2, cex = 7,
              alpha = 0.3, row = "method",col = "level",label = "", vjust = 1.1)%>%
  plot_set()+
  theme(legend.key.width = unit(.6, "cm"),
        legend.key.height = unit(2, "cm"))

CairoPNG(paste0("result/Model/CV_validation.png"),
         width = 5000,height = 4500,res=240)
cv_grid
dev.off()

##  Validation of model using raw data
xgboost_raw_cv <- read_rds(paste0("result/Model/Xgboost_raw_sample_CV.Rds"))
NO2_raw <- NO2_resample%>%filter(Raw_label == 1)
raw_cv_point <- cv_derive(NO2_raw,xgboost_raw_cv$prediction())%>%
  rename(y = truth, ypred = response)%>%
  mutate(level = factor(level, levels = c("Daily","Monthly","Annual")))%>%
  as.data.frame()

raw_cv_grid <- cv.plot.facet(raw_cv_point, nrow = 1, cex = 7,
                             alpha = 0.3, facet = "level", label = "", vjust = 1.1)%>%
  plot_set()+
  theme(legend.key.width = unit(.6, "cm"),
        legend.key.height = unit(2, "cm"))

CairoPNG(paste0("result/Model/CV_validation_raw_data.png"),
         width = 5000,height = 1800,res=240)
raw_cv_grid
dev.off()

##=======External validation of Taiwan=======
lrn_xgboost_tune <- read_rds("result/Model/Final_model.Rds")
NO2_taiwan <- 2007:2012%>%
  map_dfr(~fread(paste0("result/Data_analysis/Training_dataset/Taiwan_NO2_0.1g_",.x,".csv")))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))

prediction <- lrn_xgboost_tune$predict_newdata(NO2_taiwan)
cv_taiwan <- cv_derive(NO2_taiwan, prediction,raw = FALSE)%>%
  rename(y = truth, ypred = response)%>%
  mutate(level = factor(level, levels = c("Daily","Monthly","Annual")))%>%
  as.data.frame()

cv_taiwan_grid <- cv.plot.facet(cv_taiwan, nrow = 1, cex = 7,
                                alpha = 0.3, facet = "level", label = "", vjust = 1.1)%>%
  plot_set()+
  theme(legend.key.width = unit(.6, "cm"),
        legend.key.height = unit(2, "cm"))

CairoPNG(paste0("result/Model/External_validation_taiwan.png"),
         width = 5000,height = 1800,res=240)
cv_taiwan_grid
dev.off()


##========Spatial CV R2 by sites========
# devtools::install_github("clauswilke/relayer")
library(sf)
source("function/Sf_visualization.R")
# make aesthetics aware size scale, also use better scaling

spplot1_set <- function(base){
  base+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(axis.text = element_text(size = 28),
          legend.title = element_text(colour = "black", size = 38, angle = -90),
          legend.text = element_text(size = 28),
          plot.margin=unit(c(5,0,5,0),"mm"),
          panel.border = element_rect(size = 1,colour = "black", fill = NA),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = alpha("black", 0.2), linewidth = .2))
}
polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

cv_list <- c("sample_CV","site_CV","by_year_CV")
NO2_resample <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))
coords <- 2013:2020%>%
  map_dfr(~fread(paste0("result/Data_analysis/Training_dataset/Training_NO2_0.1g_",.x,".csv")))%>%
  distinct(GridID, X_Lon, Y_Lat)
## RMSE label
breaks = c(-Inf,4,8,12,16,20,+Inf)
labels = label_gen(breaks)
size_cut = c(2,3,4,5,8,10,12)

for (i in 1:length(cv_list)) {
  xgboost_cv <- read_rds(paste0("result/Model/Xgboost_",cv_list[i],".Rds"))
  xgboost_cv_df <- xgboost_cv$prediction()%>%as.data.table()

  xgboost_cv_dfd <- NO2_resample%>%
    as.data.table()%>%
    mutate(row_ids = 1:n())%>%
    dplyr::select(row_ids,GridID,DOY,Month,Year)%>%
    left_join(xgboost_cv_df)%>%left_join(coords)

  xgboost_cv_site <- xgboost_cv_dfd%>%
    group_by(GridID,X_Lon,Y_Lat)%>%
    summarise(R2 = cor(truth, response)^2, 
              RMSE = sqrt(mean((truth - response)^2)))%>%
    ungroup()%>%dplyr::select(X_Lon,Y_Lat,R2, RMSE)%>%
    mutate(RMSE_mod = cut(RMSE, breaks = breaks, labels = labels))%>%
    st_as_sf(coords = c("X_Lon","Y_Lat"),crs = "+proj=longlat +datum=WGS84")%>%
    st_transform(crs = "+proj=aea +lat_0=0 +lon_0=105 +lat_1=25 +lat_2=47 +x_0=0 +y_0=0 +datum=WGS84")
  
  summary(xgboost_cv_site)
  colors <- rev(RColorBrewer::brewer.pal(8,"RdBu"))

  Spatial_CVR2 <- sf_point_plot(xgboost_cv_site, value = "R2")+
    theme_set(theme_bw())+
    geom_sf(data = polygon_sf, fill = "white", linewidth = 0.4)+
    geom_sf(data = Nline_sf, fill = NA, linewidth = 0.5)+
    geom_sf(aes_string(fill = "R2", size = "RMSE_mod"), colour = "white", shape = 21, alpha = .5)+
    scale_size_manual( values = size_cut,
                        name = expression(paste("RMSE (¦Ìg/m"^"3",")")),
                        guide = guide_legend(frame.colour = "black",
                                             title.position="right",
                                             byrow = TRUE,ticks = F,
                                             keywidth = 1.2,keyheight = 1.2,
                                             override.aes = list(fill = colors[8])))+
    scale_fill_gradientn(name = expression(paste("CV R "^"2")),
                         colours = colors,
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth = .5,
                                                title.position="right",
                                                title.vjust = 3,
                                                ticks = F, title.hjust = .5,
                                                barwidth = 3, barheight = 15))+
    theme(legend.title = element_text(colour = "black", size = 30, angle = -90))
  


  CairoPNG(paste0("result/Model/Spatial_CV_accuracy_",cv_list[i],".png"), 
           width = 2500,height = 1500,res = 160)
  print(spplot1_set(Spatial_CVR2))
  dev.off()

}
