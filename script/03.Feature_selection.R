setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)

##======Read in dataset======
year = 2013
# fread(paste0("result/Data_analyais/Training_dataset/Traing_NO2_0.1g_",year,".csv"))
NO2_resample <- 2013:2020%>%
  map_dfr(~ fread(paste0("Train_Raw/Resampled_Train_Raw_",.x,".csv"))%>%
            dplyr::select(-GridID,-X_Lon,-Y_Lat,-Label))
colnames(NO2_resample)

NO2_task <- TaskRegr$new("NO2_multiyear_test", NO2_resample, target = "Site_NO2")
# split <- partition(task_NO2, ratio = 0.9)
NO2_task$data()
NO2_task$filter(which(complete.cases(NO2_task$data())))

#======Feature selection======
library(correlation)
feature_selection <- function(task, lrn){
  ## filter by information gain
  filter = flt("information_gain")
  filter$calculate(task)
  information_gain <- as.data.table(filter)%>%
    rename(information_gain = score)%>%
    mutate(IG_percent = information_gain/sum(information_gain, na.rm = T)*100)
  as.data.table(filter)[,"feature"]
  
  ## filter by pearson correlation

  cor_df <- correlation(task$data()%>%dplyr::select(task$target_names),
                        task$data()%>%dplyr::select(-task$target_names),
                        method = "pearson")%>%as.data.table()%>%
    mutate(r2 = r^2)%>%
    arrange(desc(r2),p)%>%
    rename(feature = Parameter2)
  
  ## filter by variable importance
  filter_lrn = flt("importance", learner = lrn)
  filter_lrn$calculate(task)

  importance_df <- as.data.table(filter_lrn)%>%
    rename(importance = score)%>%
    mutate(import_percent = importance/sum(importance, na.rm = T)*100)
  
  feature_sel <- left_join(cor_df,information_gain)%>%
    left_join(importance_df)
  return(feature_sel)
}

# fwrite(information_gain, "result/Data_processing/Feature_selection_information_gain.csv")

# fwrite(cor_df, "result/Data_processing/Feature_selection_pearson.csv")

## Prior feature selection
lrn_xgboost <- lrn("regr.xgboost",objective="reg:squarederror", 
                   eval_metric="rmse",
                   booster = "gbtree",
                   nthread = 10,
                   nrounds=500L, eta=0.1, gamma=0.1)
tic("Feature selection")
feature_sel <- feature_selection(NO2_task, lrn_xgboost)
toc()
fwrite(feature_sel, "result/Model/Feature_selection.csv")

fwrite(cor(NO2_task$data())%>%as.data.frame(), 
       "result/Model/Correlation_martix.csv", row.names = TRUE)

## After arbitrary selection
rmv <- c("TQI","TQL","TROPPB","TO3","U2M","V2M","Month","GridID")
feature_clean <- setdiff(NO2_task$feature_names,rmv)

NO2_task_feature <- NO2_task$clone()$select(feature_clean)

tic("Feature selection")
feature_sel_after <- feature_selection(NO2_task_feature, lrn_xgboost)
toc()
fwrite(feature_sel_after, "result/Model/Feature_selection_afterwards.csv")


##=====xgboost training======
lrn_xgboost <- lrn("regr.xgboost",objective="reg:squarederror", 
                   eval_metric="rmse",
                   booster = "gbtree",
                   nthread = 10,
                   nrounds=500L, eta=0.1, gamma=0.1)
# lrn_xgboost <- makeLearner("regr.xgboost",predict.type = "response")
lrn_xgboost$train(NO2_task)
lrn_xgboost$importance()
xgboost::xgb.importance(model = lrn_xgboost$model)

predictions <- lrn_xgboost$predict(NO2_task)
measure = msrs(c("regr.rsq","regr.rmse"))
predictions$score(measure)


p1 <- read_rds("result/Model/Xgboost_sample_CV.Rds")
p2 <- read_rds("result/Model/Xgboost_site_CV.Rds")
p3 <- read_rds("result/Model/Xgboost_by_year_CV.Rds")

measure <- msrs(c("regr.rsq","regr.rmse"))
p1$prediction()$score(measure)
p2$prediction()$score(measure)
p3$prediction()$score(measure)
