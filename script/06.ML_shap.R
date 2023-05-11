## interpretable machine learning
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(tictoc)
library(mlr3)
library(iml)
library(furrr)
library(fastshap)
library(shapper)
library(xgboost)
library(Cairo)
library(shapviz)

## Function for feature name replacement
sname_re <- function(x){
    ifelse(x == "forest_ratio","Forest Land",
    ifelse(x == "constr_ratio","Construction Land",
    ifelse(x == "water_ratio","Water Land",
    ifelse(x == "crop_ratio","Crop Land",
    ifelse(x == "grass_ratio","Grass Land",
    ifelse(x == "unused_ratio","Unused Land", 
    ifelse(x == "ColumnNO2","Column NO2",
    ifelse(x == "blh","Boundary Layer Height",
    ifelse(x == "Pop","Population",
    ifelse(x == "lai_hv"|x == "laihv", "High Vegetation Index",
    ifelse(x == "lai_lv"|x== "lailv","Low Vegetation Index",
    ifelse(x == "t2m","Temperature",
    ifelse(x == "tcw", "Total Column Water",
    ifelse(x == "DOY","Day of Year",
    ifelse(x == "e", "Evaporation",
    ifelse(x == "lsp","Large-scale Precipitation",
    ifelse(x == "sp", "Surface Pressure",
    ifelse(x == "tcc","Total Cloud Cover",
    ifelse(x == "tco3","Total Column Ozone",
    ifelse(x == "u10","U wind",
    ifelse(x == "v10","V wind",x)))))))))))))))))))))
}


fname_re <- function(data, col = NULL, trans_col = FALSE){
  if(trans_col == FALSE){
    data%>%
      mutate_at(all_of(col), .f = sname_re)
  }else{
    new_col <- colnames(data)%>%
      map_chr(.f = sname_re)
    new_data <- setnames(data, new_col)
    return(new_data)
  }
}

SHAP_plot_set <- function(base){
  base+
    theme_bw()+
    xlab(expression(paste("SHAP value: Impact on NO"[2]," concentration"," (¦Ìg/m"^"3",")")))+
    ylab(NULL)+
    scale_color_gradientn(
      name = "Feature Normalized Value\n",
      colours = colors,
      guide = guide_colorbar(frame.colour = "black",
                             frame.linewidth = .5,
                             title.position = "right",ticks = F,
                             title.hjust = .5))+
    theme(panel.background = element_rect(linewidth = 1, fill='white', colour='black'),
          panel.border = element_rect(linewidth = 1, colour = "black", fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(legend.position = "right",
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(3, "cm"),
          legend.text = element_text(colour = "black", size = rel(3.5)),
          legend.title = element_text(colour = "black",size = rel(4), angle = -90),
          axis.text = element_text(colour = "black", size = rel(3.5)),
          axis.title = element_text(colour = "black", size = rel(4)))
  
}

## read in tuning model
lrn_xgboost_tune <- read_rds("result/Model/Final_model.Rds")
modelfn <- lrn_xgboost_tune$model$feature_names

## read in training dataset
# NO2_raw <- fread("Train_Raw/Resampled_Train_Raw_2013-2020.csv")
NO2_resample <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))%>%
  dplyr::select(all_of(modelfn),Site_NO2)

xgboost_im <- xgb.importance(model = lrn_xgboost_tune$model)%>%
  fname_re(col = "Feature")

fwrite(xgboost_im, "result/Interpretation/Feature_Importance.csv")


## combine regional-specific grids during 2007-2020
SR_ID <- fread("Train_Raw/Resampled_Train_Raw_For_Specific_Region.csv")%>%
  distinct(GridID,region)
HS_data <- 2007:2012%>%
  map_dfr(~fread(paste0("Train_Raw/Historical_Sites_Dataset_",.x,".csv")))
Site_data_pre <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  filter(Label == 0)%>%dplyr::select(-Site_NO2,-Label)%>%
  full_join(HS_data)%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))%>%
  left_join(SR_ID)%>%mutate(region = ifelse(is.na(region),"Other",region))

Site_data <- Site_data_pre%>%
  dplyr::select(all_of(modelfn),region)

##======Feature Importance======
xgboost_implot <- xgb.ggplot.importance(importance_matrix = xgboost_im, n_clusters = c(0:1),
                                        measure = "Gain")+
  theme_bw()+
  scale_fill_manual(values = "#de2d26",expand = c(0.01, 0))+
  scale_x_discrete(name = "",expand = c(0.05, 0))+
  scale_y_continuous(name = "Importance (Fractional Contribution)", 
                     expand = c(0.01, 0))+
  labs(title = "")+
  theme(panel.background = element_rect(linewidth = 1,fill='white', colour='black'),
        panel.border = element_rect(linewidth = 1,colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = rel(2)),
        axis.title = element_text(colour = "black", size = rel(2.5)))

CairoPNG("result/Interpretation/Feature_importance.png", width = 2500, height = 2000, res = 150)
xgboost_implot
dev.off()

#======Calculate SHAP value in BTH, YRD and PRD======
# library(SHAPforxgboost)
library(shapviz)
library(Cairo)
library(cowplot)
library(forcats)

## Feature name sorted by importance
fn_imp <- fread("result/Interpretation/Feature_Importance.csv")
fn <- data.frame(fn = lrn_xgboost_tune$model$feature_names)%>%
  fname_re("fn")

lrn_xgboost_tune$model$feature_names <- fn$fn
SR = "BTH"
##**** SHAP calculation****
tic("SHAP value for high-polluted regions. Time-consuming warning in advance")
c("BTH","YRD","PRD")%>%
  purrr::map(function(SR){
    ## change feature names for further plot
    X = Site_data%>%filter(region == SR)%>%dplyr::select(-region)%>%
      fname_re(trans_col = TRUE)

    ## parallel shap calculation
    future::plan(multisession, workers = 18)
    shp_future <- future::future({
      library(shapviz)
      print("check")
      shapviz(lrn_xgboost_tune$model, X_pred = data.matrix(X),
              X, interactions = TRUE, approxcontrib = FALSE)
    })

    shp <- value(shp_future)
    saveRDS(shp, paste0("result/Interpretation/SHAP_value_",SR,".Rds"))
    future:::ClusterRegistry("stop")
  })
toc()

## SHAP summary plot for all predictions (daily level)
tic("Making summary plot")
colors <- rev(RColorBrewer::brewer.pal(8,"RdBu"))
c("BTH","YRD","PRD")%>%
  purrr::map(function(SR){
    shp <- read_rds(paste0("result/Interpretation/SHAP_value_",SR,".Rds"))

    # sv_importance(shp, kind = "bar")
    shap_summary <- shapviz::sv_importance(shp, kind = "beeswarm",
                                  bee_width= .3, bee_adjust = 8,
                                  max_display = 10)+
      aes(y = factor(feature,levels = rev(xgboost_im$Feature)),
          xmin = -20, xmax = 40)
  
    CairoPNG(paste0("result/Interpretation/SHAP_summary_plot_",SR,".png"),
             width = 2500, height = 2000, res = 150)
    print(SHAP_plot_set(shap_summary))
    dev.off()
  })
toc()

## SHAP waterfall plot
shap_waterfall <- shapviz::sv_waterfall(shp, shp$X$Year == 2019)
ggsave(paste0("result/Interpretation/SHAP_waterfall_plot_",SR,".png"),shap_waterfall)

##======Join SHAP in three hotspots======
## Change shp variable names for visualization
SHAP_rename <- function(shp){
  colnames(shp$S) <- colnames(shp$S)%>%data.frame(x = .)%>%
    fname_re(col = "x")%>%.[,1]
  colnames(shp$X) <- colnames(shp$X)%>%data.frame(x = .)%>%
    fname_re(col = "x")%>%.[,1]
  if(!is.null(shp$S_inter)){
    colnames(shp$S_inter) <- colnames(shp$S_inter)%>%data.frame(x = .)%>%
      fname_re(col = "x")%>%.[,1]
    dimnames(shp$S_inter)[[3]] <- dimnames(shp$S_inter)[[3]]%>%data.frame(x = .)%>%
      fname_re(col = "x")%>%.[,1]
  }
  return(shp)
}

## Group SHAP value by joined data frame (e.g. time and space)
#' @param label_df: should only contain the group variable
#' @param shp: SHAP object
#' @param filter: If TRUE add filter argument in...
#' @param interaction: Whether calculate interaction effect
SHAP_group <- function(shp, label_df, filter = FALSE,...){

  if(filter == TRUE){
    SHAP_S <- shp$S%>%cbind(label_df)%>%filter(...)
    SHAP_X <- shp$X%>%cbind(label_df)%>%filter(...)
    group_len <- distinct(label_df)%>%filter(...)%>%nrow
  }else{
    SHAP_S <- shp$S%>%cbind(label_df)
    SHAP_X <- shp$X%>%cbind(label_df)
    group_len <- distinct(label_df)%>%nrow
  }
  
  group_name <- colnames(label_df)
  fn <- colnames(shp$S_inter)
  
  SHAP_S_interg <- array(dim = c(group_len,length(fn),length(fn)),
                         dimnames = list(1:group_len,fn,fn))
  if(!is.null(shp$S_inter)){
    for (i in 1:length(fn)) {
      if(filter == TRUE){
        SHAP_S_interg[,,i] <- shp$S_inter[,,i]%>%cbind(label_df)%>%
          filter(...)%>%
          group_by_at(vars(all_of(group_name)))%>%
          summarise_all(mean)%>%ungroup()%>%
          dplyr::select(-all_of(group_name))%>%as.matrix
      }else{
        SHAP_S_interg[,,i] <- shp$S_inter[,,i]%>%cbind(label_df)%>%
          group_by_at(vars(all_of(group_name)))%>%
          summarise_all(mean)%>%ungroup()%>%
          dplyr::select(-all_of(group_name))%>%as.matrix
      }
    }
  }
  
  ## grouped shp
  shp_group <- shp
  shp_group$S <- SHAP_S%>%
    group_by_at(vars(all_of(group_name)))%>%
    summarise_all(mean)%>%ungroup()%>%
    dplyr::select(-all_of(group_name))%>%as.matrix()
  
  shp_group$X <- SHAP_X%>%
    group_by_at(vars(all_of(group_name)))%>%
    summarise_all(mean)%>%ungroup()%>%
    dplyr::select(-all_of(group_name))%>%as.data.frame()
  if(!is.null(shp$S_inter)){
    shp_group$S_inter <- SHAP_S_interg
  }
  return(shp_group)
}
shp_join <- c("BTH","YRD","PRD")%>%
  purrr::map(~read_rds(paste0("result/Interpretation/SHAP_value_",.x,".Rds")))%>%
  reduce(rbind)%>%
  SHAP_rename()

## Label for shap data frame
Label_index <- Site_data_pre%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time))%>%
  dplyr::select(all_of(modelfn),region,GridID,Year)%>%
  left_join(shp_join$X,.)%>%
  dplyr::select(GridID,Year,region)%>%
  rename(Year_label = Year)
Site_index <- Label_index%>%dplyr::select(GridID,Year_label)

shp_by_site <- SHAP_group(shp_join,Site_index)
shp_by_sitec <- SHAP_group(shp_join,Site_index,filter = TRUE, Year_label>=2013)
shp_by_siteh <- SHAP_group(shp_join,Site_index,filter = TRUE, Year_label<2013)

shp_by_site$X$Population <- shp_by_site$X$Population/1e6
shp_by_sitec$X$Population <- shp_by_sitec$X$Population/1e6
shp_by_siteh$X$Population <- shp_by_siteh$X$Population/1e6

## Summary plot of three hotspots
colors <- rev(RColorBrewer::brewer.pal(8,"RdBu"))
shap_summary_site <- shapviz::sv_importance(shp_by_site, kind = "beeswarm",
                                       bee_width= .3, bee_adjust = 8,
                                       max_display = 10)+
  aes(y = factor(feature,levels = rev(xgboost_im$Feature)),
      xmin = -10, xmax = 20)
ggsave(paste0("result/Interpretation/SHAP_summary_plot_site.png"),
       SHAP_plot_set(shap_summary_site),width = 3000, height = 2000,dpi = 150, units = "px")


list(shp_by_site,shp_by_sitec,shp_by_siteh)%>%
  map2(.x =., .y = c("site","sitec","siteh"),
       function(x,y){
        shap_contribution <- shapviz::sv_importance(x, max_display = 10)+
           aes(xmin = -10, xmax = 20)+theme_bw()
         
         ggsave(paste0("result/Interpretation/SHAP_contribution_plot_",y,".png"),
                shap_contribution)
       })

SHAP_depplot_set <- function(base,title.position = "top",...){
  base+
    theme_bw()+ylab("SHAP Value")+
    scale_color_gradientn(
      colours = colors,
      guide = guide_colorbar(frame.colour = "black",
                             frame.linewidth = .5,
                             title.position = title.position,ticks = F),...)+
    theme(panel.background = element_rect(linewidth = 1.2, fill='white', colour='black'),
          panel.border = element_rect(linewidth = 1.2, colour = "black", fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(legend.position = "right",
          legend.key.width = unit(.8, "cm"),
          legend.key.height = unit(2, "cm"),
          legend.text = element_text(colour = "black", size = rel(2)),
          legend.title = element_text(colour = "black", hjust = .5,
                                       size = rel(2.5)),
          axis.text = element_text(colour = "black", size = rel(2)),
          axis.title = element_text(colour = "black", size = rel(2.5)))
}

## Main dependence plot
colors <- rev(RColorBrewer::brewer.pal(8,"PuOr"))
list(shp_by_site,shp_by_sitec,shp_by_siteh)%>%
  map2(.x =., .y = c("site","sitec","siteh"),
       function(x,y){
         shap_dep_site <- c("Column NO2","Boundary Layer Height",
                            "Construction Land", "Population")%>%
           purrr::map(~ sv_dependence(x, v = .x, color_var = "Column NO2")%>%
                        SHAP_depplot_set()+
                        geom_smooth(color = "red")
           )%>%patchwork::wrap_plots(nrow = 2)
         
         ggsave(paste0("result/Interpretation/SHAP_dependence_plot_",y,".png"),
                shap_dep_site,width = 2500, height = 2000,dpi = 150, units = "px")
       })

## Dependence plot with other meteorological variables
list(shp_by_site,shp_by_sitec,shp_by_siteh)%>%
  map2(.x =., .y = c("site","sitec","siteh"),
       function(x,y){
         shap_dep_site <- c("Temperature","Total Column Water","U wind","V wind")%>%
           purrr::map(~ sv_dependence(x, v = .x, color_var = "Column NO2")%>%
                        SHAP_depplot_set()+
                        geom_smooth(color = "red")
           )%>%patchwork::wrap_plots(nrow = 2)
         
         ggsave(paste0("result/Interpretation/SHAP_dependence_plot_add_",y,".png"),
                shap_dep_site,width = 2500, height = 2000,dpi = 150, units = "px")
       })



## Interaction effect
## Interaction effect of ColumnNO2
shapviz::sv_interaction(shp_by_site, kind = "no")[,1]%>%
  sort
colors <- rev(RColorBrewer::brewer.pal(8,"PuOr"))
shap_dep_inter <- c("Day of Year","Boundary Layer Height",
                    "Construction Land","Population",
                    "Total Column Water","Temperature")%>%
  purrr::map(function(x){
    if(x == "Year"){
      sv_dependence(shp_by_site, v = "Column NO2", 
                    color_var = x, interactions = TRUE)%>%
        SHAP_depplot_set(title.position = "right",breaks = c(2007,2020))+
        geom_smooth(color = "red")+
        ylab(NULL)+
        theme(legend.justification = "left",
              legend.title = element_text(colour = "black", size = rel(2), angle = -90))
    }else{
      sv_dependence(shp_by_site, v = "Column NO2", 
                    color_var = x, interactions = TRUE)%>%
        SHAP_depplot_set(title.position = "right")+
        geom_smooth(color = "red")+
        ylab(NULL)+
        theme(legend.justification = "left",
              legend.title = element_text(colour = "black", size = rel(2.5), angle = -90))
      
    }
  })%>%patchwork::wrap_plots(nrow = 2, guides = 'auto')

ggsave(paste0("result/Interpretation/SHAP_dep-interaction_plot.png"),
       shap_dep_inter,width = 4000, height = 2000,dpi = 150, units = "px")


SHAP_interplot_set <- function(base){
  SHAP_plot_set(base)+
    xlab(expression(paste("SHAP value: Impact on NO"[2]," concentration"," (¦Ìg/m"^"3",")")))+
    aes(Variable1 = factor(Variable1,levels = rev(xgboost_im$Feature)),
        Variable2 = factor(Variable1,levels = rev(xgboost_im$Feature)))+
    theme(
      strip.background = element_rect(linewidth = 1.5, fill = "white",colour='black'),
      strip.text.x = element_text(size = 30, color = "black",margin = margin(.4, .1, .4, .1, "cm")),
      strip.text.y = element_text(size = 30, color = "black",margin = margin(.1, .4, .4, .1, "cm")),
      axis.text.x = element_text(size = 30, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_text(size = 40),
      legend.key.width = unit(.8, "cm"),
      legend.key.height = unit(3, "cm"),
      legend.text = element_text(colour = "black", size = 35),
      legend.title = element_text(colour = "black", vjust = .5, size = 45))
  
}
colors <- rev(RColorBrewer::brewer.pal(8,"RdBu"))
VOI <- c("Column NO2","Boundary Layer Height","Construction Land","Population")
shp_by_inter <- shp_by_site
shp_by_inter$S_inter <- shp_by_site$S_inter[,VOI,VOI]
shp_by_inter$X <- shp_by_site$X%>%dplyr::select(VOI)
shp_by_inter$S <- shp_by_site$S[,VOI]

shap_inter <- sv_interaction(shp_by_inter, alpha = .4,
                             max_display = 4, stroke = .1, size = 2)%>%
  SHAP_interplot_set()

CairoPNG(paste0("result/Interpretation/SHAP_interaction_plot.png"),
         width = 3500, height = 3200, res = 150)
print(shap_inter)
dev.off()


VOI <- c("Temperature","Total Column Water","U wind","V wind")
shp_by_inter <- shp_by_site
shp_by_inter$S_inter <- shp_by_site$S_inter[,VOI,VOI]
shp_by_inter$X <- shp_by_site$X%>%dplyr::select(VOI)
shp_by_inter$S <- shp_by_site$S[,VOI]
shap_inter <- sv_interaction(shp_by_inter, alpha = .4,
                             max_display = 10, stroke = .1, size = 3)%>%
  SHAP_interplot_set()

CairoPNG(paste0("result/Interpretation/SHAP_interaction_plot_add.png"),
         width = 3200, height = 2800, res = 150)
print(shap_inter)
dev.off()

##======Interaction data frame======
Interaction_df <- xgboost_im$Feature%>%
  map_dfr(~ shp_by_site$S_inter[,,.x]%>%
            as.data.frame()%>%
            dplyr::select(xgboost_im$Feature)%>%
            summarise_all(mean)%>%ungroup()%>%
            mutate(Feature = .x, .before = 1))

fwrite(Interaction_df, "result/Interpretation/SHAP_feature_interaction.csv")
IF <- 1:nrow(Interaction_df)%>%
  map_chr(~ ifelse(Interaction_df[.x,.x+1] > 0.1, Interaction_df[.x,"Feature"], NA))%>%
  na.omit()

Interaction_IFdf <- Interaction_df%>%
  filter(Feature%in%IF)%>%
  dplyr::select(Feature,all_of(IF))
Interaction_IFdf
fwrite(Interaction_IFdf, "result/Interpretation/SHAP_importance_feature_interaction.csv")


## validation (absolute value)
## interaction part is multiplied by 2
sv_interaction(shp_by_site, kind = "no")[1,]
Itest <- xgboost_im$Feature%>%
  map_dfr(~ shp_by_site$S_inter[,,.x]%>%
            as.data.frame()%>%
            summarise_all(.funs = function(x){mean(abs(x))})%>%ungroup()%>%
            mutate(Feature = .x, .before = 1))
Itest[1,]
Itest[2,]


(Interaction_IFdf[3,])
Mtest <- shp_by_site$S%>%
  as.data.frame()%>%
  summarise_all(mean)
Mtest
Interaction_df[1,]%>%dplyr::select(-Feature)%>%sum
sv_importance(shp_by_site, kind = "no")

##======Main SHAP by Year======
shp_by_year <- SHAP_group(shp_join,Label_index%>%dplyr::select(Year_label))
Time = shp_by_year$X$Year
SHAP_S <- shp_by_year$S%>%
  as.data.frame()%>%
  mutate(Year_label = Time, .before = 1)%>%
  pivot_longer(cols = -Year_label)%>%
  pivot_wider(names_from = "Year_label",values_from = "value")%>%
  as.data.frame()
SHAP_X <- shp_by_year$X%>%mutate(Year_label = Time, .before = 1)%>%
  pivot_longer(cols = -Year_label)%>%
  pivot_wider(names_from = "Year_label",values_from = "value")%>%
  as.data.frame()
SHAP_df <- full_join(SHAP_S,SHAP_X)
fwrite(SHAP_df, "result/Interpretation/SHAP_by_year.csv")

##=====Comparison with OLS model======
## Some variables are defined in ML SHAP part
library(stargazer)
.min_max_scale <- function(z, na.rm = TRUE) {
  r <- range(z, na.rm = na.rm)
  d <- diff(r)
  if (is.na(d) || d == 0) {
    z[!is.na(z)] <- 0.5
    return(z)
  }
  (z - r[1L]) /(r[2L] - r[1L])
}

.scale_X <- function(X) {
  X_scaled <- apply(data.matrix(X), 2L, FUN = .min_max_scale)
  X_scaled%>%as.vector()
}

# modelfnsel <- c("ColumnNO2","Pop","constr_ratio","blh","t2m","u10","v10","tcw")
LM_data <- fread(paste0("Train_Raw/Resampled_Train_Raw_2013-2020.csv"))%>%
  mutate(Month = month(as.Date("2013-01-01")+Time),
         DOY = yday(as.Date("2013-01-01")+Time),
         Pop = Pop/1e6)%>%
  dplyr::select(all_of(modelfn),Site_NO2)

sociofn <- c("constr_ratio","crop_ratio","forest_ratio","grass_ratio",
             "water_ratio","lai_hv","lai_lv","Pop")
timefn <- c("DOY","Month","Year")
metfn <- setdiff(modelfn,c(sociofn,timefn,"ColumnNO2"))
modelfn <- c("ColumnNO2",sociofn,timefn,metfn)

## Main regression
fm <- paste0("Site_NO2 ~",paste(modelfn,collapse = "+"))
fm1 <- paste0("Site_NO2 ~",paste(setdiff(modelfn,timefn),collapse = "+"))
fm2 <- paste0("Site_NO2 ~",paste(setdiff(modelfn,sociofn),collapse = "+"))
fm3 <-  paste0("Site_NO2 ~",paste(setdiff(modelfn,metfn),collapse = "+"))

LM_model <- lm(fm, LM_data)
LM_model1 <- lm(fm1, LM_data)
LM_model2 <- lm(fm2, LM_data)
LM_model3 <- lm(fm3, LM_data)


stargazer(LM_model, LM_model1, LM_model2, LM_model3,
          type = "html", single.row = TRUE, out = "result/Interpretation/LM_model.html")

LM_shp <- c("BTH","YRD","PRD")%>%
  purrr::map(function(SR){
    ## change feature names for further plot
    X <- Site_data%>%filter(region == SR)%>%dplyr::select(-region)%>%
      mutate(Pop = Pop/1e6)%>%
      dplyr::select(modelfn)
    
    shp_SR <- fastshap::explain(LM_model, X = X,
                      newdata = X, exact = TRUE)%>%
      shapviz(X = X, baseline = shp_join$baseline)
    
    return(shp_SR)
  })%>%reduce(rbind)

LM_site_shp <- SHAP_group(LM_shp,Site_index)%>%
  SHAP_rename()

## Summary plot for OLS model
colors <- rev(RColorBrewer::brewer.pal(8,"RdBu"))
LMshap_summary_site <- shapviz::sv_importance(LM_site_shp, kind = "beeswarm",
                                            bee_width= .3, bee_adjust = 8,
                                            max_display = 10)+
  aes(y = factor(feature,levels = rev(xgboost_im$Feature)),
      xmin = -10, xmax = 20)
ggsave(paste0("result/Interpretation/OLS_SHAP_summary_plot_site.png"),
       SHAP_plot_set(LMshap_summary_site),width = 3000, height = 2000,dpi = 150, units = "px")

##  Dependence plot
colors <- rev(RColorBrewer::brewer.pal(8,"PuOr"))
LMshap_dep_site <- c("Column NO2","Boundary Layer Height","Construction Land","Population")%>%
  purrr::map(~ sv_dependence(LM_site_shp, v = .x, color_var = "Column NO2")%>%
               SHAP_depplot_set()
  )%>%patchwork::wrap_plots(nrow = 2)

ggsave(paste0("result/Interpretation/OLS_SHAP_dependence_plot_site.png"),
       LMshap_dep_site,width = 2500, height = 2000,dpi = 150, units = "px")


list(LM_model,LM_model1,LM_model2,LM_model3)%>%
  purrr::map(function(model){
    
    LM_shp <- fastshap::explain(model, X = SR_data,
                                newdata = SR_data, exact = TRUE)
  })

## OLS SHAP Tabel
LM_shp_mean <- LM_site_shp$S%>%
  as.data.frame()%>%
  pivot_longer(everything())%>%
  group_by(name)%>%
  summarise(OLS = mean(value), OLS_sd = sd(value))%>%
  ungroup()

ML_shp_mean <- shp_by_site$S%>%
  as.data.frame()%>%
  pivot_longer(everything())%>%
  group_by(name)%>%
  summarise(XGBoost = mean(value), XGBoost_sd = sd(value))%>%
  ungroup()
Average_SHAP <- left_join(ML_shp_mean, LM_shp_mean)%>%
  arrange(name)
Average_SHAP
fwrite(Average_SHAP, "result/Interpretation/Average_SHAP_comparison.csv")

LM_shp%>%as.data.frame()%>%
  mutate(Year_Label = (SR_data$Year))%>%
  group_by(Year_Label)%>%
  summarise_all(mean)%>%ungroup()%>%as.data.frame()

LM_shp%>%
  mutate(GridID = LM_SR$GridID)%>%
  group_by(GridID,Year)%>%
  summarise_all(mean)%>%ungroup()
