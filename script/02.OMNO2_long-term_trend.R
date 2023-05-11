##  This is a subfile of NO2_long-term_trend.R
##  calculate long-term OMNO2 trend
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3mbo)
library(ggplot2)
library(tictoc)
library(Cairo)
library(lubridate)

source("function/Sf_visualization.R")
source("function/Grid_visualization.R")
source("function/Spatial_analysis.R")
source("function/Time_series_analysis.R")

polygon_sf <- st_read(dsn = "function/data", layer = "Province")
Nline_sf <- st_read(dsn = "function/data", layer = "Nline")

SR1 <- st_read(dsn = "data/specific_region/jingjinji", layer = "jingjinji")
SR2 <- st_read(dsn = "data/specific_region/changsanjiao", layer = "changsanjiao")
SR3 <- st_read(dsn = "data/specific_region/zhusanjiao", layer = "zhusanjiao")
OMNO2_annual <- fread("result/Data_analysis/Predictor_annual_mean.csv")%>%
  dplyr::select(X_Lon,Y_Lat,Year,ColumnNO2)

Grid_info <- fread("data/Grid_Information.csv")%>%
  dplyr::select(X_Lon,Y_Lat,GridID)

OMNO2_monthly <- fread("result/OMI_prediction_csv/OMNO2_prediction_0.1g_monthly_2007-2020.csv")%>%
  left_join(Grid_info)
OMNO2_monthly
## NO2 Average in each month
OMNO2_mm <- OMNO2_monthly%>%
  group_by(X_Lon,Y_Lat, Month)%>%
  summarise(OMNO2_mm = mean(ColumnNO2))%>%ungroup()%>%
  left_join(OMNO2_monthly,.)%>%
  mutate(OMNO2_anom = ColumnNO2 - OMNO2_mm)
OMNO2_mm
##========Temporal trend plot: Month anomalies========
## calculate national mean NO2
OMNO2_national_mm <- OMNO2_mm%>%
  mutate(start_month = (Year-2007)*12+Month)%>%
  group_by(start_month)%>%
  summarise(OMNO2_anom = mean(OMNO2_anom))%>%ungroup()%>%
  mutate(region = "China")

## extract grids within specific region
OMNO2_SR_mm <- list(SR1,SR2,SR3)%>%
  map2_dfr(.x = ., .y = c("BTH","YRD","PRD"),
           ~ OMNO2_mm%>%mutate(start_month = (Year-2007)*12+Month)%>%
             dplyr::select(X_Lon,Y_Lat,OMNO2_anom,start_month)%>%
             points_mask(x = "X_Lon", y = "Y_Lat", polygon = .x)%>%
             mutate(region = .y))%>%
  group_by(start_month,region)%>%
  summarise(OMNO2_anom = mean(OMNO2_anom))%>%ungroup()

OMNO2_ts <- full_join(OMNO2_national_mm, OMNO2_SR_mm)%>%
  mutate(period = ifelse(start_month< 72, "2007-2012","2013-2020"))%>%
  full_join(mutate(.,period = "2007-2020"))
fwrite(OMNO2_ts, "result/Time_analysis/Multiyear_OMNO2_Monthly_Anomaly.csv")

OMNO2_ts <- fread("result/Time_analysis/Multiyear_OMNO2_Monthly_Anomaly.csv")%>%
  mutate(period = factor(period, levels = c("2007-2012","2013-2020","2007-2020")),
         region = factor(region, levels = c("China","BTH","YRD","PRD")))
OMNO2_ts

plot_text <- OMNO2_ts%>%
  dplyr::group_by(region, period)%>%
  do(model = broom::tidy(lm(OMNO2_anom ~ start_month, data = .)))%>%
  unnest("model")%>%filter(term == "start_month")%>%
  mutate(signif = ifelse(p.value < 0.005, "***",
                         ifelse(p.value < 0.01, "**",
                                ifelse(p.value < 0.05, "*", ""))))%>%
  mutate(label = paste0(period,": ", expSup(estimate*12),"*\"",
                        signif,"\"*","(",expSup(12*(estimate - 2*std.error)),"-",expSup(12*(estimate + 2*std.error)),")"))%>%
  dplyr::select(region, label, period)%>%
  mutate(x = 5, y = -Inf, hjust = 0,
         vjust = ifelse(period == "2007-2012",-2,
                        ifelse(period == "2013-2020", -1, 0)))

# Time series plot
colors <-  c("#3288bd","#1a9850","#d53e4f")

OMNO2_pplot <- time_series_facet_plot(data = OMNO2_ts,
                                      x = "start_month", y = "OMNO2_anom",
                                      facet = "region", facet_size = 45)+
  scale_x_continuous(name = NULL, expand = c(0.005, 0.005), limits = c(0,168),
                     breaks = seq(0, 168, by = 12), 
                     labels = c(paste0(strrep(" ",9),2007:2020),""))+
  scale_y_continuous(name= expression(paste("NO"[2]," VCD (10"^"15"," mol/cm"^"2",")")), 
                     expand = expansion(mult = c(0.5,0)))+
  theme(plot.title = element_text(size = rel(3), hjust = 0.05, vjust = 1.05),
        plot.title.position = "plot", plot.caption.position =  "plot")

OMNO2_tsplot <- add_line(OMNO2_pplot,x = "start_month", y = "OMNO2_anom", group = "period", 
                         model_colors = colors, size = 1.2)+
  scale_color_manual(values = colors, name = "")+
  geom_text(data  = plot_text, mapping = aes(x = x, y = y, label = label, 
                                             color = period, vjust = vjust, hjust = hjust),
            parse = T, cex = 13)

CairoPNG(paste0("result/Time_analysis/Time_Series_Plot_OMNO2.png"),
         width = 3200, height = 3500, res=140)
print(OMNO2_tsplot)
dev.off()
