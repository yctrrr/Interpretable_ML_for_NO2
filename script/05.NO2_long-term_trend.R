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

Pop <- fread("data/Urban_Pop_0.1g_2000-2020.csv")%>%
  rename(X_Lon = x, Y_Lat = y,
         Pop = Pop_all, UPop = Pop_U, RPop = Pop_O)%>%
  dplyr::select(X_Lon,Y_Lat,Year,Pop,UPop,RPop)

NO2_monthly <- fread(paste0("result/NO2_prediction_csv/NO2_prediction_0.1g_monthly_2007-2020.csv"))%>%
  left_join(Pop)%>%na.omit()

## NO2 Average in each month
NO2_mm <- NO2_monthly%>%
  group_by(X_Lon, Y_Lat, month)%>%
  summarise(NO2_mm = mean(NO2_response))%>%ungroup()%>%
  left_join(NO2_monthly,.)%>%
  mutate(NO2_anom = NO2_response - NO2_mm)

##========Temporal trend plot: Month anomalies========
## calculate national mean NO2 (already done)
# NO2_national_mm <- NO2_mm%>%
#   mutate(start_month = (Year-2007)*12+month)%>%
#   group_by(start_month)%>%
#   summarise(NO2_anom = mean(NO2_anom))%>%ungroup()%>%
#   mutate(region = "China")
# 
# ## extract grids within specific region
# NO2_SR_mm <- list(SR1,SR2,SR3)%>%
#   map2_dfr(.x = ., .y = c("BTH","YRD","PRD"),
#            ~ NO2_mm%>%mutate(start_month = (Year-2007)*12+month)%>%
#              dplyr::select(X_Lon,Y_Lat,NO2_anom,start_month)%>%
#              points_mask(x = "X_Lon", y = "Y_Lat", polygon = .x)%>%
#              mutate(region = .y))%>%
#   group_by(start_month,region)%>%
#   summarise(NO2_anom = mean(NO2_anom))%>%ungroup()

# NO2_ts <- full_join(NO2_national_mm, NO2_SR_mm)%>%
#   mutate(period = ifelse(start_month< 72, "2007-2012","2013-2020"))%>%
#   full_join(mutate(.,period = "2007-2020"))
# fwrite(NO2_ts, "result/Time_analysis/Multiyear_NO2_Monthly_Anomaly.csv")

NO2_ts <- fread("result/Time_analysis/Multiyear_NO2_Monthly_Anomaly.csv")%>%
  mutate(period = factor(period, levels = c("2007-2012","2013-2020","2007-2020")),
         region = factor(region, levels = c("China","BTH","YRD","PRD")))
NO2_ts
plot_text <- NO2_ts%>%
  dplyr::group_by(region, period)%>%
  do(model = broom::tidy(lm(NO2_anom ~ start_month, data = .)))%>%
  unnest("model")%>%filter(term == "start_month")%>%
  mutate(signif = ifelse(p.value < 0.005, "***",
                  ifelse(p.value < 0.01, "**",
                  ifelse(p.value < 0.05, "*", ""))))%>%
  mutate(label = paste0(period,": ", expSup(estimate*12),"*\"",signif,"\"*","(",expSup(12*(estimate - 2*std.error)),"-",expSup(12*(estimate + 2*std.error)),")"))%>%
  dplyr::select(region, label, period)%>%
  mutate(x = 5, y = -Inf, hjust = 0,
         vjust = ifelse(period == "2007-2012",-2,
                 ifelse(period == "2013-2020", -1, 0)))

# Time series plot
colors <-  c("#3288bd","#1a9850","#d53e4f")

NO2_pplot <- time_series_facet_plot(data = NO2_ts,
                       x = "start_month", y = "NO2_anom", 
                       facet = "region", facet_size = 45)+
  scale_x_continuous(name = NULL, expand = c(0.005, 0.005), limits = c(0,168),
                     breaks = seq(0, 168, by = 12), 
                     labels = c(paste0(strrep(" ",9),2007:2020),""))+
  scale_y_continuous(name= expression(paste("NO"[2]," anomalies"," (¦Ìg/m"^"3",")")), 
                     expand = expansion(mult = c(0.3,0)))+
  theme(plot.title = element_text(size = rel(3), hjust = 0.05, vjust = 1.05),
        plot.title.position = "plot", plot.caption.position =  "plot")

NO2_tsplot <- add_line(NO2_pplot,x = "start_month", y = "NO2_anom", group = "period", 
                       model_colors = colors, size = 1.2)+
  scale_color_manual(values = colors, name = "")+
  geom_text(data  = plot_text, mapping = aes(x = x, y = y, label = label, 
                                             color = period, vjust = vjust, hjust = hjust),
            parse = T, cex = 13)

CairoPNG(paste0("result/Time_analysis/Time_Series_Plot.png"),
         width = 3200, height = 3500, res=140)
print(NO2_tsplot)
dev.off()

##========Trend plot at grid level========
spplot_set <- function(base){
  base+
    geom_sf(data = polygon_sf, fill = NA, linewidth = 0.2)+
    geom_sf(data = Nline_sf, fill = NA, linewidth = 0.5)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(strip.text.x = element_text(
      size = 35, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
      strip.background = element_rect(color = NA, fill = NA))+
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

period_lst <- list(c(2007:2012),c(2013:2020),c(2007:2020))
period_name <- c("2007-2012","2013-2020","2007-2020")
## Temporal trend calculated at grid level
# NO2_trend <- 1:length(period_lst)%>%
#   map_dfr(~ NO2_mm%>%
#             mutate(start_month = (Year-2007)*12+month)%>%
#             filter(Year%in%period_lst[[.x]])%>%
#             dplyr::group_by(X_Lon,Y_Lat)%>%
#             do(model = broom::tidy(lm(NO2_anom ~ start_month, data = .)))%>%
#             unnest("model")%>%filter(term == "start_month")%>%
#             dplyr::select(X_Lon,Y_Lat,estimate, std.error, p.value)%>%
#             mutate(period = period_name[.x]))
# 
# fwrite(NO2_trend,"result/Time_analysis/Multiyear_NO2_Trend.csv")

NO2_trend <- fread("result/Time_analysis/Multiyear_NO2_Trend.csv")
breaks <- c(-Inf,-1,-0.5,-0.1,-0.01,0.01,0.1,0.5,1,Inf)
labels <- label_gen(breaks)

## Label the trend estimate
NO2_trend_mod <- period_name%>%
  purrr::map_dfr(~NO2_trend%>%filter(period == .x)%>%
            dplyr::select(-period)%>%
            dfproj(method = "bilinear")%>%
            mutate(period = .x))%>%
  mutate(signif= factor(ifelse(p.value < 0.01, "p < 0.01",
                        ifelse(p.value < 0.05, "p < 0.05",
                        ifelse(p.value < 0.1, "p < 0.1","p >= 0.1"))),
                        levels = c("p < 0.01", "p < 0.05", "p < 0.1", "p >= 0.1")),
         trend = cut(estimate*12 , breaks = breaks, labels = labels, right = TRUE),
         period = factor(period, levels = period_name))

colors <- RColorBrewer::brewer.pal(length(labels),"RdYlBu")%>%rev()

NO2_trend_plot <- grid_facet_plot(NO2_trend_mod, value = "trend",
                                  facet = "period", facet_text = 35)+
  scale_fill_manual(name = expression(paste("NO"[2]," trend (¦Ìg/m"^"3","/year)")),
                    values = colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))

CairoPNG(paste0("result/Time_analysis/Multiyear_NO2_Trend.png"),
         width = 3500,height = 1500,res=160)
print(spplot_set(NO2_trend_plot))
dev.off()

colors <- colorspace::sequential_hcl(4, palette = "Purp")

NO2_signif_plot <- grid_facet_plot(NO2_trend_mod, value = "signif", 
                                   facet = "period", facet_text = 35)+
  scale_fill_manual(name = expression("Statistical significance"),
                    values =  colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))
CairoPNG(paste0("result/Time_analysis/Multiyear_NO2_Trend_Significance.png"),
         width = 3500,height = 1500,res=160)
print(spplot_set(NO2_signif_plot))
dev.off()

