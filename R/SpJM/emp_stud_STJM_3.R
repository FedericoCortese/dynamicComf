library(dplyr)
library(lubridate) 
library(tidyr)
library(tidyverse)
library(shiny)
library(geosphere)
library(caret)
library(leaflet)
library(reshape2)
library(scales)
# Set Singapore time zone
Sys.setenv(TZ="GMT+08")
##
utci_sing=read.table("utci_202304_17-27.txt",header=T)

#dalle 0 del 17 aprile alle 23 del 27 aprile
start_utci=as.POSIXct("2023-04-17 01:00:00",tz="UTC")
end_utci=as.POSIXct("2023-04-27 23:00:00",tz="UTC")

# Construct the time vector
time_utci=seq(from=start_utci,to=end_utci,by="hour")

# Singapore time zone
time_utci_8=format(time_utci, tz="Singapore",usetz=TRUE)
time_utci_8=as.POSIXct(time_utci_8,format="%Y-%m-%d %H:%M:%S",tz="Singapore")
data_utci=data.frame(time=time_utci_8,UTCI=utci_sing[,1]-273.15)

####

# SEE emp_stud_STJM_2 for Y_3.Rdata
#load("D:/git/dynamicComf/dynamicComf/R/SpJM/Y_3.Rdata")
load("Y_3.Rdata")
load("data_utci.Rdata")

# Merge with pre-existing data
Y_4=Y_complete[,c("t","m","time")]%>%
  left_join(Y_3,by=c("t","m"))

Y_5=Y_4%>%
  left_join(data_utci,by="time")

Y_6=Y_5[-(1:56),]
Y_6=subset(Y_6,select=-time)
###

load("Y_6.Rdata")

# leaflet map singapore ---------------------------------------------------

leaflet(data_stat_number, options = leafletOptions(zoomControl = FALSE)) %>%
  addTiles() %>%  # Add default OpenStreetMap tiles
  addCircleMarkers(
    ~longitude, ~latitude,  # Specify longitude and latitude
    radius = 6,  # Set circle size
    color = "red",  # Outline color of the circle
    fillColor = "red",  # Fill color of the circle
    fillOpacity = 0.6,  # Circle fill opacity
    label = ~paste0("S", m),  # Use m to label the stations as "S1", "S2", ...
    labelOptions = labelOptions(noHide = TRUE, textsize = "12px", direction = "top", style = list("color" = "red"))  # Set label options
  ) %>%
  addMarkers(
    ~longitude, ~latitude, 
    label = ~paste0("S", m), 
    popup = ~paste("Station:", station), 
    labelOptions = labelOptions(noHide = TRUE, textsize = "12px", direction = "top", style = list("color" = "red"))
  ) %>%
  setView(lng = mean(data_stat_number$longitude), lat = mean(data_stat_number$latitude), zoom = 12) %>%
  addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE))  # Add a scale bar in kilometers (metric)

summary(Y_6)
round(cor(Y_6[complete.cases(Y_6),c(3:6,17)]),2)
library(ppcor)
pcor(Y_6[complete.cases(Y_6), c(3:6,17)])
Y_6%>%group_by(as.factor(m))%>%
  select(air_temp,rh,rainfall,wind_speed)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# STJM fit ----------------------------------------------------------------

source("Utils.R")
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

#lambda=0.05
#gamma=0.05

lambda=0.05
gamma=0.05

fit=STjumpDist(Y_6,3,D,
               jump_penalty=lambda,
               spatial_penalty=gamma,
               initial_states=NULL,
               max_iter=10, n_init=10, 
               tol=1e-4, 
               verbose=T,timeflag=T)

# State characterization
M=length(unique(Y_6$m))
State=c(t(fit$best_s))
State=order_states_condMean(Y_6$air_temp,State)

S_est=matrix(State,ncol=M,byrow = T)

wdn=4
times=Y_5$time[-(1:56)]
Y_res=data.frame(Y_6,State,times)
table(Y_res$State)/length(times)*100

# tapply(Y_res$air_temp,Y_res$State,mean,na.rm=T)
# tapply(Y_res$rh,Y_res$State,mean,na.rm=T)
# tapply(Y_res$rainfall,Y_res$State,mean,na.rm=T)
# tapply(Y_res$wind_speed,Y_res$State,mean,na.rm=T)
# tapply(Y_res$UTCI,Y_res$State,mean,na.rm=T)
tapply(Y_res$air_temp,Y_res$State,median,na.rm=T)
tapply(Y_res$rh,Y_res$State,median,na.rm=T)
tapply(Y_res$rainfall,Y_res$State,median,na.rm=T)
tapply(Y_res$wind_speed,Y_res$State,median,na.rm=T)
tapply(Y_res$UTCI,Y_res$State,median,na.rm=T)

tapply(Y_res$windy,Y_res$State,Mode)
tapply(Y_res$hour,Y_res$State,Mode)

TY=unique(Y_6$t)
timesY=unique(times)
TT_6=length(TY)

#load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/emp_stud_STJM_singapore/cozie_compare.Rdata")

# Strano che ws_timestamp_start ha SOLO orari dalle 2 di notte alle 11 di mattina, forse meglio prendere time
q_data=cozie_compare[,c("q_thermal_preference","time","ws_longitude","ws_latitude")]
#q_data=cozie_compare[,c("q_thermal_preference","ws_timestamp_start","ws_longitude","ws_latitude")]
colnames(q_data)[2:4]=c("time","longitude","latitude")
# wdn="1 hour"
# q_data$time=round_date(q_data$time,wdn)
q_data$State=q_data$q_thermal_preference
q_data$State=recode(q_data$State, "Warmer" = 1, "No change" = 2, "Cooler" = 3)

cozie_dist=distm(q_data[,c("longitude","latitude")],
                 data_stat_number[,c("longitude","latitude")], 
                 fun = distGeo)/1000

summary(cozie_dist)

# Radius for the neighbourhood
R=max(apply(cozie_dist, 1, min ))

# Which stations are in the neighbourhood of each feedback location
neigh=cozie_dist<=R
colnames(neigh)=loc_weath$id
neigh=data.frame(time=q_data$time,neigh)
colnames(neigh)[2:ncol(neigh)]=data_stat_number$m
sum(apply(neigh[,-1],1,sum))

S_est=data.frame(time=timesY,S_est)
colnames(S_est)[2:ncol(S_est)]=paste("S", data_stat_number$m)
indx=left_join(neigh,S_est,by="time")
indx=indx[,2:15]*indx[,16:29]
indx[indx==0]=NA

modes=factor(apply(indx,1,Mode,na.rm=T),levels=1:3)
comparison=factor(q_data$State)
confusionMatrix(modes,comparison)


# state plot --------------------------------------------------------------

## COMFORT LEVEL NON VA BENE, DOBBIAMO TROVARE UN NOME MENO SPECIALIZZATO

long_df <- S_est %>%
  pivot_longer(cols = starts_with("S "), names_to = "Station", values_to = "ComfortLevel")

# Step 2: Calculate the proportion of each ComfortLevel for each time point
S_summary <- long_df %>%
  group_by(time, ComfortLevel) %>%
  summarise(Proportion = n() / n_distinct(long_df$Station), .groups = 'drop')  # Calculate proportion correctly

# Step 3: Extract hour from the time column
S_summary <- S_summary %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S"),
         Hour = hour(time))  # Extract hour as integer

S_summary$ComfortLevel <- factor(S_summary$ComfortLevel, levels = c("1", "2", "3"))

S_summary_complete <- S_summary %>%
  complete(time, ComfortLevel = factor(1:3), fill = list(Proportion = 0))

S_summary_complete <- S_summary_complete %>%
  mutate(Hour = hour(time))

avg_weath <- Y_res%>%
  group_by(times)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Create the initial plot object for the area plot
plot_base <- ggplot(S_summary_complete, aes(x = time, y = Proportion, fill = as.factor(ComfortLevel))) +
  geom_area(alpha = 0.6, size = 1, color = "black") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime",
                    labels=c("Cool","Neutral","Hot")) +
  labs(
    #title = "Temporal Evolution of Thermal Comfort Levels and Wind Speed",
    # x = "Date and Hour",
    x=" ",
    y = "Prop of Locations",
    fill = "Comfort Regime") +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%Y-%m-%d %H:%M") +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.8, "cm"),  # controls the box size
    legend.text = element_text(size = 10),  # controls the text size
    text = element_text(size = sz),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y.right = element_blank(),
    axis.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "top"
  )
  # theme(text = element_text(size = 12),
  #       axis.text.x = element_text(angle = 45, hjust = 1),
  #       legend.position = "top")
# + 
library(scales)
temp_state_plot <- plot_base +
  geom_line(
    data = avg_weath,
    aes(x = times, y = rescale(air_temp, to = c(0, 1))),  # Riscalata tra 0 e 1
    color = "red", size = 1.2, inherit.aes = FALSE, linetype = 6
  ) +
  scale_y_continuous(
    name = "Prop of Locations",
    sec.axis = sec_axis(
      ~ rescale(., from = c(0, 1), to = range(avg_weath$air_temp, na.rm = TRUE)),
      name = "Average Air Temp (°C)",
      labels = number_format(accuracy = 0.1)
    )
  )

rh_state_plot <- plot_base +
  geom_line(
    data = avg_weath,
    aes(x = times, y = rescale(rh, to = c(0, 1))),
    color = "red", size = 1.2, inherit.aes = FALSE, linetype = 6
  ) +
  scale_y_continuous(
    name = "Prop of Locations",
    sec.axis = sec_axis(
      ~ rescale(., from = c(0, 1), to = range(avg_weath$rh, na.rm = TRUE)),
      name = "Average Rel Humidity (%)",
      labels = number_format(accuracy = 0.1)
    )
  )

# temp_state_plot <- temp_state_plot + theme(axis.title.x = element_blank(),
#                                            axis.text.x = element_blank(),
#                                            axis.ticks.x = element_blank())
# 
# rh_state_plot=rh_state_plot+ theme(axis.title.x = element_blank(),
#                                    axis.text.x = element_blank(),
#                                    axis.ticks.x = element_blank())

# Combine the plots
# combined_plot <- (temp_state_plot / rh_state_plot) + 
#   plot_layout(guides = 'collect') & 
#   theme(legend.position = "top")
# 
# # Display the plot
# combined_plot
# 
# png(width = 800, height = 600,filename="temp_rh_state_plot.png")
# combined_plot
# dev.off()

rainfall_state_plot <- plot_base +
  geom_line(
    data = avg_weath,
    aes(x = times, y = rescale(rainfall, to = c(0, 1))),
    color = "red", size = 1.2, inherit.aes = FALSE, linetype = 6
  ) +
  scale_y_continuous(
    name = "Prop of Locations",
    sec.axis = sec_axis(
      ~ rescale(., from = c(0, 1), to = range(avg_weath$rainfall, na.rm = TRUE)),
      name = "Average Rainfall (mm)",
      labels = number_format(accuracy = 0.1)
    )
  )
#rainfall_state_plot

wind_state_plot <- plot_base +
  geom_line(
    data = avg_weath,
    aes(x = times, y = rescale(wind_speed, to = c(0, 1))),
    color = "red", size = 1.2, inherit.aes = FALSE, linetype = 6
  ) +
  scale_y_continuous(
    name = "Prop of Locations",
    sec.axis = sec_axis(
      ~ rescale(., from = c(0, 1), to = range(avg_weath$wind_speed, na.rm = TRUE)),
      name = "Average Wind Speed (m/s)",
      labels = number_format(accuracy = 0.1)
    )
  )

rainfall_state_plot <- rainfall_state_plot + theme(axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank())

# Combine the plots
# combined_plot_2 <- (rainfall_state_plot / wind_state_plot) + 
#   plot_layout(guides = 'collect') & 
#   theme(legend.position = "top")
# 
# # Display the plot
# combined_plot_2
# 
# png(width = 800, height = 600,filename="rf_windspeed_state_plot.png")
# combined_plot_2
# dev.off()

library(patchwork)
temp_state_plot <- temp_state_plot + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

rh_state_plot <- rh_state_plot + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

rainfall_state_plot <- rainfall_state_plot + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# wind_state_plot mantiene l’asse x

# Combina i plot verticalmente
combined_plot <- (
  temp_state_plot /
    rh_state_plot /
    rainfall_state_plot /
    wind_state_plot
) + plot_layout(guides = "collect") & 
  theme(legend.position = "top")

# Visualizza il plot
combined_plot

png(width = 500, height = 900,filename="all_state_plot.png")
combined_plot
dev.off()

#
plot_base <- ggplot(S_summary_complete, aes(x = time, y = Proportion, fill = as.factor(ComfortLevel))) +
  geom_area(alpha = 0.6, size = 1, color = "black") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime",
                    labels=c("Cool","Neutral","Hot")) +
  labs(
    #title = "Temporal Evolution of Thermal Comfort Levels and Wind Speed",
    # x = "Date and Hour",
    x=" ",
    y = "Prop of Locations",
    fill = "Comfort Regime") +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%Y-%m-%d %H:%M") +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.8, "cm"),  # controls the box size
    legend.text = element_text(size = 15),  # controls the text size
    text = element_text(size = sz),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y.right = element_blank(),
    axis.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "top"
  )
utci_state_plot <- plot_base +
  geom_line(
    data = avg_weath,
    aes(x = times, y = rescale(UTCI, to = c(0, 1))),
    color = "red", size = 1.5, inherit.aes = FALSE, linetype = 6
  ) +
  scale_y_continuous(
    name = "Prop of Locations",
    sec.axis = sec_axis(
      ~ rescale(., from = c(0, 1), to = range(avg_weath$UTCI, na.rm = TRUE)),
      name = "UTCI (°C)",
      labels = number_format(accuracy = 0.1)
    )
  )



png(width = 800, height = 600,filename="utci_state_plot.png")
utci_state_plot
dev.off()


# hourly bar plot ---------------------------------------------------------

S_summary$time <- as.POSIXct(S_summary$time, format = "%Y-%m-%d %H:%M:%S")

# Extract hour from the 'time' column
S_summary$Hour <- hour(S_summary$time)

# Summarize the proportion of each comfort level at each hour
hourly_distribution <- S_summary %>%
  group_by(Hour, ComfortLevel) %>%
  summarise(Proportion = sum(Proportion)) %>%
  ungroup()

# Create the stacked bar plot
barplot_state=ggplot(hourly_distribution, aes(x = factor(Hour), 
                                              y = Proportion, 
                                              fill = as.factor(ComfortLevel))) +
  geom_bar(alpha=.6,stat = "identity", position = "fill", color = "black", size = .7) +  # Add solid black lines between bars
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime",
                    labels=c("Cool","Neutral","Hot")) +
  labs(
    #title = "Comfort Level Distribution by Hour",
    x = "Hour of the Day",
    y = "Proportion of Locations") +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.8, "cm"),  # controls the box size
    legend.text = element_text(size = 18),  # controls the text size
    text = element_text(size = sz),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y.right = element_blank(),
    axis.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "top"
  )
barplot_state

png(width = 800, height = 600,filename="barplot_state.png")
barplot_state
dev.off()

# station bar plot (better) -----------------------------------------------

# Convert data into long format for ggplot
S_long <- melt(S_est, id.vars = "time", variable.name = "Station", value.name = "ComfortLevel")

# Calculate proportions for each Comfort Level per station
proportion_data <- S_long %>%
  group_by(Station, ComfortLevel) %>%
  summarise(Proportion = n() / nrow(S_est)) %>%
  ungroup()

# Create stacked bar plot
stat_bar_plot=ggplot(proportion_data, aes(x = Station, y = Proportion, fill = as.factor(ComfortLevel))) +
  geom_bar(alpha=.6,stat = "identity" ,position = "fill", color = "black", size = 0.7) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime",
                    labels = c("Cool", "Neutral", "Hot")) +
  theme_minimal() +
  labs(
    #title = "Comfort Level Proportions for Each Station",
       x = "Station",
       y = "Proportion of Comfort Regimes") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



png(width = 800, height = 600,filename="stat_bar_plot.png")
stat_bar_plot
dev.off()

# Entropy -----------------------------------------------------------------

# Calculate the entropy for each station
calculate_entropy <- function(proportions) {
  -sum(proportions * log(proportions + 1e-10))  # Adding a small value to avoid log(0)
}

# Calculate proportions for each Comfort Level per station
proportion_data_entropy <- S_long %>%
  group_by(Station, ComfortLevel) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Station) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  summarise(Entropy = calculate_entropy(Proportion))

# Create bar plot of entropy for each station
x11()
ggplot(proportion_data_entropy, aes(x = Station, y = Entropy, fill = Station)) +
  geom_bar(alpha=.6,stat = "identity", color = "black") +
  scale_fill_manual(values = rep("pink",M)) +
  theme_minimal() +
  labs(title = "Entropy of Comfort Levels for Each Station",
       x = "Station",
       y = "Entropy") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # No legend needed as each bar is unique to a station


# heatmap -----------------------------------------------------------------

S_est_heat=S_est[,-1]
S_est_heat$t=timesY

S_long <- melt(S_est_heat, id.vars = "t", 
               variable.name = "Location", value.name = "ComfortLevel")

heatmap_stat_time <- ggplot(S_long, aes(x = t, y = Location, fill = as.factor(ComfortLevel))) +
  geom_tile(alpha=.6) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime") +
  scale_x_datetime(breaks = seq(as.POSIXct("2023-04-18 20:00:00"), 
                                as.POSIXct("2023-04-26 07:00:00"), 
                                by = "24 hours")
                   ,
                   labels = scales::date_format("%Y-%m-%d %H:%M")
                   ) +
  # Add vertical dotted lines at each break point
  geom_vline(xintercept = seq(as.POSIXct("2023-04-18 12:00:00"), 
                                         as.POSIXct("2023-04-26 07:00:00"), 
                                         by = "24 hours"), 
             color = "grey20", linetype = "dotted") +
  labs(
    x = "Time",
    y = "Station") +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
heatmap_stat_time

png(width = 800, height = 600,filename="heatmap_stat_time.png")
heatmap_stat_time
dev.off()

# Heatmap girata
sz=15
heatmap_time_stat <- ggplot(S_long, aes(x = Location, y = t, fill = as.factor(ComfortLevel))) +
  geom_tile(alpha = 0.6) +
  scale_fill_manual(
    values = c("lightblue", "lightgreen", "orange"), 
    name = "Comfort Regime",
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  scale_y_datetime(
    breaks = seq(as.POSIXct("2023-04-18 20:00:00"), 
                 as.POSIXct("2023-04-26 07:00:00"), 
                 by = "24 hours"),
    labels = scales::date_format("%Y-%m-%d %H:%M"),
    sec.axis = dup_axis(name = "Time")
  ) +
  # Add horizontal dotted lines at each break point
  geom_hline(yintercept = seq(as.POSIXct("2023-04-18 12:00:00"), 
                              as.POSIXct("2023-04-26 07:00:00"), 
                              by = "24 hours"), 
             color = "grey20", linetype = "dotted") +
  labs(
    x = "Station",
    y = NULL) +  # Remove y-axis label on the left
  theme_minimal() +
  theme(
    legend.key.size = unit(0.8, "cm"),  # controls the box size
    legend.text = element_text(size = 10),  # controls the text size
    text = element_text(size = sz),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y.right = element_blank(),
    axis.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "top"
  )


png(width = 800, height = 600,filename="heatmap_time_stat.png")
heatmap_time_stat
dev.off()


# leaflet with colored areas ----------------------------------------------

data_stat_number2=data_stat_number[,-1]
data_stat_number2$Station=paste("S",data_stat_number2$m)
Sys.setenv(TZ="GMT")
S_long <- S_est %>%
  pivot_longer(cols = starts_with("S"), names_to = "Station", values_to = "ComfortLevel")

# Merge the spatial data with comfort levels data
S_long <- S_long %>%
  left_join(data_stat_number2, by = "Station")

S_long$time=as.POSIXct(S_long$time,tz="GMT")+hours(8)

# Define a color palette for the comfort levels
comfort_colors <- list("1" = "blue", "2" = "green", "3" = "red")


# Shiny UI
library(htmlwidgets)
ui <- fluidPage(
  titlePanel("Dynamic Comfort Level Map"),
  
  sidebarLayout(
    sidebarPanel(
      # Slider bar to select a specific time point
      sliderInput("time", "Select Time:", 
                  min = min(S_long$time), 
                  max = max(S_long$time), 
                  value = min(S_long$time), 
                  timeFormat = "%Y-%m-%d %H:%M:%S", 
                  step = 3600, 
                  animate = animationOptions(interval = 1000, loop = TRUE)
                  ,
                  timezone = "GMT"
                  )
    ),
    mainPanel(
      # Display the selected time as text above the map
      # textOutput("selectedTime", container = tags$h3, align = "center"), 
      
      # Output map
      leafletOutput("comfortMap", height = 600)
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  
  # Reactive function to filter data based on the selected time
  filtered_data <- reactive({
    S_long %>% filter(time == input$time)
  })
  
  # Render the selected time as text
  output$selectedTime <- renderText({
    paste("Selected Time: ", format(input$time, "%Y-%m-%d %H:%M:%S"
                                    ,tz = "GMT"
                                    )
          )
  })
  
  # Create a leaflet map that updates based on the selected time
  output$comfortMap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%  # Add default tiles
      setView(lng = mean(S_long$longitude), lat = mean(S_long$latitude), zoom = 12)
  })
  
  # Observe the selected time and update the map accordingly
  observe({
    leafletProxy("comfortMap", data = filtered_data()) %>%
      clearMarkers() %>%
      clearControls() %>%
      addCircleMarkers(
        lng = ~longitude, lat = ~latitude,
        color = ~unname(comfort_colors[as.character(ComfortLevel)]),
        fillColor = ~unname(comfort_colors[as.character(ComfortLevel)]),
        fillOpacity = 1,
        radius = 10,
        stroke = TRUE,
        weight = 1,
        opacity = 10,  
        label = ~paste0("Station: ", Station, "<br>Comfort Level: ", ComfortLevel),
        labelOptions = labelOptions(style = list("color" = "black"), direction = "auto")
      ) %>%
      addLegend(
        position = "topright",
        colors = unlist(comfort_colors),
        labels = c("Cool", "Neutral", "Hot"),
        title = "Comfort Levels",
        opacity = 0.7
      )
  })
}

# Run the Shiny App
app <- shinyApp(ui, server)

app

# Load svf and gvf --------------------------------------------------------

# Load geo data
load("geo_feat.Rdata")

D_geo_weath=distm(loc_weath[,c("longitude","latitude")], 
                  geo_feat[,c("longitude","latitude")], 
                  fun = distGeo)/1000
indx=apply(D_geo_weath, 1,which.min)

# 21 is teh closest with no 0 entries for the indexes
indx8=sort(D_geo_weath[8,])[21]
i8=which(D_geo_weath[8,]==indx8)
geo_feat[i8,]

indx[8]=i8

# Maximum distance between locations of SVF and GVF and weather station is less than 800 m
max(diag(D_geo_weath[,indx]))

data_geo_weath=geo_feat[indx,]
data_geo_weath$station=loc_weath$id

# Does the order of the stations match? YES!
plot(data_geo_weath$longitude,data_geo_weath$latitude, xlab = "Long", 
     ylab = "Lat", main = " ", pch = 19, col = as.numeric(data_geo_weath$station))
points(loc_weath$longitude,loc_weath$latitude, col = as.numeric(loc_weath$id), pch = 2)


data_SVF_GVF=merge(data_geo_weath,data_stat_number[,1:2], by= "station")
data_SVF_GVF=select(data_SVF_GVF, -c(station))
rownames(data_SVF_GVF)=data_geo_weath$m
data_SVF_GVF[order(data_SVF_GVF$green_view_mean),]
data_SVF_GVF[order(data_SVF_GVF$sky_view_mean),]

# Box plot -------------------------------------------------------------

library(ggplot2)

# Assuming df is your dataframe containing the data
bp_air=ggplot(Y_res, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), y = air_temp, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels=c("Cool","Neutral","Hot")) +
  labs(x = " ", y = "Air Temperature (°C)", fill = "Comfort Regime", color = "State") +
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey40") +  # Add grey dotted lines between stations
  theme_minimal() +
  theme(text = element_text(size = 13),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

bp_rh=ggplot(Y_res, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), 
                  y = rh, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels=c("Cool","Neutral","Hot")) +
  labs(x = " ", y = "Rel Humidity (%)", fill = "Comfort Regime", color = "State") +
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey40") +  # Add grey dotted lines between stations
  theme_minimal() +
  theme(text = element_text(size = 13),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

bp_rf=ggplot(Y_res, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), 
                  y = rainfall, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels=c("Cool","Neutral","Hot")) +
  labs(x = " ", y = "Rainfall (mm)", fill = "Comfort Regime", color = "State") +
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey40") +  # Add grey dotted lines between stations
  theme_minimal() +
  theme(text = element_text(size = 13),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

bp_ws=ggplot(Y_res, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), 
                  y = wind_speed, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels=c("Cool","Neutral","Hot")) +
  labs(x = "Station", y = "Wind Speed (m/s)", fill = "Comfort Regime", color = "State") +
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey40") +  # Add grey dotted lines between stations
  theme_minimal() +
  theme(text = element_text(size = 13),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

library(gridExtra)
library(ggpubr)
png(width = 900, height = 500,filename="bp_air_rh.png")
ggarrange(bp_air, bp_rh, ncol = 1, nrow = 2,common.legend = T)
dev.off()

png(width = 900, height = 500,filename="bp_rf_ws.png")
ggarrange(bp_rf, bp_ws, ncol = 1, nrow = 2,common.legend = T)
dev.off()

png(width = 500, height = 900,filename="bp_all.png")
ggarrange(bp_air, bp_rh,bp_rf, bp_ws, ncol = 1, nrow = 4,common.legend = T)
dev.off()

ggarrange(bp_rf, bp_ws, ncol = 1, nrow = 2,common.legend = T)

head(data_SVF_GVF)
Y_res2=Y_res%>%
  left_join(data_SVF_GVF, by = "m")%>%
  mutate(SVF=round(sky_view_mean,2),GVF=round(green_view_mean,2))

# Create the boxplot
# Create the boxplot
ggplot(Y_res2, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), y = air_temp, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels = c("1" = "Cool", "2" = "Neutral", "3" = "Hot")) +
  
  # Add SVF and GVF labels using annotate()
  annotate("text", x = 1:length(unique(Y_res2$m)), y = max(Y_res2$air_temp,na.rm=T) + 1, 
           label = paste0("SVF: ", round(data_SVF_GVF$sky_view_mean, 2), "\nGVF: ", round(data_SVF_GVF$green_view_mean, 2)),
           size = 3, fontface = "bold", color = "black") +
  
  labs(x = "Station", y = "Air Temperature (°C)", fill = "State") +
  
  # Add dotted vertical lines between each station group
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey30") +  
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off")  # Prevents the labels from being clipped

ggplot(Y_res2, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), y = rh, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels = c("1" = "Cool", "2" = "Neutral", "3" = "Hot")) +
  
  # Add SVF and GVF labels using annotate()
  annotate("text", x = 1:length(unique(Y_res2$m)), y = max(Y_res2$rh,na.rm=T) + 5, 
           label = paste0("SVF: ", round(data_SVF_GVF$sky_view_mean, 2), "\nGVF: ", round(data_SVF_GVF$green_view_mean, 2)),
           size = 3, fontface = "bold", color = "black") +
  
  labs(x = "Station", y = "Relative Humidity (%)", fill = "State") +
  
  # Add dotted vertical lines between each station group
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey30") +  
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off")  # Prevents the labels from being clipped

ggplot(Y_res2, aes(x = factor(m, labels = paste0("S", 1:length(unique(m)))), y = wind_speed, fill = factor(State))) +
  geom_boxplot(show.legend = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("1" = "lightblue", "2" = "lightgreen", "3" = "orange"),
                    labels = c("1" = "Cool", "2" = "Neutral", "3" = "Hot")) +
  
  # Add SVF and GVF labels using annotate()
  annotate("text", x = 1:length(unique(Y_res2$m)), y = max(Y_res2$wind_speed,na.rm=T) + 5, 
           label = paste0("SVF: ", round(data_SVF_GVF$sky_view_mean, 2), "\nGVF: ", 
                          round(data_SVF_GVF$green_view_mean, 2)),
           size = 3, fontface = "bold", color = "black") +
  
  labs(x = "Station", y = "Wind Speed (m/s)", fill = "State") +
  
  # Add dotted vertical lines between each station group
  geom_vline(xintercept = seq(1.5, length(unique(Y_res$m)) - 0.5, by = 1), 
             linetype = "dotted", color = "grey30") +  
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off")  # Prevents the labels from being clipped


S_long <- melt(S_est, id.vars = "time", variable.name = "Station", value.name = "ComfortLevel")

# Calculate proportions for each Comfort Level per station
proportion_data <- S_long %>%
  group_by(Station, ComfortLevel) %>%
  summarise(Proportion = n() / nrow(S_est)) %>%
  ungroup()

# Create stacked bar plot
stat_bar_plot_svf_gvf=ggplot(proportion_data, aes(x = Station, y = Proportion, fill = as.factor(ComfortLevel))) +
  geom_bar(alpha=.6,stat = "identity" ,position = "fill", color = "black", size = 0.7) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Regime",
                    labels = c("Cool", "Neutral", "Hot")) +
  annotate("text", x = 1:length(unique(Y_res2$m)), y = 1.1, 
           label = paste0("SVF: ", round(data_SVF_GVF$sky_view_mean, 2), "\nGVF: ", 
                          round(data_SVF_GVF$green_view_mean, 2)),
           size = 3, fontface = "bold", color = "black") +
  theme_minimal() +
  labs(
    #title = "Comfort Level Proportions for Each Station",
    x = "Station",
    y = "Proportion of Comfort Regimes") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



png(width = 800, height = 600,filename="stat_bar_plot_svf_gvf.png")
stat_bar_plot_svf_gvf
dev.off()