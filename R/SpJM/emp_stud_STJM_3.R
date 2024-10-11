library(dplyr)
library(lubridate) 
library(tidyr)
library(tidyverse)
library(shiny)
library(geosphere)
library(caret)
library(leaflet)
library(reshape2)
# Set Singapore time zone
Sys.setenv(TZ="UTC")
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

load("D:/git/dynamicComf/dynamicComf/R/SpJM/Y_3.Rdata")

# Merge with pre-existing data
Y_4=Y_complete[,c("t","m","time")]%>%
  left_join(Y_3,by=c("t","m"))

Y_5=Y_4%>%
  left_join(data_utci,by="time")

Y_6=Y_5[-(1:56),]
Y_6=subset(Y_6,select=-time)
###

# leaflet map singapore ---------------------------------------------------

leaflet(data_stat_number) %>%
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
  setView(lng = mean(data_stat_number$longitude), lat = mean(data_stat_number$latitude), zoom = 12)


# STJM fit ----------------------------------------------------------------



source("Utils.R")
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

lambda=.05
gamma=.05

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

tapply(Y_res$air_temp,Y_res$State,mean,na.rm=T)
tapply(Y_res$rh,Y_res$State,mean,na.rm=T)
tapply(Y_res$rainfall,Y_res$State,mean,na.rm=T)
tapply(Y_res$wind_speed,Y_res$State,mean,na.rm=T)
#tapply(Y_res$wind_dir,Y_res$State,mean,na.rm=T)
#tapply(Y_res$rainy,Y_res$State,Mode)
tapply(Y_res$windy,Y_res$State,Mode)
tapply(Y_res$hour,Y_res$State,Mode)
tapply(Y_res$UTCI,Y_res$State,mean,na.rm=T)

TY=unique(Y_6$t)
timesY=unique(times)
TT_6=length(TY)

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
  geom_area(alpha = 0.6, size = .7, color = "black") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Level") +
  labs(
    #title = "Temporal Evolution of Thermal Comfort Levels and Wind Speed",
    x = "Date and Hour",
    y = "Proportion of Locations",
    fill = "Comfort Level") +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%Y-%m-%d %H:%M") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
# +
#   geom_vline(xintercept = x_breaks, color = "grey30", linetype = "dotted", size = 0.5)

wind_state_plot <- plot_base +
  geom_line(data = avg_weath, aes(x = times, y = wind_speed 
                                  / 10
  ),  # Scale wind speed for better visualization
  color = "red", size = 1, inherit.aes = FALSE,linetype=6) +
  #geom_point(data = avg_weath, aes(x = times, y = wind_speed / 10),  # Scale wind speed for better visualization
  # color = "blue", size = 2, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Proportion of Locations",
    sec.axis = sec_axis(~ . * 10, name = "Average Wind Speed (m/s)", labels = scales::number_format(accuracy = 0.1))
  )

wind_state_plot

rainfall_state_plot <- plot_base +
  geom_line(data = avg_weath, aes(x = times, y = rainfall 
                                  #/ 10
  ),  
  color = "red", size = 1, inherit.aes = FALSE,linetype=6) +
  #geom_point(data = avg_weath, aes(x = times, y = wind_speed / 10),  # Scale wind speed for better visualization
  # color = "blue", size = 2, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Proportion of Locations",
    sec.axis = sec_axis(~ . * 1, name = "Average Rainfall (mm)", 
                        labels = scales::number_format(accuracy = 0.1))
  )
rainfall_state_plot

temp_state_plot <- plot_base +
  geom_line(data = avg_weath, aes(x = times, y = air_temp 
                                  / 40
  ),  
  color = "red", size = 1, inherit.aes = FALSE,linetype=6) +
  #geom_point(data = avg_weath, aes(x = times, y = wind_speed / 10),  # Scale wind speed for better visualization
  # color = "blue", size = 2, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Proportion of Locations",
    sec.axis = sec_axis(~ . * 40, name = "Average Air Temp (°C)", 
                        labels = scales::number_format(accuracy = 0.1))
  )
temp_state_plot

rh_state_plot <- plot_base +
  geom_line(data = avg_weath, aes(x = times, y = rh 
                                  / 100
  ),  #
  color = "red", size = 1, inherit.aes = FALSE,linetype=6) +
  #geom_point(data = avg_weath, aes(x = times, y = wind_speed / 10),  # Scale wind speed for better visualization
  # color = "blue", size = 2, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Proportion of Locations",
    sec.axis = sec_axis(~ . * 100, name = "Relative Humidity (%)", 
                        labels = scales::number_format(accuracy = 0.1))
  )
rh_state_plot

utci_state_plot <- plot_base +
  geom_line(data = avg_weath, aes(x = times, y = UTCI 
                                  / 40
  ),  #
  color = "red", size = 1, inherit.aes = FALSE,linetype=6) +
  #geom_point(data = avg_weath, aes(x = times, y = wind_speed / 10),  # Scale wind speed for better visualization
  # color = "blue", size = 2, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Proportion of Locations",
    sec.axis = sec_axis(~ . * 40, name = "UTCI (°C)", 
                        labels = scales::number_format(accuracy = 0.1))
  )
utci_state_plot


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
  geom_bar(stat = "identity", position = "fill", color = "black", size = 0.5) +  # Add solid black lines between bars
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Level") +
  labs(
    #title = "Comfort Level Distribution by Hour",
    x = "Hour of Day",
    y = "Proportion of Locations") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top")
barplot_state


# station bar plot (better) -----------------------------------------------

# Convert data into long format for ggplot
S_long <- melt(S_est, id.vars = "time", variable.name = "Station", value.name = "ComfortLevel")

# Calculate proportions for each Comfort Level per station
proportion_data <- S_long %>%
  group_by(Station, ComfortLevel) %>%
  summarise(Proportion = n() / nrow(S_est)) %>%
  ungroup()

# Create stacked bar plot
ggplot(proportion_data, aes(x = Station, y = Proportion, fill = as.factor(ComfortLevel))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Level",
                    labels = c("Cool", "Neutral", "Hot")) +
  theme_minimal() +
  labs(title = "Comfort Level Proportions for Each Station",
       x = "Station",
       y = "Proportion of Comfort Levels") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


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
ggplot(proportion_data_entropy, aes(x = Station, y = Entropy, fill = Station)) +
  geom_bar(stat = "identity", color = "black") +
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

ggplot(S_long, aes(x = Location, y = t, fill = as.factor(ComfortLevel))) +
  geom_tile() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange"), 
                    name = "Comfort Level"
                    # ,
                    # labels = c("Hot", "Neutral", "Cool")
  ) +
  labs(
    #title = "Heatmap",
    x = "Spatial Location",
    y = "Time (Hourly)") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



# leaflet with colored areas ----------------------------------------------

data_stat_number2=data_stat_number[,-1]
data_stat_number2$Station=paste("S",data_stat_number2$m)

S_long <- S_est %>%
  pivot_longer(cols = starts_with("S"), names_to = "Station", values_to = "ComfortLevel")

# Merge the spatial data with comfort levels data
S_long <- S_long %>%
  left_join(data_stat_number2, by = "Station")

# Define a color palette for the comfort levels
comfort_colors <- list("1" = "lightblue", "2" = "lightgreen", "3" = "orange")

# Shiny UI
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
                  animate = animationOptions(interval = 1000, loop = TRUE))
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
    paste("Selected Time: ", format(input$time, "%Y-%m-%d %H:%M:%S"))
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
        labels = c("Cool (1)", "Neutral (2)", "Hot (3)"),
        title = "Comfort Levels",
        opacity = 0.7
      )
  })
}

# Run the Shiny App
shinyApp(ui, server)
