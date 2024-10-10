library(dplyr)
library(lubridate) 
library(tidyr)
library(tidyverse)
library(shiny)
library(geosphere)
library(caret)
library(leaflet)

# Set Singapore time zone
Sys.setenv(TZ="GMT+8")

source("Utils.R")

load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/cleaned data/cleaned_data.RData")

# Dataset pulito da Antonio
load("C:/Users/federico/OneDrive - CNR/General - Comfort/Kaggle/singapore/q_outdoor.RData")

cozie_q.wd <- cozie_q.wd %>%
  mutate(date = as.Date(time))  # Extract the date (year-month-day) from the time column
wdn="1 hour"
cozie_q.wd$time=round_date(cozie_q.wd$time,wdn)
cozie_q.wd$ws_timestamp_start=round_date(cozie_q.wd$time,wdn)
cozie_q.wd$ws_timestamp_location=round_date(cozie_q.wd$ws_timestamp_location,wdn)

# Retrieve data for future comparison with STJM estimates
time_range=c(as.Date('2023-04-18'),as.Date('2023-04-25'))
cozie_compare=cozie_q.wd %>%
  group_by(date) %>%
  filter(date >= time_range[1] & date <= time_range[2])

range(cozie_compare$time)
table(cozie_compare$q_thermal_preference)
# # Count the number of observations for each specific day
# daily_counts <- cozie_q.wd %>%
#   group_by(date) %>%
#   summarise(observation_count = n())
# 
# # Print the result
# print(daily_counts)
# plot(daily_counts$date, daily_counts$observation_count, type = "l", xlab = "Date", ylab = "Number of observations", main = "Number of observations per day")
# 
# max_date <- daily_counts$date[
# which.max(daily_counts$observation_count)]  # Find the date with the most observations
# 
# # Select only the data for the date with the most observations
# max_date_data <- cozie_q.wd %>%
#   filter(date == max_date)

stats=intersect(colnames(air_temp),colnames(rainfall))
stats=intersect(stats,colnames(RH))
stats=intersect(stats,colnames(wind_speed))
stats=intersect(stats,colnames(wind_dir))

data_air_temp <- air_temp%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
apply(data_air_temp, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop S102
data_air_temp=subset(data_air_temp, select = -c(S102))

# Summarize by hourly mean
wdn="1 hour"
data_air_temp=data_air_temp%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

data_air_temp=data.frame(data_air_temp)
max(diff(data_air_temp$time))
range(data_air_temp$time)
# For RH, retrieve the data for the same time range as the comfort data

data_rh <- RH%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
apply(data_rh, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop S102
data_rh=subset(data_rh, select = -c(S102))

# Summarize by hourly mean
wdn="1 hour"
data_rh=data_rh%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_rh$time))

# For rainfall, retrieve the data for the same time range as the comfort data

data_rainfall <- rainfall%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
sts=apply(data_rainfall, 2, function(x) 100*sum(is.na(x))/length(x))

# Drop stations with 100% NAs
drp=which(sts==100)
data_rainfall=subset(data_rainfall, select = -drp)

# Summarize by hourly mean
wdn="1 hour"
data_rainfall=data_rainfall%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_rainfall$time))

# For wind_speed, retrieve the data for the same time range as the comfort data

data_wind_speed <- wind_speed%>%
  mutate(date = as.Date(time)) %>%
  filter(date >= time_range[1] & date <= time_range[2])%>%
  dplyr::select(-date)

# Count NAs for each station
sts=apply(data_wind_speed, 2, function(x) 100*sum(is.na(x))/length(x))
sts

# Drop stations with 100% NAs
drp=which(sts==100)
data_wind_speed=subset(data_wind_speed, select = -drp)

# Summarize by hourly mean
wdn="1 hour"
data_wind_speed=data_wind_speed%>%
  group_by(time=floor_date(time,wdn))%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

max(diff(data_wind_speed$time))

# For wind_dir, retrieve the data for the same time range as the comfort data

# data_wind_dir <- wind_dir%>%
#   mutate(date = as.Date(time)) %>%
#   filter(date >= time_range[1] & date <= time_range[2])%>%
#   dplyr::select(-date)
# 
# # Count NAs for each station
# sts=apply(data_wind_dir, 2, function(x) 100*sum(is.na(x))/length(x))
# sts
# 
# # Drop stations with 100% NAs
# drp=which(sts==100)
# data_wind_dir=subset(data_wind_dir, select = -drp)
# 
# # Summarize by hourly mean
# wdn="1 hour"
# data_wind_dir=data_wind_dir%>%
#   group_by(time=floor_date(time,wdn))%>%
#   summarise_if(is.numeric, mean, na.rm = TRUE)
# 
# max(diff(data_wind_dir$time))

TT=nrow(data_air_temp)
data_air_temp$t=1:TT


# Introduce temporal gaps
pGap=.05
time_set=data_air_temp$time
time_set=time_set[-which(time_set%in%unique(cozie_compare$time))]
set.seed(1234)
gaps=sort(sample(time_set,round(TT*pGap),replace=F))
gaps=which(data_air_temp$time%in%gaps)

times_complete=data_air_temp$time

data_air_temp=data.frame(data_air_temp)
data_air_temp$t=1:TT
data_air_temp=data_air_temp[-gaps,]

data_rh=data.frame(data_rh)
data_rh$t=1:TT
data_rh=data_rh[-gaps,]

data_rainfall=data.frame(data_rainfall)
data_rainfall$t=1:TT
data_rainfall=data_rainfall[-gaps,]

data_wind_speed=data.frame(data_wind_speed)
data_wind_speed$t=1:TT
data_wind_speed=data_wind_speed[-gaps,]

# data_wind_dir=data.frame(data_wind_dir)
# data_wind_dir$t=1:TT
# data_wind_dir=data_wind_dir[-gaps,]


# Keep only relevant stations
stats=intersect(colnames(data_air_temp),colnames(data_rainfall))
stats=intersect(stats,colnames(data_rh))
stats=intersect(stats,colnames(data_wind_speed))
#stats=intersect(stats,colnames(data_wind_dir))

data_air_temp=data_air_temp[,stats]
data_rainfall=data_rainfall[,stats]
data_rh=data_rh[,stats]
data_wind_speed=data_wind_speed[,stats]
#data_wind_dir=data_wind_dir[,stats]

# Long format
data_air_temp_long <- data_air_temp %>%
  gather(key = "station", value = "air_temp", -c(t,time))
data_rh_long <- data_rh %>%
  gather(key = "station", value = "rh", -c(t,time))
data_rainfall_long <- data_rainfall %>%
  gather(key = "station", value = "rainfall", -c(t,time))
data_wind_speed_long <- data_wind_speed %>%
  gather(key = "station", value = "wind_speed", -c(t,time))
# data_wind_dir_long <- data_wind_dir %>%
#   gather(key = "station", value = "wind_dir", -c(t,time))

# Merge all the weather data
data_weather_long <- data_air_temp_long %>%
  left_join(data_rh_long, by = c("t","time", "station")) %>%
  left_join(data_rainfall_long, by = c("t","time", "station")) %>%
  left_join(data_wind_speed_long, by = c("t","time", "station")) 

# Retrieve relevant locations
loc_weath=locations %>%
  filter(id %in% stats)


# Save the data
save(cozie_compare,loc_weath,
     data_air_temp, data_rh, data_rainfall, data_wind_speed, 
     data_weather_long, 
     file = "data_weather_v2.RData")

# HERE --------------------------------------------------------------------


# Load weather data
load("data_weather_v2.RData")

M=dim(loc_weath)[1]
data_stat_number=data.frame(station=loc_weath$id, m=1:M,
                            longitude=loc_weath$longitude,
                            latitude=loc_weath$latitude)
data_stat_number=data_stat_number[order(data_stat_number$m),]

# data_stat_number=merge(data_stat_number,data_geo_weath[,c("station",                                                           "green_view_mean",                                                           "sky_view_mean",
#                                                           "building_view_mean")],
#                        by="station")
# 
# Y=data_weather_long %>%
#   left_join(data_stat_number[,c("station","m","green_view_mean", "sky_view_mean")], by = c("station"))

Y=data_weather_long %>%
  left_join(data_stat_number[,c("station","m")], by = c("station"))



# Order based on t and m
Y=Y[order(Y$t,Y$m),]
rownames(Y)=NULL

# Add categorical features
#Y$rainy=ifelse(Y$rainfall>0,1,0)+1
#Y$windy=ifelse(Y$wind_speed>3,1,0)+1

# Construct rainy feature according to Table 1 of Marsico 2021
# Y$rainy=factor(
#   cut(Y$rainfall, breaks = c(-Inf, 2.5, 10, 50, Inf), labels = c(1, 2, 3, 4)),
#   levels = c(1, 2, 3, 4)
# )
# Y$rainy=droplevels(Y$rainy)
# 
# # Same for Beaufort scale, remember to start from light air
Y$windy=factor(
  cut(Y$wind_speed, breaks = c(-Inf, 1.5,
                               3.3,5.4,7.9,10.7,Inf), labels = 1:6,
      levels = 1:6
  ))
Y$windy=droplevels(Y$windy)

# NOTA BENE: remember to say that only one observation falls into the 7th category so we
# include it into the 6th one

Y$hour=hour(Y$time)

# Change NaN in NA
for(i in 2:ncol(Y)){
  Y[,i]=ifelse(is.nan(Y[,i]),NA,Y[,i])
}

Y$hour=as.factor(Y$hour)
#Y$rainy=as.factor(Y$rainy)
Y$windy=as.factor(Y$windy)

# Drop time and station columns
# times=Y[,1]
# Y=Y[,-c(1,3)]
# head(Y)
# str(Y)

summary(Y)
cor(Y[complete.cases(Y),4:7])

Y_complete=Y
Y=subset(Y,select=-c(time,station))

wdn=5
library(zoo)
Y_rollmean=Y[,2:6]%>%group_by(m)%>%
  mutate_if(is.numeric, funs(rollapply(.,width=wdn,FUN=mean,align="right",fill=NA)))
colnames(Y_rollmean)[1:4]=paste0(colnames(Y_rollmean)[1:4],"_rollmean")
Y_rollmean$t=Y$t

# Y_rollmean%>%group_by(m)%>%
#   ggplot(aes(x=t,y=air_temp_rollmean,color=as.factor(m)))+
#   geom_line()

Y_rollsd=Y[,2:6]%>%group_by(m)%>%
  mutate_if(is.numeric, funs(rollapply(.,width=wdn,FUN=sd,align="right",fill=NA)))
colnames(Y_rollsd)[1:4]=paste0(colnames(Y_rollsd)[1:4],"_rollsd")
Y_rollsd$t=Y$t

# Y_rollsd%>%group_by(m)%>%
#   ggplot(aes(x=t,y=air_temp_rollsd,color=as.factor(m)))+
#   geom_line()

# NO CORRELATION WITH RAINFALL SINCE THERE ARE A LOT OF ZEROS
# Y_corr <- Y %>%
#   group_by(m) %>%
#   mutate(
#     airtemp_rh_corr = rollapply(
#       data = air_temp, 
#       width = wdn, 
#       FUN = function(x) cor(x, rh[seq_along(x)]), 
#       fill = NA, 
#       align = "right"
#     ),
#     airtemp_windspeed_corr = rollapply(
#       data = air_temp, 
#       width = wdn, 
#       FUN = function(x) cor(x, wind_speed[seq_along(x)]), 
#       fill = NA, 
#       align = "right"
#     ),
#     rh_windspeed_corr = rollapply(
#       data = rh, 
#       width = wdn, 
#       FUN = function(x) cor(x, wind_speed[seq_along(x)]), 
#       fill = NA, 
#       align = "right"
#     )
#   )
# Y_corr_2=Y_corr[,c(1,6,9:11)]

Y_2=merge(Y,Y_rollmean,by=c("m","t"))
Y_2=merge(Y_2,Y_rollsd,by=c("m","t"))
#Y_2=merge(Y_2,Y_corr_2,by=c("m","t"))

Y_2=Y_2[order(Y_2$t,Y_2$m),]
rownames(Y_2)=NULL

Y_3=Y_2[-(1:((wdn-1)*M)),]
TT_3=length(unique(Y_3$t))
apply(Y_3,2,function(x) sum(is.na(x)))

# Graphical analysis

Y_3%>%group_by(m)%>%
  ggplot(aes(x=t,y=rainfall_rollsd,color=as.factor(m)))+
  geom_line()

Y_3%>%group_by(m)%>%
  ggplot(aes(x=t,y=rainfall_rollmean,color=as.factor(m)))+
  geom_line()


# Summ stat ---------------------------------------------------------------

summary(Y_3)
table(diff(timesY))/length(timesY)*100
cor(Y_3[complete.cases(Y_3),3:6])
library(ppcor)
pcor(Y_3[complete.cases(Y_3), 3:6])

# Leaflet map of Singapore ------------------------------------------------
# Create a leaflet map of Singapore
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


# Y_3%>%group_by(m)%>%
#   ggplot(aes(x=t,y=airtemp_rh_corr,color=as.factor(m)))+
#   geom_line()
# 
# Y_3%>%group_by(m)%>%
#   ggplot(aes(x=t,y=airtemp_windspeed_corr,color=as.factor(m)))+
#   geom_line()

source("Utils.R")
D=distm(data_stat_number[,c("longitude","latitude")], 
        data_stat_number[,c("longitude","latitude")], 
        fun = distGeo)/1000

lambda=.05
gamma=.05

fit=STjumpDist(Y_3,3,D,
               jump_penalty=lambda,
               spatial_penalty=gamma,
               initial_states=NULL,
               max_iter=10, n_init=10, 
               tol=1e-4, 
               verbose=T,timeflag=T)

# State characterization

State=c(t(fit$best_s))
State=order_states_condMean(Y_3$air_temp,State)

S_est=matrix(State,ncol=M,byrow = T)

times=Y_complete$time[-(1:((wdn-1)*M))]
Y_res=data.frame(Y_3,State,times)
table(Y_res$State)/length(times)*100

tapply(Y_res$air_temp,Y_res$State,mean,na.rm=T)
tapply(Y_res$rh,Y_res$State,mean,na.rm=T)
tapply(Y_res$rainfall,Y_res$State,mean,na.rm=T)
tapply(Y_res$wind_speed,Y_res$State,mean,na.rm=T)
#tapply(Y_res$wind_dir,Y_res$State,mean,na.rm=T)
#tapply(Y_res$rainy,Y_res$State,Mode)
tapply(Y_res$windy,Y_res$State,Mode)
tapply(Y_res$hour,Y_res$State,Mode)
# tapply(Y_res$green_view_mean,Y_res$State,mean,na.rm=T)
# tapply(Y_res$sky_view_mean,Y_res$State,mean,na.rm=T)
# tapply(Y_res$ALLSKY_SFC_SW_DWN,Y_res$State,mean,na.rm=T)
# tapply(Y_res$ALLSKY_SFC_LW_DWN,Y_res$State,mean,na.rm=T)
# tapply(Y_res$airtemp_rh_corr,Y_res$State,mean,na.rm=T)
# tapply(Y_res$airtemp_windspeed_corr,Y_res$State,mean,na.rm=T)
# tapply(Y_res$rh_windspeed_corr,Y_res$State,mean,na.rm=T)


TY=unique(Y_3$t)
timesY=unique(times)
TT_3=length(TY)

# Strano che ws_timestamp_start ha SOLO orari dalle 2 di notte alle 11 di mattina, forse meglio prendere time
q_data=cozie_compare[,c("q_thermal_preference","time","ws_longitude","ws_latitude")]
#q_data=cozie_compare[,c("q_thermal_preference","ws_timestamp_start","ws_longitude","ws_latitude")]
colnames(q_data)[2:4]=c("time","longitude","latitude")
# wdn="1 hour"
# q_data$time=round_date(q_data$time,wdn)
q_data$State=q_data$q_thermal_preference
q_data$State=recode(q_data$State, "Warmer" = 1, "No change" = 3, "Cooler" = 2)

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


# Shiny plot of the map ---------------------------------------------------
# Hourly average of the data
data_hour_av=Y_res%>%
  group_by(times)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

S_est_matrix=as.matrix(S_est[,-1])
# Define UI
ui <- fluidPage(
  titlePanel("Dynamic Plot for Varying t"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("t", "Select t:", 
                  min = min(TY), 
                  max = max(TY), 
                  value = min(TY), 
                  step = NULL)  # step is NULL to allow discrete values
    ),
    
    mainPanel(
      plotOutput("dynamicPlot")
    )
  )
)

server <- function(input, output, session) {
  # Update slider to only allow values from TY
  observe({
    updateSliderInput(session, "t", 
                      min = 1, 
                      max = length(TY), 
                      value = TY[1], 
                      step = 1)
  })
  
  output$dynamicPlot <- renderPlot({
    # Ensure t is one of the values from TY
    t_value <- input$t
    
    # Define color mapping based on selected t_value
    #par(mfrow=c(1,2))
    # plot(result$spatial_points, 
    #      col = S_true[t_value, ], 
    #      pch = 19, cex = 1.5, 
    #      main = paste("True"))
    plot(data_stat_number[,c("longitude","latitude")], 
         col = S_est_matrix[t_value, ], 
         pch = 17, cex = 1.5, 
         #main = bquote("ST-JM " ~ lambda == .(lambda) ~ " and " ~ gamma == .(gamma))
         main = paste("t = ", timesY[t_value])
    )
    if(any(q_data$time==timesY[t_value])){
      points(q_data$longitude[which(q_data$time==timesY[t_value])], 
             q_data$latitude[which(q_data$time==timesY[t_value])], 
             col = q_data$State[which(q_data$time==timesY[t_value])], pch = 19)
    }
    #legend("topright", legend = 1:3, fill = 1:3)
    legend("topright", legend = c("Cool","Neutral","Hot"), fill = c(1,3,2))
  })
}

# Run the app
shinyApp(ui = ui, server = server)


# Time evol of weather vars and state seq decod ---------------------------

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
  scale_fill_manual(values = c("lightblue", "orange", "lightgreen"), 
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
        legend.position = "top")+
  geom_vline(xintercept = x_breaks, color = "grey30", linetype = "dotted", size = 0.5)

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

# Hourly bar plot ----------------------------------------------------

# Ensure 'time' column is in POSIXct format
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
  scale_fill_manual(values = c("lightblue", "orange", "lightgreen"), 
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


# Station bar plot --------------------------------------------------------

# Step 1: Convert the data frame to long format
S_est_long <- S_est %>%
  pivot_longer(cols = starts_with("S "), names_to = "Station", values_to = "ComfortLevel")

# Step 2: Summarize to get the count of each comfort level for each station
station_summary <- S_est_long %>%
  group_by(Station, ComfortLevel) %>%
  summarise(Count = n()/length(timesY), .groups = 'drop')

station_summary$Station <- factor(station_summary$Station, levels = paste0("S ", 1:14))

# Step 4: Create a bar plot showing the distribution of comfort levels for each station
ggplot(station_summary, aes(x = Station, y = Count, fill = as.factor(ComfortLevel))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("lightblue", "orange", "lightgreen"), name = "Comfort Level") +
  labs(
    #title = "Distribution of Comfort Levels for Each Station",
       x = "Station",
       y = "% Comfort Levels") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))



# Heatmap -----------------------------------------------------------------

S_est_heat=S_est[,-1]
S_est_heat$t=timesY

S_long <- melt(S_est_heat, id.vars = "t", 
               variable.name = "Location", value.name = "ComfortLevel")

ggplot(S_long, aes(x = Location, y = t, fill = as.factor(ComfortLevel))) +
  geom_tile() +
  scale_fill_manual(values = c("lightblue", "orange", "lightgreen"), 
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
comfort_colors <- list("1" = "lightblue", "2" = "orange", "3" = "lightgreen")

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
        radius = 8,
        stroke = TRUE,
        weight = 1,
        opacity = 1,  
        label = ~paste0("Station: ", Station, "<br>Comfort Level: ", ComfortLevel),
        labelOptions = labelOptions(style = list("color" = "black"), direction = "auto")
      ) %>%
      addLegend(
        position = "topright",
        colors = unlist(comfort_colors),
        labels = c("Cool (1)", "Hot (2)", "Neutral (3)"),
        title = "Comfort Levels",
        opacity = 0.7
      )
  })
}

# Run the Shiny App
shinyApp(ui, server)

# leaflet map with most prob comfort level --------------------------------

# Most_prob_lev=station_summary%>%group_by(Station)%>%
#   filter(Count==max(Count))
# 
# Most_prob_lev=merge(Most_prob_lev,data_stat_number2,by="Station")
# 
# Most_prob_lev <- Most_prob_lev %>%
#   mutate(ComfortLevel = as.character(ComfortLevel))
# 
