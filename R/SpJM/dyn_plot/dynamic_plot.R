load("C:/Users/federico/OneDrive - CNR/Comfort - HMM/simres_STJM_dist/emp_res_STJMdist_singapore_including_UTCI(2).RData")

data_stat_number2=data_stat_number[,-1]
data_stat_number2$Station=paste("S",data_stat_number2$m)
Sys.setenv(TZ="GMT")
S_long <- S_est %>%
  pivot_longer(cols = starts_with("S"), names_to = "Station", values_to = "ComfortLevel")

# Merge the spatial data with comfort levels data
S_long <- S_long %>%
  left_join(data_stat_number2, by = "Station")

S_long$time=as.POSIXct(S_long$time,tz="GMT")+hours(8)

new_times=as.POSIXct(timesY,tz="GMT")+hours(8)

for (input in 1:length(new_times)){
  
  filtered_data <- S_long %>%
    filter(time %in% new_times[input])
  
  safe_timestamp <- format(new_times[input], "%Y-%m-%d_%H-%M-%S")
  
  # Create the leaflet map
  pp <- leaflet(data = filtered_data, options = leafletOptions(zoomControl = FALSE)) %>%
    addTiles() %>%
    setView(
      lng = mean(filtered_data$longitude, na.rm = TRUE), 
      lat = mean(filtered_data$latitude, na.rm = TRUE) + 0.03, 
      zoom = 11.4
    ) %>%
    addCircleMarkers(
      lng = ~longitude, lat = ~latitude,
      color = ~unname(comfort_colors[as.character(ComfortLevel)]),
      fillColor = ~unname(comfort_colors[as.character(ComfortLevel)]),
      fillOpacity = 1,
      radius = 10,
      stroke = TRUE,
      weight = 1,
      opacity = 1,  
      label = ~paste0("Station: ", Station, "<br>Comfort Level: ", ComfortLevel),
      labelOptions = labelOptions(style = list("color" = "black"), direction = "auto")
    ) %>%
    addLegend(
      position = "topright",
      colors = unlist(comfort_colors),
      labels = c("Cool", "Neutral", "Hot"),
      title = "Comfort Levels",
      opacity = 0.7
    ) %>%
    addControl(
      html = paste0("<h3 style='text-align:center;'>Time: ", new_times[input], "</h3>"),
      position = "topleft"
    )
  
  # Step 1: Save the leaflet map as an HTML file with a safe file name
  html_file <- paste0("dynamic_plot_", safe_timestamp, ".html")
  saveWidget(pp, file = html_file)
  png_file <- paste0(input,"_dynamic_plot_", safe_timestamp, ".png")
  webshot2::webshot(url = html_file, file = png_file, vwidth = 600, vheight = 500)
  
}

