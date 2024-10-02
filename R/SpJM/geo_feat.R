library(geojsonR)
file_js = FROM_GeoJson(url_file_string = "Singapore_nodes_100m.geojson")
str(file_js)
file_js_1=file_js[[1]]

coords_temp=lapply(file_js_1,function(x)x$geometry$coordinates)
coords=matrix(unlist(coords_temp),ncol=2,byrow=TRUE)

building_count_temp=lapply(file_js_1,function(x)x$properties$`Building Count`)
building_count=unlist(building_count_temp)

building_view_mean_temp=lapply(file_js_1,function(x)x$properties$`Building View Mean`)
building_view_mean=unlist(building_view_mean_temp)

complexity_mean_temp=lapply(file_js_1,function(x)x$properties$`Complexity Mean`)
complexity_mean=unlist(complexity_mean_temp)

green_view_mean_temp=lapply(file_js_1,function(x)x$properties$`Green View Mean`)
green_view_mean=unlist(green_view_mean_temp)

road_view_mean_temp=lapply(file_js_1,function(x)x$properties$`Road View Mean`)
road_view_mean=unlist(road_view_mean_temp)

sky_view_mean_temp=lapply(file_js_1,function(x)x$properties$`Sky View Mean`)
sky_view_mean=unlist(sky_view_mean_temp)

geo_feat=data.frame(longitude=coords[,1],
                    latitude=coords[,2],
                    building_count,
                    building_view_mean,
                    complexity_mean,
                    green_view_mean,
                    road_view_mean,
                    sky_view_mean)

save(geo_feat,file="geo_feat.RData")

