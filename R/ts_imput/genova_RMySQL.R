library(RMySQL)
library(dplyr)
library(tidyr)
library(geojsonio)
library(leaflet)
library(leaflet.extras)
library(lubridate)
library(htmltools)
Sys.setlocale("LC_TIME")

# read db -----------------------------------------------------------------
# create a MySQL connection object
mysqlconnection = dbConnect(RMySQL::MySQL(),
                            dbname='raise_meteo',
                            host='slytherin.mi.imati.cnr.it',
                            #port=3306,
                            user='federico',
                            password='raise_cortese?482Kq'
                          )

dbListTables(mysqlconnection) # displays the tables available in this database.


# result2 = dbSendQuery(mysqlconnection, "select 
# TEMPO_A, LUOGO,PRECPBIWC1,PATMBRMWCL,
# TEMPTRMWC0,RSTORDTWC1,
# UMREIGRWCO,VVMDTACWAN from ARPAL_FULL 
#                      where LUOGO like 'GENOVA%' order by TEMPO_A,LUOGO")

result2 = dbSendQuery(mysqlconnection, "select * from ARPAL_FULL where LUOGO like 'GENOVA%'")

dat = fetch(result2,n=Inf)
dat$TEMPO_DA=dmy_hm(dat$TEMPO_DA)
dat$TEMPO_A=dmy_hm(dat$TEMPO_A)

which(is.na(dat$TEMPO_DA))
which(is.na(dat$TEMPO_A))


dat[,-(1:6)]=dat[,-(1:6)]%>%mutate_all(as.numeric)
dat$LUOGO=as.factor(dat$LUOGO)
# sort by TEMPO_DA
dat = dat[order(dat$TEMPO_DA),]

str(dat)
summary(dat)

which(is.na(dat$TEMPO_DA))

# All hours
eqdist_times=seq(from=dat$TEMPO_DA[1],
                 to=tail(dat$TEMPO_DA,1),
                 by="1 hour")

# Merge data
dat2=dat%>% 
  right_join(data.frame(TEMPO_DA=eqdist_times),
             by="TEMPO_DA")%>%
  arrange(TEMPO_DA)

# Save
save(dat,file="genova.RData")


# Genova ------------------------------------------------------------------

load("genova.RData")

locations=readxl::read_xlsx("locations_genova.xlsx",sheet=2)
colnames(locations)[5:7]=c("Height","Longitude","Latitude")

map=leaflet(data=locations) %>% 
  addTiles() %>%
  addCircleMarkers(data=locations,
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height
                   )
  )

#windows()
map

# map=map %>% addLegend("bottomright", 
#                       colors = c("red","black"), 
#                       labels = c("Stations","Participants"), title = "Legend")

# Wide format only for TEMPTRMWC0
dat_temp_wide = dat%>%select(TEMPO_DA,LUOGO,TEMPTRMWC0)%>%
pivot_wider(names_from = LUOGO,values_from = TEMPTRMWC0)

summary(dat_temp_wide)

# Percentage of NA for each station
na=apply(dat_temp_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

# Remove stations with more than 50% of NA
dat_temp_wide=dat_temp_wide[,c(1,which(na<0.5)+1)]

Amelia::missmap(dat_temp_wide,main="Temperature")

# Function to count gaps
count_gaps=function(x){
  gaps=list()
  gaps[[1]]=x[,1]
  for(i in 2:ncol(x)){
    d_na <- as.numeric(is.na(x[,i]))
    Nal=rle(d_na)
    tosave=Nal$values==1
    Nals=Nal$lengths[tosave]
    gaps[[i]]=Nals
  }
  names(gaps)=colnames(x)
  return(gaps)
}

Na_count=count_gaps(dat_temp_wide)
Na_count

# Plot these residual stations as in "map" above
locations2=locations[locations$`NOME STAZIONE` %in% colnames(dat_temp_wide),]

map2=leaflet(data=locations) %>% 
  addTiles() %>%
  addCircleMarkers(data=locations2,
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height, "<br>Station:", `NOME STAZIONE`
                   )
  )

map2

# Wide format for UMREIGRWCO
dat_umreig_wide = dat%>%select(TEMPO_DA,LUOGO,UMREIGRWCO)%>%
pivot_wider(names_from = LUOGO,values_from = UMREIGRWCO)

summary(dat_umreig_wide)

# Keep stations with less than 50% of NA
na=apply(dat_umreig_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

dat_umreig_wide=dat_umreig_wide[,c(1,which(na<0.5)+1)]

colnames(dat_umreig_wide)
colnames(dat_temp_wide)

# Amelia::missmap(dat_umreig_wide,main="Relative humidity")

# Wide format for PRECPBIWC1
dat_prec_wide = dat%>%select(TEMPO_DA,LUOGO,PRECPBIWC1)%>%
pivot_wider(names_from = LUOGO,values_from = PRECPBIWC1)

# Keep stations with less than 50% of NA
na=apply(dat_prec_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

dat_prec_wide=dat_prec_wide[,c(1,which(na<0.5)+1)]

summary(dat_prec_wide)

colnames(dat_temp_wide)
colnames(dat_prec_wide)

# Wide format for PATMBRMWCL
dat_patm_wide = dat%>%select(TEMPO_DA,LUOGO,PATMBRMWCL)%>%
pivot_wider(names_from = LUOGO,values_from = PATMBRMWCL)

# keep stations with less than 50% of NA
na=apply(dat_patm_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

dat_patm_wide=dat_patm_wide[,c(1,which(na<0.5)+1)]

summary(dat_patm_wide)

colnames(dat_patm_wide)
colnames(dat_temp_wide)

# Wide format for RSTORDTWC1
dat_rstord_wide = dat%>%select(TEMPO_DA,LUOGO,RSTORDTWC1)%>%
pivot_wider(names_from = LUOGO,values_from = RSTORDTWC1)

# keep stations with less than 50% of NA
na=apply(dat_rstord_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

dat_rstord_wide=dat_rstord_wide[,c(1,which(na<0.5)+1)]

summary(dat_rstord_wide)

colnames(dat_rstord_wide)
colnames(dat_temp_wide)

# Wide format for VVMDTACWAN
dat_vvmdtacwan_wide = dat%>%select(TEMPO_DA,LUOGO,VVMDTACWAN)%>%
pivot_wider(names_from = LUOGO,values_from = VVMDTACWAN)

# keep stations with less than 50% of NA
na=apply(dat_vvmdtacwan_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100

dat_vvmdtacwan_wide=dat_vvmdtacwan_wide[,c(1,which(na<0.5)+1)]

summary(dat_vvmdtacwan_wide)

colnames(dat_vvmdtacwan_wide)
colnames(dat_temp_wide)

# Scatterplot dat_temp_wide vs dat_prec_wide for each station

plot(dat_temp_wide$`GENOVA - PONTEDECIMO`,dat_prec_wide$`GENOVA - PONTEDECIMO`,
     xlab="Temperature",ylab="Precipitation",main="GENOVA - PONTEDECIMO")
cor(dat_temp_wide$`GENOVA - PONTEDECIMO`,dat_prec_wide$`GENOVA - PONTEDECIMO`, 
    use="complete.obs")
