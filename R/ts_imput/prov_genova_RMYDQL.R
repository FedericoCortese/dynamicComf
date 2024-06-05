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

result2 = dbSendQuery(mysqlconnection, "select * from ARPAL_FULL where PV like 'GE'")
dat = fetch(result2,n=Inf)

dat$TEMPO_DA=dmy_hm(dat$TEMPO_DA)
dat$TEMPO_A=dmy_hm(dat$TEMPO_A)

which(is.na(dat$TEMPO_DA))
which(is.na(dat$TEMPO_A))

dat[,-(1:6)]=dat[,-(1:6)]%>%mutate_all(as.numeric)

# Remove space and -
dat$LUOGO=gsub(" ","",dat$LUOGO)
dat$LUOGO=gsub("-","",dat$LUOGO)

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
dat=dat%>% 
  right_join(data.frame(TEMPO_DA=eqdist_times),
             by="TEMPO_DA")%>%
  arrange(TEMPO_DA)

datGE=dat

# Save
save(datGE,file="prov_genova.RData")

# Provincia di Genova -----------------------------------------------------

load("prov_genova.RData")

dat=datGE

locations=readxl::read_xlsx("locations_genova.xlsx",sheet=3)
colnames(locations)[5:7]=c("Height","Longitude","Latitude")
locations$Longitude=as.numeric(locations$Longitude)

# Merge duplicate stations
locations=locations%>%group_by(`NOME STAZIONE`)%>%
  summarise(Longitude=mean(Longitude),
            Latitude=mean(Latitude),
            Height=mean(Height))

# Remove space and -
locations$`NOME STAZIONE`=gsub(" ","",locations$`NOME STAZIONE`)
locations$`NOME STAZIONE`=gsub("-","",locations$`NOME STAZIONE`)

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


# Air Temperature ---------------------------------------------------------

# Wide format only for TEMPTRMWC0
dat_temp_wide = dat%>%select(TEMPO_DA,LUOGO,TEMPTRMWC0)%>%
  pivot_wider(names_from = LUOGO,values_from = TEMPTRMWC0)

# summary(dat_temp_wide)

# Percentage of NA for each station
na=apply(dat_temp_wide[,-1],2,function(x){sum(is.na(x))/length(x)})

# Remove stations with more than 50% of NA
dat_temp_wide=dat_temp_wide[,c(1,which(na<0.5)+1)]

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

# windows()
# Amelia::missmap(dat_temp_wide,main="Temperature",margins=c(10,1),x.cex=.5)
# 
# windows()
# Amelia::missmap(dat_temp_wide[,c("S. MARGHERITA LIGURE","MONTE PORTOFINO")],
#                 main="Temperature",margins=c(10,1),x.cex=.5)
# 
# windows()
# Amelia::missmap(dat_temp_wide[,c("MONTE PENNELLO")],
#                 main="Temperature",margins=c(10,5),x.cex=.5)

# Remove stations with many consecutive NAs

dat_temp_wide2=select(dat_temp_wide,subset=-c(`ALPEVOBBIA`,`MONTEPORTOFINO`,`S.MARGHERITALIGURE`,
                                              `BUSALLA`,`FIORINO`))
dat_temp_wide2=dat_temp_wide2[(1:27000),]

windows()
Amelia::missmap(dat_temp_wide2,
                main="Temperature",margins=c(10,5),x.cex=.5)

# Plot these residual stations as in "map" above
locations2=locations[locations$`NOME STAZIONE` %in% colnames(dat_temp_wide2),]

map2=leaflet() %>% 
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

source("vers2/Utils2.R")
# locations2=locations2[,c("NOME STAZIONE","Longitude","Latitude")]
# 
# x_data=dat_temp_wide2

dat_temp_wide3=weightdist_imp(dat_temp_wide2,locations2[,c("NOME STAZIONE","Longitude","Latitude")])
na=apply(dat_temp_wide3[,-1],2,function(x){sum(is.na(x))/length(x)})
na*100

# Relative humidity -------------------------------------------------------

dat_rh_wide = dat%>%select(TEMPO_DA,LUOGO,UMREIGRWCO)%>%
  pivot_wider(names_from = LUOGO,values_from = UMREIGRWCO)
na_rh=apply(dat_rh_wide[,-1],2,function(x){sum(is.na(x))/length(x)}); na*100
dat_rh_wide=dat_rh_wide[,c(1,which(na_rh<0.5)+1)]
dat_rh_wide2=dat_rh_wide[(1:27000),]

commonstat=intersect(colnames(dat_rh_wide2),colnames(dat_temp_wide3))
dat_rh_wide2=dat_rh_wide2[,commonstat]
dat_rh_wide2=data.frame(time=dat_temp_wide3$time,dat_rh_wide2)
dat_temp_wide3=dat_temp_wide3[,commonstat]
dat_temp_wide3=data.frame(time=dat_rh_wide2$time,dat_temp_wide3)

windows()
Amelia::missmap(dat_rh_wide2)

locations3=locations2[locations2$`NOME STAZIONE` %in% commonstat,]

map3=leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(data=locations3,
                   radius = 4,
                   color = 'red',
                   stroke = FALSE, fillOpacity = 1,
                   popup = ~paste("<br>Lat:",Latitude, "<br>Lon:", Longitude
                                  , "<br>Height (m):", Height, "<br>Station:", `NOME STAZIONE`
                   )
  )
map3

dat_rh_wide3=weightdist_imp(dat_rh_wide2,locations3[,c("NOME STAZIONE","Longitude","Latitude")])

# Sort columns in alphabetical order
dat_rh_wide3=dat_rh_wide3[,order(colnames(dat_rh_wide3))]
dat_rh_wide3=dat_rh_wide3%>% relocate(time, .before = everything())

dat_temp_wide3=dat_temp_wide3[,order(colnames(dat_temp_wide3))]
dat_temp_wide3=dat_temp_wide3%>% relocate(time, .before = everything())

save(locations3,dat_rh_wide3,dat_temp_wide3,file="prov_genova_tempANDrh.RData")


