# valuto leave-one-out cross validation su pezzi di serie
# spazio-temporale

library(dplyr)
library(lubridate)
library(gstat)
library(sp)

# wd weather data (NB: non usare per wind_speed) 
load("cleaned_data.RData")
wd = air_temp
wname = "air_temp"

#source("weather_summ.R") # default media ogni 15 minuti
source("Utils2.R")
tmp = weather_summ(wd)
wd_summ = tmp$mean_x
n_summ = data.frame(tmp$n_x) # numero di valori validi per la media in ciascun intervallo di tempo
rm(tmp)

# scelta periodo di lavoro
init = 90
len = 7
freq = 96
subs = (freq*(init-1)+1) : (freq*(init-1+len))
nsub = len*freq
# riduzione al sottoinsieme
wd_summ_sub = wd_summ[subs,]
n_summ_sub = n_summ[subs,]

# preparazione oggetto ST
locations_wd = locations
coordinates(locations_wd) = ~ longitude + latitude
locations_wd = locations_wd[match(names(wd_summ_sub)[-1], locations_wd$id),]
proj4string(locations_wd) = "+proj=longlat +datum=WGS84"

row.names(locations_wd) = locations_wd$id
# creo oggetto ST
space = list(value = names(wd_summ_sub)[-1])
wd_summ_ST = stConstruct (wd_summ_sub[,-1], space, wd_summ_sub$time, SpatialObj= locations_wd)

# devo creare un oggetto solo spaziale con un indicatore per le date
# dell'intervallo prescelto, partendo da wd_summ_ST che è l'intero
# oggetto spazio-tempo
# creo una lista con i dati spaziali per ciascun tempo e
# la converto in data frame

spdf.lst <- lapply(1:nsub,
                   function(i) {
                     x = wd_summ_ST[,i]
                     x$ti = i # ti = time index
                     x$id = wd_summ_ST@sp$id
                     return(x)}
)

spdf <- do.call(rbind, spdf.lst)
# summary(spdf) # Object of class SpatialPointsDataFrame

# variogramma raggruppato (trattasi di variogramma sui
# residui dalle medie per ogni istante temporale
# che sono diverse, ma spazialmente no). Devo eliminare
# le misurazioni NA da ogni sezione temporale (dati spaziali
# a un dato tempo), perché vgm non le gestisce
# la specifica di dX permette la stima raggruppata (pooled)
# del variogramma nella quale vengono considerate le coppie
# entro ciascun gruppo per il calcolo della semivarianza
# complessiva: dX è la distanza massima tra regressori
# che definisce ciascun gruppo. Il nostro regressore è ti=1, 2, 3 ecc
# Funziona con dX=0 anche se nell'help due misure sono
# nello stesso gruppo se a distanza < dX. 
# NB: l'effetto è lo stesso definendo
# ti come fattore: viene stimata una media spaziale per ciascun
# intervallo di tempo 

vl <- variogram(values ~ ti, spdf[!is.na(spdf$values),], cutoff=25, dX=0)
# plot.numbers mostra il numero di coppie per ogni data distanza
plot(vl, plot.numbers=T, main=paste(wname," ", nsub, " 15-min intervals lumped: ", len, "days"))

# NB: la distanza minima (geodesica) tra due punti per le
# stazioni di Singapore è di 5 km, quindi niente variogramma
# per distanze inferiori (come si vede nel grafico)
# dists = spDists(locations_wd,longlat=T) # matrice distanze
# range(dists[upper.tri(dists)]) # escludo diagonale

# inizializzo variogramma esponenziale con sill 
# e range effettivo. fit.method 6 non pesato sembra meno esigente.
isill = quantile(vl$gamma,0.95)
ihrange = vl$dist[match(T,vl$gamma>=isill)]
vlm <- fit.variogram(vl, model=vgm(isill, "Exp", ihrange/3),fit.method=6) # non assegnando nugget nel modello viene stimato zero
plot(vl, plot.numbers=T, model=vlm, main=paste(wname," ", nsub, " 15-min intervals lumped: ", len, "days"))

# verifica kriging 1: sui 96 quarti d'ora dell'ultimo giorno
# del periodo usato per la stima del variogramma. Per ogni
# quarto d'ora cv leave-one-out di ogni stazione
# raccolgo i risultati in una lista

k.lump.cv = vector("list",freq)
for(i in 1:freq){
  k.lump.cv[[i]] = krige.cv(values ~ 1,spdf[!is.na(spdf$values) & spdf$ti==(nsub-i+1),],model=vlm)
}

# questo comando crea un data frame impilando i data frame
# di ciascun elemento della lista, che hanno tutti la stessa
# struttura
k.lump.cv = do.call(rbind,k.lump.cv)

# osservati e previsti
plot(k.lump.cv$observed,k.lump.cv$var1.pred)
abline(0,1)

# rmse
sqrt(mean(k.lump.cv$residual^2))
plot(k.lump.cv$observed,k.lump.cv$residual)
# sui dati del vento il residuo cresce al crescere del
# valore osservato.

# residuals vs observations

# coordinate dove è stata eseguita la cv e confronto con quelle
# del dataset
plot(k.lump.cv@coords)
points(locations_wd@coords,pch=2)

# la funzione non esegue le previsioni sulle stazioni 
# con tutti i valori mancanti, per le quali occorrono delle
# chiamate a parte di krige()

# come individuare queste stazioni
apply(wd_summ_sub,2,FUN = function(x) all(is.na(x)))




