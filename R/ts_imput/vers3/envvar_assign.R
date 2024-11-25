# assegno valore variabili ambientali ai feedback con kriging ordinario
# raggruppato
library(dplyr)
library(lubridate)
library(gstat)
library(sp)
load("cleaned_data.RData")

wd = air_temp
wname = "air_temp"

#source("weather_summ.R") # default media ogni 15 minuti
source("Utils2.R")
tmp = weather_summ(wd, window="15 minutes")
wd_summ = tmp$mean_x
n_summ = data.frame(tmp$n_x) # numero di valori validi per la media in ciascun intervallo di tempo
rm(tmp)

# riduco il file a misure quando il partecipante risponde al questionario e si trova fuori
# ws_latitude non è assegnata quando non è interrogato, altre variabili individiuali vengono misurate lo stesso
#
cozie_q = cozie_train[!is.na(cozie_train$ws_latitude) & cozie_train$q_location=="Outdoor",]
# le coordinate delle posizioni del test sono in ws_longitude e
# ws_latitude, come tempo prendiamo ws_timestamp_start

# per ogni timestamp entro un dato quarto d'ora devo creare
# un oggetto con le coordinate da assegnare a "newdata"
# nella funzione krige
# trovo gli intervalli di 15 minuti in wd_summ$time
# il default della funzione è chiuso a sinistra e aperto a destra
# chiudo ultimo intervallo (tt = testtime)
ttindex = findInterval(cozie_q$ws_timestamp_start,wd_summ$time,rightmost.closed=T)

# per assegnare i valori di temperatura
# 1) per ogni valore unico di ttindex (che identifica data e int. di tempo)
# raccolgo le coordinate delle prove in quell'intervallo
# 2) stimo variogramma usando tot giorni di cui uno contiene
# l'intervallo, per esempio 7 giorni all'indietro partendo da lì
# 3) faccio kriging sulle coordinate

# numero di intervalli al giorno
freq = 96 # 96 quarti d'ora al giorno
n = nrow(wd_summ)

# periodo minimo iniziale per la stima del variogramma
# (per esempio 7 giorni con media ogni 15 minuti)
len = 7
subs = 1 : (freq*len)
nsub = freq*len

# preparazione per oggetto ST, con le posizioni delle stazioni
# presenti in wd_summ e nell'ordine con cui compaiono
# in wd_summ
locations_wd = locations
coordinates(locations_wd) = ~ longitude + latitude
locations_wd = locations_wd[match(names(wd_summ)[-1], locations_wd$id),]
proj4string(locations_wd) = "+proj=longlat +datum=WGS84"
row.names(locations_wd) = locations_wd$id

# raccolta temperature da assegnare ai feedback con varianze di previsione
y = matrix(NA,nrow(cozie_q),2)

for(ind in unique(ttindex)){
  print(ind)
  # finestra mobile appena possibile
  if(ind>freq*len) subs = (ind-freq*len+1):ind
  # riduzione al sottoinsieme
  wd_summ_sub = wd_summ[subs,]
  n_summ_sub = n_summ[subs,]
  
  # creo oggetto ST
  space = list(value = names(wd_summ_sub)[-1])
  wd_summ_ST = stConstruct (wd_summ_sub[,-1], space, wd_summ_sub$time, SpatialObj= locations_wd)
  
  # devo creare un oggetto solo spaziale con un indicatore per le date
  # dell'intervallo prescelto, partendo da wd_summ che è l'intero
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
  # tempo ma è molto più veloce che definendo ti come fattore.
  
  # la prima stima del variogramma è usata per tutti i valori di ind
  # fino a freq* len
  if(ind==ttindex[1] & ind<=freq*len){
    # cutoff basato su esame di dists = spDists(locations_wd,longlat=T)
    # e visualizzazione di alcuni variogrammi
    vl <- variogram(values ~ ti, spdf[!is.na(spdf$values),], cutoff=25, dX=0)
    isill = quantile(vl$gamma,0.95)
    ihrange = vl$dist[match(T,vl$gamma>=isill)]
    vlm <- fit.variogram(vl, model=vgm(isill, "Exp", ihrange/3),fit.method=6) # non assegnando nugget nel modello viene stimato zero
    # plot(vl, plot.numbers=T, model=vlm, main=paste(wname, " ", nsub, " 15-min intervals lumped: ", len, "days"))
  }
  if(ind>freq*len){ # devo rifare la stima quando la finestra inizia a muoversi
    vl <- variogram(values ~ ti, spdf[!is.na(spdf$values),], cutoff=25, dX=0)
    isill = quantile(vl$gamma,0.95)
    ihrange = vl$dist[match(T,vl$gamma>=isill)]
    vlm <- fit.variogram(vl, model=vgm(isill, "Exp", ihrange/3),fit.method=6) # non assegnando nugget nel modello viene stimato zero
  }
  # ora faccio il kriging su tutte le posizioni degli individui
  # nell'intervallo di tempo definito da "ind"
  # usando SpatialPoints(), che richiede le coordinate 
  # e il metodo di calcolo delle distanze geodetiche
  slocs <- SpatialPoints(cozie_q[which(ttindex==ind),c("ws_longitude","ws_latitude")],
                         proj4string=CRS(proj4string(spdf)))
  # gridded(sing.bbox.grid) <- TRUE # non ho capito se occorre o no
  
  # kriging
  # nei primi freq*len intervalli di tempo il sottoinsieme di wd_summ
  # selezionato da wd_summ_sub sono le prime freq*len righe 
  # allora quando spdf$ti==ind sto selezionando il quarto d'ora in wd_summ 
  # corrispondente a quello del feedback
  # quando la finestra comincia a scorrere, il quarto d'ora del feedback
  # in wd_summ_sub è sempre l'ultimo
  if(ind <=freq*len) 
    k.lump = krige(values ~ 1,spdf[!is.na(spdf$values) & spdf$ti==ind,], newdata=slocs,model=vlm)
  else
    k.lump = krige(values ~ 1,spdf[!is.na(spdf$values) & spdf$ti==nsub,], newdata=slocs,model=vlm)
  # posiziono le temperature trovate in y per poi attaccarle a cozie_q
  y[which(ttindex==ind),] = as.matrix(k.lump@data)
}

stop("termine kriging, eseguire il resto a mano")

# aggiungo variabile amb. e varianze kriging a cozie_q
# scrivere nome a mano
cozie_q.wd = data.frame(cozie_q.wd,RH=y[,1],RH.v=y[,2])

# controllo se marginalmente la var. amb ha effetto sulla preferenza
prefs = aggregate(cozie_q.wd$q_thermal_preference,by=list(temperature= round(cozie_q.wd$air_temp,digits = 0)), FUN="table")
barplot(height = t(prefs$x),beside=FALSE, names.arg=prefs$temperature,legend.text = c("cooler","ok","warmer"),args.legend = list(x = "topleft"),xlab="temperature")
title("absolute frequencies")
prefs$x[,1:3] = 100 * prefs$x[,1:3] / apply(prefs$x[,1:3],1,sum)
barplot(height = t(prefs$x),beside=FALSE, names.arg=prefs$temperature,xlab="temperature")
title("percentages")


# grafico delle temperature imputate vs timestamp
plot(wd_summ$time[ttindex], cozie_q.wd$RH,type="l")

# grafico imputate vs stazione
plot(wd_summ$S24[ttindex],cozie_q.wd$RH)
