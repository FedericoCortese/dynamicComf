library(RMySQL)

# create a MySQL connection object
mysqlconnection = dbConnect(RMySQL::MySQL(),
                            dbname='raise_meteo',
                            host='slytherin.mi.imati.cnr.it',
                            #port=3306,
                            user='federico',
                            password='raise_cortese?482Kq'
                          )

# hostname: slytherin.mi.imati.cnr.it
# username: federico
# password: raise_cortese?482Kq (K maiuscolo)
# db_name: raise_meteo (tutto minuscolo)
# tabella: ARPAL_FULL (tutto maiuscolo)

dbListTables(mysqlconnection) # displays the tables available in this database.

# result = dbSendQuery(mysqlconnection, "select * from ARPAL_FULL")
# dat1 = fetch(result,n=Inf)
# colnames(dat1)

result2 = dbSendQuery(mysqlconnection, "select 
TEMPO_A, LUOGO,PRECPBIWC1,PATMBRMWCL,
TEMPTRMWC0,RSTORDTWC1,
UMREIGRWCO,VVMDTACWAN from ARPAL_FULL 
                     where LUOGO like 'GENOVA%' order by TEMPO_A,LUOGO")

dat = fetch(result2,n=Inf)
dat$TEMPO_A=as.POSIXct(dat$TEMPO_A, format = "%d/%m/%Y %H:%M")

dat[,-(1:2)]=dat[,-(1:2)]%>%mutate_all(as.numeric)
dat$LUOGO=as.factor(dat$LUOGO)

# sort by TEMPO_A
dat = dat[order(dat$TEMPO_A),]

diff_t=diff(dat$TEMPO_A)



str(dat)
length(unique(dat$LUOGO))

range(dat$TEMPO_A,na.rm = T)

which(is.na(dat$TEMPO_A))

# NAs in TEMPO_A
# dat[311175:311189,]
# dat[323178:323192,]
# dat[335173:335187,]
# dat[369376:369387,]

dat[374274,]
