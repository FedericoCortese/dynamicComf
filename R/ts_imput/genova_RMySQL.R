library(RMySQL)

# create a MySQL connection object
mysqlconnection = dbConnect(RMySQL::MySQL(),
                            dbname='raise_meteo',
                            host='slytherin.mi.imati.cnr.it',
                            #port=3306,
                            user='federico',
                            password='raise_cortese?482Kq')

# hostname: slytherin.mi.imati.cnr.it
# username: federico
# password: raise_cortese?482Kq (K maiuscolo)
# db_name: raise_meteo (tutto minuscolo)
# tabella: ARPAL_FULL (tutto maiuscolo)

dbListTables(mysqlconnection) # displays the tables available in this database.

result = dbSendQuery(mysqlconnection, "select * from ARPAL_FULL_gennaio")
dat1 = fetch(result,n=Inf)
colnames(dat1)

result2 = dbSendQuery(mysqlconnection, "select TEMPO_A, LUOGO,TEMPTRMWC0,UMREIGRWCO from ARPAL_FULL 
                     where LUOGO like 'GENOVA%' order by TEMPO_A,LUOGO")

dat = fetch(result2,n=Inf);dat
