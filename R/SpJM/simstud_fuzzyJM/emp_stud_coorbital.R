
# 2002AA29 ----------------------------------------------------------------


load("2002AA29_cleaned.Rdata")

# domega=diff(df2002AA29$omega)
# Y$domega=tail(domega,dim(Y)[1])

source("Utils_fuzzyJM.R")
lambda=0.6

str(Y)

fit_2002AA29=fuzzy_jump(Y[,3:4],2,lambda,verbose=TRUE,tol=1e-6)

table(fit_2002AA29$MAP,ground_truth)
table(fit_2002AA29$MAP)


summary(fit_2002AA29$best_S)

res_2002AA29=data.frame(Date=tail((timesY),dim(Y)[1]),
                        tail(df2002AA29,dim(Y)[1]),
                        marg_prob_1=fit_2002AA29$best_S[,1],
                        marg_prob_2=fit_2002AA29$best_S[,2],
                        MAP=as.factor(fit_2002AA29$MAP))

tapply(res_2002AA29$theta, res_2002AA29$MAP, sd)

colors <- rgb(1 - res_2002AA29$marg_prob_1, res_2002AA29$marg_prob_2, 0)
# Plot base con linea grigia
plot(res_2002AA29$Date, res_2002AA29$theta, type = "l", col = "grey", lwd = 2, 
     xlab = "Date", ylab = " ", main = " ")

# Aggiunta dei punti colorati in base al valore di Bull
points(res_2002AA29$Date, res_2002AA29$theta, col = colors, pch = 16, cex = 1)
abline(h=0)


# 2016HO3 --------------------------------------------------------------------

load("2016HO3_cleaned.Rdata")

# domega=diff(df2016HO3$omega)
# Y$domega=tail(domega,dim(Y)[1])

source("Utils_fuzzyJM.R")
lambda=.5
str(Y)

fit_2016HO3=fuzzy_jump(Y,2,lambda,verbose=TRUE,tol=1e-6)

table(fit_2016HO3$MAP,ground_truth)


# 100006174 ---------------------------------------------------------------

load("100006174_cleaned.Rdata")

source("Utils_fuzzyJM.R")
lambda=10

str(Y)

# Convert indicators into factors
Y[,5]=as.factor(Y[,5])
Y[,6]=as.factor(Y[,6])
Y[,7]=as.factor(Y[,7])
Y[,8]=as.factor(Y[,8])
Y[,9]=as.factor(Y[,9])

str(Y)

fit_100006174=fuzzy_jump(Y,4,lambda,verbose=TRUE,tol=1e-6)

summary(fit_100006174$best_S)

res_1000016174=data.frame(Date=timesY,
                          tail(df100006174,dim(Y)[1]),
                          marg_prob_1=fit_100006174$best_S[,1],
                          marg_prob_2=fit_100006174$best_S[,2],
                          marg_prob_3=fit_100006174$best_S[,3],
                          marg_prob_4=fit_100006174$best_S[,4],
                          MAP=as.factor(fit_100006174$MAP))




# Definizione colori personalizzati
color_palette <- c("1" = "black", 
                   "2" = "red", 
                   "3" = "blue", 
                   "4" = "green")

# Plot interattivo
plot_ly(data = res_1000016174,
        text = ~MAP) %>%
  add_trace(x = ~Date, 
            y = ~theta, 
            type = 'scatter', 
            mode = 'lines', 
            line = list(color = 'lightgray', width = 1),
            name = 'Trend') %>%
  add_trace(x = ~Date, 
            y = ~theta, 
            type = 'scatter', 
            mode = 'markers', 
            color = ~MAP
            # ,
            # colors = color_palette
            ) %>%
  layout(title = " ",
         xaxis = list(title = "Date"),
         yaxis = list(title = " "),
         legend = list(title = list(text = "Colori MAP")))
