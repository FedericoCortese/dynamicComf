# simulazione dal modello di Potts (generalizzazione a più colori del modello di Ising)

require("potts")

set.seed(2024)

# Inizializzazione a caso dell'oggetto (e trasformazione in formato raw) che rappresenta la matrice degli 
# stati (4 colori) in ciascuna posizione (pixel, immagine nr*nc)

ncolor = as.integer(4) # transizione di fase continua per 1 <= ncolor <= 4
nr = 100
nc = 100
init <- matrix(sample(ncolor, nr * nc, replace = TRUE), nrow = nr, ncol=nc)
init <- packPotts(init, ncol = ncolor)

# valori dei parametri del modello
# il parametro beta regola la similitudine tra pixel vicini prossimi, più beta è alto più sono simili
# il valore specifico log(1 + sqrt(ncolor)) rappresenta il punto critico della transizione di fase, cioè
# con valori inferiori si ottengono configurazioni disordinate, con valori superiori configurazioni monocolore
beta <- log(1 + sqrt(ncolor))

# il modello implementato è quello della famiglia esponenziale con statistica di dimensione ncolor + 1
# e logverosimiglianza sum_i t_i(x)*theta_i - m(theta), con x la matrice dei colori e, posto
# con d=ncolor +1, theta_d coincide con beta e moltiplica la componente di interazione t_d(x) che vale 0.5 * n. coppie di vicini con lo stesso colore
# gli altri t_i(x) sono il numero totale di pixel di colore i, inizializzando a zero i corrispondenti parametri non si privilegia nessun
# colore in particolare, altrimenti mettere valori non nulli 

theta <- c(rep(0, ncolor), beta)

# output della simulazione: il numero di iterazioni è dato dal numero di batches per la lunghezza di ciascun batch
# con il parametro nspac=1, un batch è un gruppo di iterazioni consecutive. L'output out$batch è di default (ma è modificabile)
# la media delle statistiche canoniche della famiglia esponenziale t_i(x)
out <- potts(init, param=theta, nbatch = 200 , blen=1, nspac=1)
image(out$final) # NB è solo lo stato finale della catena di Markov che simula il modello di Potts!

# verificare la stabilizzazione della statistica t_d(x) per vedere se il sistema ha raggiunto lo stato stazionario
plot(out$batch[,ncolor+1],type="l")

# simulazione MCMC di una sequenza di immagini
nit=10
sim= vector("list",nit)
sim[[1]] = out$final
for(it in 2:nit){
  out = potts(sim[[it-1]], param=theta, nbatch = 1 , blen=1, nspac=1)
  sim[[it]] = out$final
}


# il grafico dei colori dovrebbe mostrare una equiproporzione con i primi ncolor parametri di theta pari a zero
matplot(out$batch[,-(ncolor+1)],type="l")


