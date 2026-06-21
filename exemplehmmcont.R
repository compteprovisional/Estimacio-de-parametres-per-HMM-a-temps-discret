# exemple d'estimació de paràmetres de un model de Markov Ocult amb observacions contínues

library(mhsmm)

Npasses <- c(1000, 10000, 20000, 30000, 40000, 50000) # diferents longituds de cadena
N <- length(Npasses)

# en aquestes llistes guardarem els resultats
modelsestimats <- vector(mode = "list", length = N) # models (paràmetres de la matriu d'emissió i transició)
cadests <- vector(mode = "list", length = N) # cadenes màximes versemblants pel model estimat
props <- numeric(N) # proporció d'estats de la cadena d'estats estimada que són correctes

inici <- c(0, 1) # pi
matriutrans <- matrix(c(0.9, 0.2, 0.1, 0.8), ncol = 2) # A
paramshmm <- list(mu = c(120, 140), sigma = c(20, 25)) # paràmetres de les distribucions d'emissió


inici0 <- c(0, 1) # estimació inicial de la distribució inicial
A0 <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2) # "estimació" inicial de la matriu de transició


for (i in 1:N) {
  modelreal <- hmmspec(init = inici, trans = matriutrans, parms.emission = paramshmm, dens.emission = dnorm.hsmm)
  
  cadobs <- simulate(modelreal, nsim = Npasses[i], rand.emission = rnorm.hsmm, seed = 123)
  
  ON <- cadobs$x # cadena d'observacions
  passat <- cadobs$s # cadena d'estats
  
  alcmitj <- mean(ON) # les estimacions inicials de les observacions que emprarem són la mitjana mostral i la desviació estàndard pels dos estats possibles 
  sdinicial <- sd(ON)
  paramshmm0 <- list(mu = c(alcmitj, alcmitj), sigma = c(sdinicial, sdinicial))
  
  modelinicial <- hmmspec(init = inici0, trans = A0, parms.emission = paramshmm0, dens.emission = dnorm.hsmm)
  
  modelestimat <- hmmfit(cadobs, modelinicial, mstep = mstep.norm, tol = 1e-10, maxit = 100)
  
  modelsestimats[[i]] <- modelestimat$model
  cadests[[i]] <- modelestimat$yhat
  props[[i]] <- sum(cadests[[i]] == passat)/Npasses[i]
}

modelsestimats

plot(Npasses, props, type = "l")

cosa <- Npasses*(1 - props)









