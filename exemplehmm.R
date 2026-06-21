# exemple on tiram 3 monedes i vem si és cara o creu

library(HMM)

set.seed(1)

Npasses <- 1000 # quantitat d'observacions/estats

Niters <- 1000 # quantitat d'iteracions de EM

estats <- c(1, 2, 3) # nombre de la moneda
observacions <- c(1, 2) # cara, creu

trans <- matrix(c(0.7, 0.15, 0.1, 0.2, 0.75, 0, 0.1, 0.1, 0.9), nrow = 3) # A
probas <- matrix(c(0.9, 0.5, 0.1, 0.1, 0.5, 0.9), ncol = 2) # B
inici <- c(0.2, 0.3, 0.5) # pi

hmmreal <- initHMM(estats, observacions, startProbs = inici, transProbs = trans, emissionProbs = probas)

cadobs <- simHMM(hmmreal, Npasses)

ON <- cadobs$observation

passat <- cadobs$states

sum(passat == 3)/1000

forwardsreal <- exp(forward(hmmreal, ON))

POLambdareal <- sum(forwardsreal[, Npasses])


# Primer estimarem variant les probabilitats d'emissió. Sempre emprarem les següents probabilitats de transició i la distribució inicial:

A0 <- matrix(1/3, nrow = 3, ncol = 3)
pi0 <- c(0.2, 0.3, 0.5)

# Aquestes són les matrius de transició no aleatòries que emprarem:

B0 <- matrix(1/2, nrow = 3, ncol = 2)
B1 <- matrix(c(replicate(3, sum(ON == 1)/Npasses), replicate(3, sum(ON == 2)/Npasses)), ncol = 2)
B2 <- matrix(c(0.6, 0.5, 0.4, 0.4, 0.5, 0.6), ncol = 2)

hmm0B0 <- initHMM(estats, observacions, startProbs = pi0, transProbs = A0, emissionProbs = B0)
P0B0 <- sum(exp(forward(hmm0B0, ON)[, Npasses]))

hmm0B1 <- initHMM(estats, observacions, startProbs = pi0, transProbs = A0, emissionProbs = B1)
P0B1 <- sum(exp(forward(hmm0B1, ON)[, Npasses]))

hmm0B2 <- initHMM(estats, observacions, startProbs = pi0, transProbs = A0, emissionProbs = B2)
P0B2 <- sum(exp(forward(hmm0B2, ON)[, Npasses]))


estimarmodel <- function(OBS, NiterEM, estatsmodel, observacionsmodel, paramsinicials) {
  hmm0 <- initHMM(estatsmodel, observacionsmodel, startProbs = paramsinicials[[1]], transProbs = paramsinicials[[2]], emissionProbs = paramsinicials[[3]])
  
  est <- baumWelch(hmm0, OBS, maxIterations = NiterEM, delta = 1E-100)
  
  est_inici <- est$hmm$startProbs
  est_trans <- est$hmm$transProbs
  est_emission <- est$hmm$emissionProbs
  
  return(initHMM(estatsmodel, observacionsmodel, startProbs = est_inici, transProbs = est_trans, emissionProbs = est_emission))
}

modelestimatB0 <- estimarmodel(ON, Niters, estats, observacions, list(pi0, A0, B0))

modelestimatB1 <- estimarmodel(ON, Niters, estats, observacions, list(pi0, A0, B1))

modelestimatB2 <- estimarmodel(ON, Niters, estats, observacions, list(pi0, A0, B2))

PO0 <- sum(exp(forward(modelestimatB0, ON)[, Npasses]))

PO1 <- sum(exp(forward(modelestimatB1, ON)[, Npasses]))

PO2 <- sum(exp(forward(modelestimatB2, ON)[, Npasses]))





generarmatriuB <- function() {
  cosa <- sort(runif(3), decreasing = TRUE) # generam 3 uniformes(0, 1), que seran les probabilitat de cara de les 2 monedes
  cosa2 <- 1 - cosa # probabilitat de creu
  
  vect <- cbind(cosa, cosa2)
  
  B <- matrix(vect, ncol = 2, nrow = 3)
  
  return (B)
}


Bmostres <- 50

aleatoriesB <- replicate(Bmostres, generarmatriuB())

modelsaleatorisB <- vector(mode = "list", length = Bmostres)

for (i in 1:Bmostres) {
  modelsaleatorisB[[i]] <- estimarmodel(ON, Niters, estats, observacions, list(pi0, A0, aleatories[, , i]))
}


POaleatorisB <- numeric(Bmostres)

for (i in 1:Bmostres) {
  POaleatorisB[i] <- sum(exp(forward(modelsaleatorisB[[i]], ON)[, Npasses]))
}


hist(POaleatorisB, breaks = 4, main = "Versemblances per matrius B aleatòries", xlab = "Versemblança", ylab = "Freqüència")
abline(v = mean(POaleatorisB), col = "red", lwd = 3)


vers <- which.max(POaleatorisB)
aleatoriesB[, , vers]
POaleatorisB[vers]
modelsaleatorisB[[vers]]


modelstipus1B <- modelsaleatorisB[POaleatorisB < 1e-262]
modelstipus2B <- Filter(function(x) x$emissionProbs[2,1] > 0.5, modelstipus1B)
modelstipus1B <- Filter(function(x) x$emissionProbs[2,1] <= 0.5, modelstipus1B)

modelstipus3B <- modelsaleatorisB[POaleatorisB >= 1e-262]




# passem ara a variar la matriu de transició

generarmatriuA <- function() {
  cosa <- runif(9)
  A <- matrix(cosa, ncol = 3, nrow = 3)
  
  for (i in 1:3) {
    A[i, ] <- A[i, ]/sum(A[i, ])
  }
  
  return (A)
}

A1 <- matrix(c(0.6, 0.2, 0.2, 0.2, 0.6, 0.2, 0.2, 0.2, 0.6), ncol = 3, nrow = 3)

hmm0A1 <- initHMM(estats, observacions, startProbs = pi0, transProbs = A1, emissionProbs = B0)
P0A1 <- sum(exp(forward(hmm0A1, ON)[, Npasses]))

modelestimatA1 <- estimarmodel(ON, Niters, estats, observacions, list(pi0, A1, B0))
PA1 <- sum(exp(forward(modelestimatA1, ON)[, Npasses]))

# Ara aleatòries

set.seed(1)

Amostres <- 50

aleatoriesA <- replicate(Amostres, generarmatriuA())

modelsaleatorisA <- vector(mode = "list", length = Amostres)

for (i in 1:Amostres) {
  modelsaleatorisA[[i]] <- estimarmodel(ON, Niters, estats, observacions, list(pi0, aleatoriesA[, , i], B0))
}


POaleatorisA <- numeric(Amostres)

for (i in 1:Amostres) {
  POaleatorisA[i] <- sum(exp(forward(modelsaleatorisA[[i]], ON)[, Npasses]))
}


hist(POaleatorisA, breaks = 4, main = "Versemblances per matrius A aleatòries", xlab = "Versemblança", ylab = "Freqüència")
abline(v = mean(POaleatorisA), col = "red", lwd = 3)

vers <- which.max(POaleatorisA)
aleatoriesA[, , vers]
POaleatorisA[vers]
modelsaleatorisA[[vers]]




# model restringint la transició de 3 a 2 com 0

Arestringida <- matrix(c(1/2, 1/4, 1/4, 1/4, 1/2, 0, 1/4, 1/4, 3/4), ncol = 3, nrow = 3)

modelArestringida <- estimarmodel(ON, Niters, estats, observacions, list(pi0, Arestringida, B0))

versemblança <- sum(exp(forward(modelArestringida, ON)[, Npasses]))



set.seed(123)

cadobs2 <- simHMM(hmmreal, Npasses)

passat2 <- cadobs2$states
ON2 <- cadobs2$observation


modelestimatalter <- estimarmodel(ON2, Niters, estats, observacions, list(pi0, A0, B0))

versalter <- sum(exp(forward(modelestimatalter, ON2)[, Npasses]))



