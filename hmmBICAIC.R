

library(HMM)

set.seed(1)

estimarmodel <- function(OBS, NiterEM, estatsmodel, observacionsmodel, paramsinicials) {
  hmm0 <- initHMM(estatsmodel, observacionsmodel, startProbs = paramsinicials[[1]], transProbs = paramsinicials[[2]], emissionProbs = paramsinicials[[3]])
  
  est <- baumWelch(hmm0, OBS, maxIterations = NiterEM, delta = 1E-6)
  
  est_inici <- est$hmm$startProbs
  est_trans <- est$hmm$transProbs
  est_emission <- est$hmm$emissionProbs
  
  return(initHMM(estatsmodel, observacionsmodel, startProbs = est_inici, transProbs = est_trans, emissionProbs = est_emission))
}

Npasses <- 1000 # quantitat d'observacions/estats

Niters <- 300 # quantitat d'iteracions de EM

Nestats <- 2:5

observacions <- c(1, 2) # cara, creu
estatsreals <- c(1, 2, 3)

Npars <- Nestats*(Nestats - 1) + Nestats*(length(observacions) - 1)

trans <- matrix(c(0.7, 0.15, 0.1, 0.2, 0.75, 0, 0.1, 0.1, 0.9), nrow = 3) # A
probas <- matrix(c(0.9, 0.5, 0.1, 0.1, 0.5, 0.9), ncol = 2) # B
inici <- c(0.2, 0.3, 0.5) # pi

hmmreal <- initHMM(estatsreals, observacions, startProbs = inici, transProbs = trans, emissionProbs = probas)

cadobs <- simHMM(hmmreal, Npasses)

ON <- cadobs$observation

passat <- cadobs$states

models <- vector(mode ="list", length = length(Nestats))
loglike <- numeric(length(Nestats))

for (i in seq_along(Nestats)) {
  k <- Nestats[i]
  estats <- as.vector(1:k)
  inici0 <- c(1, rep(0, k - 1))
  trans0 <- matrix(rep(1/k, k^2), nrow = k, ncol = k) 
  
  cosa1 <- k:1/(k + 1)
  
  probas0 <- matrix(cbind(cosa1, 1 - cosa1), ncol = 2)
  
  models[[i]] <- estimarmodel(ON, Niters, estats, observacions, list(inici0, trans0, probas0))
  loglike[[i]] <- log(sum(exp(forward(models[[i]], ON)[, Npasses])))
}

AIC <- 2*Npars - 2*loglike
BIC <- log(Npasses)*Npars - 2*loglike



# També ho podem veure en el cas d'un sol estat:

loglike1 <- sum(ON == 1)*log(sum(ON == 1)/length(ON)) + sum(ON == 2)*log(sum(ON == 2)/length(ON))
AIC1 <- 2 - 2*loglike1
BIC1 <- log(Npasses) - 2*loglike1











