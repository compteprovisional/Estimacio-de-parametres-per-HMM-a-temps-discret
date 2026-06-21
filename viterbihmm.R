library(HMM)

set.seed(1)

estats <- c(1, 2, 3)
observacions <- c(1, 2)
N <- 10

inici <- c(1, 0, 0) # pi
probas <- matrix(c(1, 0.6, 0.5, 0, 0.4, 0.5), ncol = 2) # B
probas

trans <- matrix(c(0, 0.2, 0.2, 0.6, 0.8, 0, 0.4, 0, 0.8), nrow = 3) # A
trans

ON <- numeric(N)
passat <- numeric(N)
estat_actual <- sample(estats, size = 1, prob = inici)

for (i in 1:N) {
  ON[i] <- sample(observacions, size = 1, prob = probas[estat_actual, ])
  passat[i] <- estat_actual
  estat_actual <- sample(estats, size = 1, prob = trans[estat_actual, ])
}

print(ON)
print(passat)


hmm <- initHMM(estats, observacions, startProbs = inici, transProbs = trans, emissionProbs = probas)

hmm

forwards <- exp(forward(hmm, ON))

backwards <- exp(backward(hmm, ON))

forwards
backwards

den <- numeric(N)

for (i in 1:N) {
  den[i] <- sum(forwards[, i]%*%backwards[, i])
}


gamma <- matrix(0, ncol = N, nrow = length(estats))

for (t in 1:N) {
  for (i in estats) {
    gamma[i, t] <- forwards[i, t]*backwards[i, t]/den[t]
    print(forwards[i, t]*backwards[i, t])
  }
}

gamma

Qestimada <- numeric(N)

for (i in 1:N) {
  Qestimada[i] <- which.max(gamma[, i])
}

Qestimada

Qbona <- viterbi(hmm, ON)

