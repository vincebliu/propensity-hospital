design <- function(x, t) {
  n <- length(t)+1
  t <- c(0, t, length(x))
  X <- NULL
  for (i in 1:n) {
    xi <- cbind(1, x[(t[i]+1):t[i+1]])
    if (i > 1) {
      for (j in 1:(i-1)) {
        xi <- cbind(0, 0, xi)
      }
    }
    if (i < n) {
      for (j in 1:(n-i)) {
        xi <- cbind(xi, 0, 0)
      }
    }
    X <- rbind(X, xi)
  }
  return(X)
  
}

Logpost <- function(X, Y, beta) {
  return (sum(-exp(X%*%beta)+(X%*%beta * Y)))
}

mode <- function(x) {
  return(as.numeric(names(sort(table(x), decreasing=TRUE))[1]))
}

dmvnlog <- function(x, mu, sig) {
  - 1/2 * t(x-mu) %*% solve(sig) %*% (x-mu)
}

DIC <- function(x, Y, params, CHGPT) {
  P <- 0
  for (i in 1:nrow(params)) {
    P <- P + -2 * Logpost(design(x, params[i,((CHGPT+1)*2+1):ncol(params)]), Y, params[i,1:((CHGPT+1)*2)])
  }
  modes <- NULL
  for (i in 1:CHGPT) {
    modes[i] <- mode(params[,(CHGPT+1)*2+i])
  }
  G <- -2* Logpost(design(x, modes), Y, colMeans(params[,1:((CHGPT+1)*2)]))
  return (2 / nrow(params) * P - G)
}
mcmc <- function(x, Y, M, beta, t) {
  ac <- rep(0, ncol(beta) + length(t))
  out <- matrix(0, M, ncol(beta)*2 + length(t))
  for (iter in 1:M) {
    #update beta
    X <- design(x, t)
    for (i in seq_len(ncol(beta))) {
      proposal <- mvrnorm(1, beta[,i], tuning)
      prime <- beta
      prime[,i] <- proposal
      A <- Logpost(X, Y, matrix(prime)) - Logpost(X, Y, matrix(beta)) + dmvnlog(proposal, beta0, Sig) - dmvnlog(beta[,i], beta0, Sig) + dmvnlog(beta[,i], proposal, tuning) - dmvnlog(proposal, beta[,i], tuning)
      if (log(runif(1)) < A) {
        beta <- prime
        ac[i] <- ac[i]+1
      }
    }
    
    #update sig
    S <- matrix(0, 2, 2)
    for (i in 1:ncol(beta)) {
      S <- S + (beta[,i] - beta0) %*% t(beta[,i] - beta0)
    }
    Sig <- solve(rWishart(1, 6, solve(diag(0, 2) + S))[,,1])
    
    #update beta0
    V <- solve(solve(Sig + iC))
    beta0 <- mvrnorm(1, V%*%(solve(Sig) %*% rowSums(beta) + iC %*% mu0), V)
    
    if (iter %% 1000 == 0) {
      print(list(Sig=Sig, 
                 beta0=beta0))
    }
    
    #update t
    for (i in 1:length(t)) {
      low <- t[i-1]+2
      if (i==1) {
        low <- 2
      }
      high <- t[i+1]-2
      if (i==length(t)) {
        high <- length(x)-2
      }
      if (low != high) {
        p <- matrix(0, high-low+1, high-low+1)
        for (j in 1:(high-low+1)) {
          p[j,] <- exp(-.15 * abs(low:high - (j+low-1)))
          p[j,]<-p[j,]/sum(p[j,])
        }
        proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
        prime <- t
        prime[i] <- proposed
        A <- Logpost(design(x, prime), Y, matrix(beta)) - Logpost(design(x, t), Y, matrix(beta)) + log(p[proposed-low+1,t[i]-low+1]) - log(p[t[i]-low+1, proposed-low+1])
        if (log(runif(1)) < A) {
          ac[i+ncol(beta)] <- ac[i+ncol(beta)]+1
          t[i] <- proposed
        }
      }
    }
    out[iter,]<-c(unlist(beta), t)
  }
  print(ac)
  return(out)
}

# Toy Example
Y <- NULL
for (i in 1:50) {
  Y[i] <- rpois(1, exp(1 + 0.1 * i))
}
for (i in 51:100) {
  Y[i] <- rpois(1, exp(7 - 0.02 * i))
}


tuning <- matrix(c(.1, 0, 0, .001), 2)
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
toy <- mcmc(seq_along(Y), Y, 20000, matrix(0, 2, 2), 80)
save(toy, file="toy.RData")
load("toy.RData")
toya <- mcmc(seq_along(Y), Y, 20000, matrix(c(-1, .1, 1, -.1), 2), 40)
toyb <- mcmc(seq_along(Y), Y, 20000, matrix(c(-2, .2, 2, -.2), 2), 60)

par(mfrow=c(2,2))
for (i in c(1,2,3,5)) {
  plot(toy[,i], type="l", ylim=range(c(toy[,i], toya[,i], toyb[,i])), xlab="", ylab="")
  lines(toya[,i], col="green")
  lines(toyb[,i], col="pink")
}

par(mfrow=c(1,1))
plot(Y, xlab="", ylab="")
t <- mode(toy[,5])
b <- colMeans(toy[,1:4])
lines(1:t, exp(b[1] + b[2] * (1:t)), lty=2, lwd=2, col="cyan")
lines((t+1):length(Y), exp(b[3] + b[4] * ((t+1):length(Y))), lty=2, lwd=2, col="cyan")

data <- read.csv("owid-covid-data.csv")
us <- data[data$iso_code=="USA",]
us <- us[, c("date", "new_cases")]
us <- us[complete.cases(us),]
us$date <- as.Date(us$date)

# Covid: single Change point
Y <- us$new_cases[250:425]
tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
Cov1 <- mcmc(seq_along(Y), Y, 100000, matrix(0, 2, 2), 50)
save(Cov1, file="Cov1.RData")
load("Cov1.RData")

Cov1a <- mcmc(seq_along(Y), Y, 100000, matrix(c(2, 0.2, 5, -0.2), 2), 100)
Cov1b <- mcmc(seq_along(Y), Y, 100000, matrix(c(3, 0.1, 7, -0.3), 2), 140)
save(Cov1a, file="Cov1a.RData")
save(Cov1b, file="Cov1b.RData")

par(mfrow=c(2,2))
for (i in c(1, 3, 4, 5)) {
  plot(Cov1[,i], type="l", ylim=range(c(Cov1[,i], Cov1a[,i], Cov1b[,i])), xlab="", ylab="")
  lines(Cov1a[,i], col="green")
  lines(Cov1b[,i], col="pink")
}

par(mfrow=c(1,1))
plot(us$date[250:425], Y, xlab="", ylab="")
t <- mode(Cov1[50000:100000,5])
b <- colMeans(Cov1[50000:100000,1:4])
lines(us$date[250:425][1:t], exp(b[1] + b[2] * (1:t)), lty=2, lwd=2, col="purple")
lines(us$date[250:425][(t+1):length(Y)], exp(b[3] + b[4] * ((t+1):length(Y))), lty=2, lwd=2, col="purple")

#Covid 2 change points
Y <- us$new_cases[250:600]
tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
Cov2 <- mcmc(seq_along(Y), Y, 100000, matrix(0, 2, 3), c(50, 150))
save(Cov2, file="Cov2.RData")
load("Cov2.RData")

plot(Y)
t1 <- mode(Cov2[50000:100000,7])
t2 <- mode(Cov2[50000:100000,8])
b <- colMeans(Cov2[50000:100000,1:6])
lines(1:t1, exp(b[1] + b[2] * (1:t1)))
lines((t1+1):t2, exp(b[3] + b[4] * ((t1+1):t2)))
lines((t2+1):length(Y), exp(b[5] + b[6] * ((t2+1):length(Y))))
 
Cov2a <- mcmc(seq_along(Y), Y, 100000, matrix(c(1,.1,2,-.1,1,.1), 2, 3), c(200,300))
save(Cov2a, file="Cov2a.RData")
plot(Y)
t1 <- mode(Cov2a[50000:100000,7])
t2 <- mode(Cov2a[50000:100000,8])
b <- colMeans(Cov2a[50000:100000,1:6])
lines(1:t1, exp(b[1] + b[2] * (1:t1)))
lines((t1+1):t2, exp(b[3] + b[4] * ((t1+1):t2)))
lines((t2+1):length(Y), exp(b[5] + b[6] * ((t2+1):length(Y))))

Cov2b <- mcmc(seq_along(Y), Y, 100000, matrix(c(-1,.2,3,-.2,-1,.2), 2, 3), c(150,250))
save(Cov2b, file="Cov2b.RData")
plot(Y)
t1 <- mode(Cov2b[50000:100000,7])
t2 <- mode(Cov2b[50000:100000,8])
b <- colMeans(Cov2b[50000:100000,1:6])
lines(1:t1, exp(b[1] + b[2] * (1:t1)))
lines((t1+1):t2, exp(b[3] + b[4] * ((t1+1):t2)))
lines((t2+1):length(Y), exp(b[5] + b[6] * ((t2+1):length(Y))))

DIC(seq_along(Y), Y, Cov2[50000:10000,],2 )
DIC(seq_along(Y), Y, Cov2a[50000:10000,],2 )
DIC(seq_along(Y), Y, Cov2b[50000:10000,],2 )

load("Cov2a.RData"); load("Cov2b.RData")

par(mfrow=c(2,2))
for (i in 5:8) {
  plot(Cov2[,i], type="l", ylim=range(c(Cov2[,i], Cov2a[,i], Cov2b[,i])), xlab="", ylab="")
  lines(Cov2a[,i], col="green")
  lines(Cov2b[,i], col='pink')
}

# Four Change points
mcmc <- function(x, Y, M, beta, t) {
  ac <- rep(0, ncol(beta) + length(t))
  out <- matrix(0, M, ncol(beta)*2 + length(t))
  for (iter in 1:M) {
    #update beta
    X <- design(x, t)
    for (i in seq_len(ncol(beta))) {
      proposal <- mvrnorm(1, beta[,i], tuning)
      if (i==1) {
        while (proposal[1] < 5 | proposal[1]>15) {
          proposal <- mvrnorm(1, beta[,i], tuning)
        }
      } else if (i==2) {
        while (proposal[2] > 0 | proposal[1] < 9 | proposal[1]>19) {
          proposal <- mvrnorm(1, beta[,i], tuning)
        }
      } else if (i==4) {
        while (proposal[2] >0 | proposal[1] < 12 | proposal[1] > 22) {
          proposal <- mvrnorm(1, beta[,i], tuning)
        }
      }
      prime <- beta
      prime[,i] <- proposal
      A <- Logpost(X, Y, matrix(prime)) - Logpost(X, Y, matrix(beta)) + dmvnlog(proposal, beta0, Sig) - dmvnlog(beta[,i], beta0, Sig) + dmvnlog(beta[,i], proposal, tuning) - dmvnlog(proposal, beta[,i], tuning)
      if (log(runif(1)) < A) {
        beta <- prime
        ac[i] <- ac[i]+1
      }
    }
    
    #update sig
    S <- matrix(0, 2, 2)
    for (i in 1:ncol(beta)) {
      S <- S + (beta[,i] - beta0) %*% t(beta[,i] - beta0)
    }
    Sig <- solve(rWishart(1, 6, solve(diag(0, 2) + S))[,,1])
    
    #update beta0
    V <- solve(solve(Sig + iC))
    beta0 <- mvrnorm(1, V%*%(solve(Sig) %*% rowSums(beta) + iC %*% mu0), V)
    
    if (iter %% 1000 == 0) {
      print(list(Sig=Sig, 
                 beta0=beta0))
    }
    
    #update t
    for (i in 1:length(t)) {
      low <- t[i-1]+2
      if (i==1) {
        low <- 2
      }
      high <- t[i+1]-2
      if (i==length(t)) {
        high <- length(x)-2
      }
      if (low != high) {
        p <- matrix(0, high-low+1, high-low+1)
        for (j in 1:(high-low+1)) {
          p[j,] <- exp(-.15 * abs(low:high - (j+low-1)))
          p[j,]<-p[j,]/sum(p[j,])
        }
        proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
        if (i == 2) {
          while (proposed < 150-10 | proposed > 150+10) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==3) {
          while ( proposed < 200-10 | proposed > 200+10) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        }
        prime <- t
        prime[i] <- proposed
        A <- Logpost(design(x, prime), Y, matrix(beta)) - Logpost(design(x, t), Y, matrix(beta)) + log(p[proposed-low+1,t[i]-low+1]) - log(p[t[i]-low+1, proposed-low+1])
        if (log(runif(1)) < A) {
          ac[i+ncol(beta)] <- ac[i+ncol(beta)]+1
          t[i] <- proposed
        }
      }
    }
    out[iter,]<-c(unlist(beta), t)
  }
  print(ac)
  return(out)
}

tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 2)
beta0 <- rowSums(beta0)
mu0 <- beta0
iC <- diag(0, 2)
Cov4 <- mcmc(seq_along(Y), Y, 100000, matrix(c(10, 0, 14, 0, 0, 0, 17, 0, 0, 0), 2), c(50, 150, 200, 250))
save(Cov4, file="Cov4.RData")
load("Cov4.RData")

Cov4a <- mcmc(seq_along(Y), Y, 100000, matrix(c(5, .1, 9, -.1, 1, .1, 12, -.1, 1, .1), 2), c(100, 140, 190, 300))
save(Cov4a, file="Cov4a.RData")

Cov4b <- mcmc(seq_along(Y), Y, 100000, matrix(c(15, .05, 19, -.2, 2, .05, 22, -.2, 2, .05), 2), c(75, 160, 210, 275))
save(Cov4b, file="Cov4b.RData")
load("Cov4b.RData")

DIC(seq_along(Y), Y, Cov4[50000:100000,], 4)
DIC(seq_along(Y), Y, Cov4a[50000:100000,], 4)
DIC(seq_along(Y), Y, Cov4b[50000:100000,], 4)

par(mfrow=c(1,2))
plot(us$date[250:600], Y, xlab="", ylab="", pch=16, cex=.5, col="grey", main="(a)")
t1 <- mode(Cov2a[50000:100000,7])
t2 <- mode(Cov2a[50000:100000,8])
b <- colMeans(Cov2a[50000:100000,1:6])
lines(us$date[250:600][1:t1], exp(b[1] + b[2] * (1:t1)), lty=2, lwd=2, col="brown")
lines(us$date[250:600][(t1+1):t2], exp(b[3] + b[4] * ((t1+1):t2)), lty=2, lwd=2, col="brown")
lines(us$date[250:600][(t2+1):length(Y)], exp(b[5] + b[6] * ((t2+1):length(Y))), lty=2, lwd=2, col="brown")
plot(us$date[250:600], Y, xlab="", ylab="", pch=16, cex=.5, main="(b)", col="grey")
t1 <- mode(Cov4b[50000:100000,11])
t2 <- mode(Cov4b[50000:100000,12])
t3 <- mode(Cov4b[50000:100000,13])
t4 <- mode(Cov4b[50000:100000,14])
b <- colMeans(Cov4b[50000:100000,1:10])
lines(us$date[250:600][1:t1], exp(b[1] + b[2] * (1:t1)), lty=2, lwd=2, col="orange")
lines(us$date[250:600][(t1+1):t2], exp(b[3] + b[4] * ((t1+1):t2)), lty=2, lwd=2, col="orange")
lines(us$date[250:600][(t2+1):t3], exp(b[5] + b[6] * ((t2+1):t3)), lty=2, lwd=2, col="orange")
lines(us$date[250:600][(t3+1):t4], exp(b[7] + b[8] * ((t3+1):t4)), lty=2, lwd=2, col="orange")
lines(us$date[250:600][(t4+1):length(Y)], exp(b[9] + b[10] * ((t4+1):length(Y))), lty=2, lwd=2, col="orange")

#GBR
uk <- data[data$iso_code=="GBR",]
uk <- uk[, c("date", "new_cases")]
uk <- uk[complete.cases(uk),]
uk$date <- as.Date(uk$date)

#Three change points
mcmc <- function(x, Y, M, beta, t) {
  ac <- rep(0, ncol(beta) + length(t))
  out <- matrix(0, M, ncol(beta)*2 + length(t))
  for (iter in 1:M) {
    #update beta
    X <- design(x, t)
    for (i in seq_len(ncol(beta))) {
      proposal <- mvrnorm(1, beta[,i], tuning)
      if (i==2 | i==4) {
        while (proposal[2] > 0) {
          proposal <- mvrnorm(1, beta[,i], tuning)
        }
      }
      prime <- beta
      prime[,i] <- proposal
      A <- Logpost(X, Y, matrix(prime)) - Logpost(X, Y, matrix(beta)) + dmvnlog(proposal, beta0, Sig) - dmvnlog(beta[,i], beta0, Sig) + dmvnlog(beta[,i], proposal, tuning) - dmvnlog(proposal, beta[,i], tuning)
      if (log(runif(1)) < A) {
        beta <- prime
        ac[i] <- ac[i]+1
      }
    }
    
    #update sig
    S <- matrix(0, 2, 2)
    for (i in 1:ncol(beta)) {
      S <- S + (beta[,i] - beta0) %*% t(beta[,i] - beta0)
    }
    Sig <- solve(rWishart(1, 6, solve(diag(0, 2) + S))[,,1])
    
    #update beta0
    V <- solve(solve(Sig + iC))
    beta0 <- mvrnorm(1, V%*%(solve(Sig) %*% rowSums(beta) + iC %*% mu0), V)
    
    if (iter %% 1000 == 0) {
      print(list(Sig=Sig, 
                 beta0=beta0))
    }
    
    #update t
    for (i in 1:length(t)) {
      low <- t[i-1]+2
      if (i==1) {
        low <- 2
      }
      high <- t[i+1]-2
      if (i==length(t)) {
        high <- length(x)-2
      }
      if (low != high) {
        p <- matrix(0, high-low+1, high-low+1)
        for (j in 1:(high-low+1)) {
          p[j,] <- exp(-.15 * abs(low:high - (j+low-1)))
          p[j,]<-p[j,]/sum(p[j,])
        }
        proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
        if (i == 1) {
          while (proposed < 50-10 | proposed > 50+10) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==3) {
          while ( proposed < 200 | proposed > 250) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        }
        prime <- t
        prime[i] <- proposed
        A <- Logpost(design(x, prime), Y, matrix(beta)) - Logpost(design(x, t), Y, matrix(beta)) + log(p[proposed-low+1,t[i]-low+1]) - log(p[t[i]-low+1, proposed-low+1])
        if (log(runif(1)) < A) {
          ac[i+ncol(beta)] <- ac[i+ncol(beta)]+1
          t[i] <- proposed
        }
      }
    }
    out[iter,]<-c(unlist(beta), t)
  }
  print(ac)
  return(out)
}

Y <- uk$new_cases[300:550]
Y[which(Y<0)] <- 0

coef <- matrix(0, 2, 4)
W <- 1:50
coef[,1] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 51:150
coef[,2] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 151:230
coef[,3] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 231:length(Y)
coef[,4] <- glm(Y[W] ~ W, family=poisson())$coef

tuning <- matrix(c(3, 0, 0, .05), 2) / 10000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
Cov3 <- mcmc(seq_along(Y), Y, 100000, coef, c(50, 150, 200))

plot(uk$date[300:550], Y, xlab="", ylab="", main="UK", pch=5)
t1 <- mode(Cov3[50000:100000,9])
t2 <- mode(Cov3[50000:100000,10])
t3 <- mode(Cov3[50000:100000,11])
b <- colMeans(Cov3[50000:100000,1:8])
lines(uk$date[300:550][1:t1], exp(b[1] + b[2] * (1:t1)), lwd=2, lty=2, col="pink")
lines(uk$date[300:550][(t1+1):t2], exp(b[3] + b[4] * ((t1+1):t2)), lwd=2, lty=2, col="pink")
lines(uk$date[300:550][(t2+1):t3], exp(b[5] + b[6] * ((t2+1):t3)), lwd=2, lty=2, col="pink")
lines(uk$date[300:550][(t3+1):length(Y)], exp(b[7] + b[8] * ((t3+1):length(Y))), lwd=2, lty=2, col="pink")

save(Cov3, file="Cov3.RData")
load("Cov3.RData")

#India
ind <- data[data$iso_code=="IND",]
ind <- ind[, c("date", "new_cases")]
ind <- ind[complete.cases(ind),]
ind$date <- as.Date(ind$date)

Y <- ind$new_cases[1:590]
coef <- matrix(0, 2, 6)
W <- 1:150
coef[,1] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 151:240
coef[,2] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 241:375
coef[,3] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 376:475
coef[,4] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 476:510
coef[,5] <- glm(Y[W] ~ W, family=poisson())$coef
W <- 511:length(Y)
coef[,6] <- glm(Y[W] ~ W, family=poisson())$coef

mcmc <- function(x, Y, M, beta, t) {
  ac <- rep(0, ncol(beta) + length(t))
  out <- matrix(0, M, ncol(beta)*2 + length(t))
  for (iter in 1:M) {
    #update beta
    X <- design(x, t)
    for (i in seq_len(ncol(beta))) {
      proposal <- mvrnorm(1, beta[,i], tuning)
      if (i==3 | i==5) {
        while (proposal[2] > 0) {
          proposal <- mvrnorm(1, beta[,i], tuning)
        }
      }
      prime <- beta
      prime[,i] <- proposal
      A <- Logpost(X, Y, matrix(prime)) - Logpost(X, Y, matrix(beta)) + dmvnlog(proposal, beta0, Sig) - dmvnlog(beta[,i], beta0, Sig) + dmvnlog(beta[,i], proposal, tuning) - dmvnlog(proposal, beta[,i], tuning)
      if (log(runif(1)) < A) {
        beta <- prime
        ac[i] <- ac[i]+1
      }
    }
    
    #update sig
    S <- matrix(0, 2, 2)
    for (i in 1:ncol(beta)) {
      S <- S + (beta[,i] - beta0) %*% t(beta[,i] - beta0)
    }
    Sig <- solve(rWishart(1, 6, solve(diag(0, 2) + S))[,,1])
    
    #update beta0
    V <- solve(solve(Sig + iC))
    beta0 <- mvrnorm(1, V%*%(solve(Sig) %*% rowSums(beta) + iC %*% mu0), V)
    
    if (iter %% 1000 == 0) {
      print(list(Sig=Sig, 
                 beta0=beta0,
                 iter=iter))
    }
    
    #update t
    for (i in 1:length(t)) {
      low <- t[i-1]+2
      if (i==1) {
        low <- 2
      }
      high <- t[i+1]-2
      if (i==length(t)) {
        high <- length(x)-2
      }
      if (low != high) {
        p <- matrix(0, high-low+1, high-low+1)
        for (j in 1:(high-low+1)) {
          p[j,] <- exp(-.15 * abs(low:high - (j+low-1)))
          p[j,]<-p[j,]/sum(p[j,])
        }
        proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
        if (i == 1) {
          while (proposed < 100 | proposed > 200) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==2) {
          while ( proposed < 200 | proposed > 250) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==3) {
          while ( proposed < 400-20 | proposed > 400+20) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==4) {
          while ( proposed < 475-20 | proposed > 475+20) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        } else if (i==5) {
          while ( proposed < 500-20 | proposed > 500+20) {
            proposed <- sample(low:high, 1, prob=p[t[i]-low+1,])
          }
        }
        prime <- t
        prime[i] <- proposed
        A <- Logpost(design(x, prime), Y, matrix(beta)) - Logpost(design(x, t), Y, matrix(beta)) + log(p[proposed-low+1,t[i]-low+1]) - log(p[t[i]-low+1, proposed-low+1])
        if (log(runif(1)) < A) {
          ac[i+ncol(beta)] <- ac[i+ncol(beta)]+1
          t[i] <- proposed
        }
      }
    }
    out[iter,]<-c(unlist(beta), t)
  }
  print(ac)
  return(out)
}

tuning <- matrix(c(3, 0, 0, .05), 2) / 10000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
Cov5 <- mcmc(seq_along(Y), Y, 30000, coef, c(150, 240, 375, 475, 500))

par(mfrow=c(1,1))
plot(ind$date[1:590], Y, xlab="", ylab="", main="India", pch=4)
t1 <- mode(Cov5[10000:30000,13])
t2 <- mode(Cov5[10000:30000,14])
t3 <- mode(Cov5[10000:30000,15])
t4 <- mode(Cov5[10000:30000,16])
t5 <- mode(Cov5[10000:30000,17])
b <- colMeans(Cov5[10000:30000,1:12])
lines(ind$date[1:590][1:t1], exp(b[1] + b[2] * (1:t1)), lwd=2, lty=2, col="blue")
lines(ind$date[1:590][(t1+1):t2], exp(b[3] + b[4] * ((t1+1):t2)), lwd=2, lty=2, col="blue")
lines(ind$date[1:590][(t2+1):t3], exp(b[5] + b[6] * ((t2+1):t3)), lwd=2, lty=2, col="blue")
lines(ind$date[1:590][(t3+1):t4], exp(b[7] + b[8] * ((t3+1):t4)), lwd=2, lty=2, col="blue")
lines(ind$date[1:590][(t4+1):t5], exp(b[9] + b[10] * ((t4+1):t5)), lwd=2, lty=2, col="blue")
lines(ind$date[1:590][(t5+1):length(Y)], exp(b[11] + b[12] * ((t5+1):length(Y))), lwd=2, lty=2, col="blue")

save(Cov5, file="Cov5.RData")

load("Cov5.RData")

## Omicron
data <- read.csv("owid-covid-data1.csv")
three <- data[data$iso_code %in% c("USA", "GBR", "ZAF"),]
three <- three[, c("iso_code", "date", "new_cases")]
three <- three[complete.cases(three),]
three$date <- as.Date(three$date)
us <- three[three$iso_code=="USA",]
uk <- three[three$iso_code=="GBR",]
sa <- three[three$iso_code=="ZAF",]

Y1 <- us$new_cases[600:701]
Y2 <- uk$new_cases[600:nrow(uk)]
Y3 <- sa$new_cases[600:nrow(sa)]

#use the first mcmc function

#us
tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
coef <- matrix(0, 2, 2)
W <- 1:50
coef[,1] <- glm(Y1[W] ~ W, family=poisson())$coef
W <- 51:length(Y1)
coef[,2] <- glm(Y1[W] ~ W, family=poisson())$coef
Cova <- mcmc(seq_along(Y1), Y1, 100000, coef, 50)

#uk
coef <- matrix(0, 2, 2)
W <- 1:80
coef[,1] <- glm(Y2[W] ~ W, family=poisson())$coef
W <- 81:length(Y2)
coef[,2] <- glm(Y2[W] ~ W, family=poisson())$coef
tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
Covb <- mcmc(seq_along(Y2), Y2, 100000, coef, 80)

#sa
tuning <- matrix(c(3, 0, 0, .05), 2) / 1000
Sig <- matrix(c(3, 0, 0, .1), 2)
beta0 <- c(0, 0)
mu0 <- beta0
iC <- diag(0, 2)
coef <- matrix(0, 2, 2)
W <- 1:30
coef[,1] <- glm(Y3[W] ~ W, family=poisson())$coef
W <- 31:length(Y3)
coef[,2] <- glm(Y3[W] ~ W, family=poisson())$coef
Covc <- mcmc(seq_along(Y3), Y3, 100000, coef, 30)

save(Cova, file="Cova.RData")
save(Covb, file="Covb.RData")
save(Covc, file="Covc.RData")

load("Cova.RData")
load("Covb.RData")
load("Covc.RData")

par(mfrow=c(1, 3))

plot(us$date[600:701], Y1, xlab="", ylab="", main="(a)")
t <- mode(Cova[50000:100000,5])
b <- colMeans(Cova[50000:100000,1:4])
lines(us$date[600:701][1:t], exp(b[1] + b[2] * (1:t)), lwd=2, lty=2, col="red")
lines(us$date[600:701][(t+1):length(Y1)], exp(b[3] + b[4] * ((t+1):length(Y1))), lwd=2, lty=2, col="red")

plot(uk$date[600:nrow(uk)], Y2, xlab="", ylab="", main="(b)")
t <- mode(Covb[50000:100000,5])
b <- colMeans(Covb[50000:100000,1:4])
lines(uk$date[600:nrow(uk)][1:t], exp(b[1] + b[2] * (1:t)), lwd=2, lty=2, col="blue")
lines(uk$date[600:nrow(uk)][(t+1):length(Y2)], exp(b[3] + b[4] * ((t+1):length(Y2))), lwd=2, lty=2, col="blue")

plot(sa$date[600:nrow(sa)], Y3, xlab="", ylab="", main="(c)")
t <- mode(Covc[50000:100000,5])
b <- colMeans(Covc[50000:100000,1:4])
lines(sa$date[600:nrow(sa)][1:t], exp(b[1] + b[2] * (1:t)), lwd=2, lty=2, col="green")
lines(sa$date[600:nrow(sa)][(t+1):length(Y3)], exp(b[3] + b[4] * ((t+1):length(Y3))), lwd=2, lty=2, col="green")
