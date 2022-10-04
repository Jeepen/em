rm(list=ls())
set.seed(13072020)
library(rstudioapi)
library(data.table)
library(haven)
library(tidyverse)
library(survival)
library(doParallel)
cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)
setwd(dirname(getSourceEditorContext()$path))
source("functions.R")

# Parametre -------------------------------------------------------------------------------------------------------
n <- 1e3                                                  # Antal obs
lambda <- 1e-4                                            # Konstant baseline hazard 
HR <- 2                                                   # HR af behandling
HRsex <- 2
tau <- rexp(n, rate=1/1000)                               # Censorering
## tau <- 1000
lambdaPois <- 1                                           # Parameter til poissonfordeling
D <- 100                                                  # Behandlingslængde
M <- 4

## Simulation
starttime <- Sys.time()
ests <- foreach(i = 1:50, .combine = "rbind", #.options.RNG = 27102020, 
                .packages = c("data.table", "survival")) %dorng% {
  cat("outer iteration: ", i, "\n")
  # W <- rpois(n, lambda = lambdaPois) + 1                         # Antal behandlinger
  sex <- rbinom(n, 1, .5)                  
  W <- rpois(n, lambda = lambdaPois) + 1
  W[sex == 1] <- rpois(sum(sex == 1), lambda = 3) + 1
  # W <- sample(1:4, n, replace = TRUE)
  UW <- (W == 1) + 1.5 * (W == 2) + 4 * (W == 3) + 6 * (W >= 4)     # HR for forskellige antal behandlinger
  # UW <- pmin(exp(W - 4), 1)
  u <- runif(n)
  br <- 1 - exp(-lambda * W * D * UW * HR * HRsex^sex)
  X <- numeric(n)
  X[u < br] <- -log(1 - u[u < br]) / (lambda * HR * UW[u < br] * HRsex^sex[u < br])
  X[u >= br] <- (-log(1 - u[u >= br]) + lambda * UW[u >= br] * HRsex^sex[u >= br] * W[u >= br] * D * (1 - HR)) /
    (lambda * UW[u >= br] * HRsex^sex[u >= br])
  tau <- runif(n, 0, 1000)                                       # Censorering
  T <- pmin(X, tau)
  status <- as.numeric(X <= tau)
  obsW <- pmin(W, ceiling(T / D))                                # Vi observerer ikke W for alle
  treat <- as.numeric(T <= (D * W))                              # Behandling eller ej ved event tidspunkt
  statusW <- 1-treat                                             # Fortæller os om vi observerer W
  statusW[obsW >= M] <- 1                                        # Vi kender niveau hvis de behandles 4 eller flere gange
  tmp <- do.call("rbind", lapply(1:n, function(i){
    if(treat[i] == 1) data.table(id = i, time1 = 0, time2 = T[i], status = status[i], treat = 0, W = min(obsW[i],M),
                                 statusW = statusW[i], sex = sex[i])
    else data.table(id = i, time1 = c(0, D*W[i]), time2 = c(D*W[i], T[i]), status = c(0, status[i]),
                    treat = c(0,1), W = min(obsW[i],M), statusW = statusW[i], sex = sex[i])
  }))
  cox1 <- coxph(Surv(time1,time2,status) ~ treat + sex, data = tmp)
  cox2 <- coxph(Surv(time1,time2,status) ~ treat + strata(factor(W)) + sex, data = tmp)
  model <- em(Surv(time1,time2,status) ~ treat + strata(factor(W)) + sex, data = tmp, M = M)
  ests <- coef(model$model)[1]
  SE <- emSE(model)[1]
  cat("estimate", i, ": ", ests, "\n")
  cat("standard error", i, ": ", SE, "\n")
  # nn <- floor((T-1e-10) / D)
  # tmp2 <- do.call("rbind", lapply(1:n, function(i){
  #   sek <- seq(0, D * nn[i], D)
  #   data.table(id = i, time1 = sek, time2 = c(sek[-1], T[i]), status = c(rep(0, length(sek)-1), status[i]),
  #              sex = sex[i], treat = c(rep(0, obsW[i]), rep(1, length(sek) - obsW[i])), W = c(1:obsW[i], rep(obsW[i], length(sek) - obsW[i])))
  # }))
  # tmp2$W[tmp2$W >= M] <- M
  # cox3 <- coxph(Surv(time1,time2,status) ~ treat + sex + strata(factor(W)), data = tmp2)
  # cat("newCox", i, ": ", coef(cox3)[1], "\n")
  # c(ests, SE, coef(cox1)[1], coef(cox2)[1], coef(cox3)[1])
  c(ests, SE, coef(cox1)[1], coef(cox2)[1])
  ## c(i, ests)
}
endtime <- Sys.time()
difftime(endtime,starttime)



mean(-ests[,1])
exp(mean(-ests[,1]))
mean(ests[,2])
sd(ests[,1])
mean(ests[,ncol(ests)])
exp(mean(ests[,ncol(ests)]))
mean(-ests[,1]-1.96*ests[,2]<log(2) & -ests[,1]+1.96*ests[,2]>log(2))
sum(is.nan(ests))
apply(ests,2,mean)

## saveRDS(ests, "~/Dropbox/phd/Project2/R/RESULTSbigsim220615randCens.rds")