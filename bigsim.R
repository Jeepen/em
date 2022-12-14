rm(list = ls())
set.seed(13072020)
# library(rstudioapi)
library(data.table)
library(mets)
# library(haven)
library(tidyverse)
library(survival)
library(doRNG)
library(doParallel)
cl <- parallel::makeCluster(40)
doParallel::registerDoParallel(cl)
# setwd(dirname(getSourceEditorContext()$path))
# setwd("ucph/hdir/SundKonsolidering_BioStatHome/Documents/")
source("functions.R")

# Parametre -------------------------------------------------------------------------------------------------------
n <- 2000                                                 # Antal obs
lambda <- 1e-4                                            # Konstant baseline hazard
HR <- 1.5                                                 # HR af behandling
HRsex <- 1
tau <- 1000
D <- 100                                                  # Behandlingslængde
M <- 4
beta <- -log(HR)

## Simulation
starttime <- Sys.time()
ests <- foreach(i = 1:10000, .combine = "rbind", .options.RNG = 13102022,  # 05102022
                .packages = c("data.table", "survival")) %dorng% {
                  cat("outer iteration: ", i, "\n")
                  sex <- rbinom(n, 1, .5)
                  W <- rpois(n, lambda = 1.5) + 1
                  W[sex == 1] <- rpois(sum(sex), lambda = 3.5) + 1
                  UW <- (W == 1) + (W == 2) * 2 + (W == 3) * 3 + (W >= 4) * 4
                  u <- runif(n)
                  br <- 1 - exp(-lambda * W * D * UW * HR * HRsex ^ sex)
                  X <- numeric(n)
                  X[u < br] <- -log(1 - u[u < br]) / (lambda * HR * UW[u < br] * HRsex ^ sex[u < br])
                  X[u >= br] <- (-log(1 - u[u >= br]) + lambda * UW[u >= br] * HRsex ^ sex[u >= br] * W[u >= br] * D * (1 - HR)) /
                    (lambda * UW[u >= br] * HRsex ^ sex[u >= br])
                  T <- pmin(X, tau)
                  status <- as.numeric(X <= tau)
                  obsW <- pmin(W, ceiling(T / D))                                # Vi observerer ikke W for alle
                  treat <- as.numeric(T <= (D * W))                              # Behandling eller ej ved event tidspunkt
                  statusW <- 1 - treat                                           # Fortæller os om vi observerer W
                  statusW[obsW >= M] <- 1                                        # Vi kender niveau hvis de behandles 4 eller flere gange
                  tmp <- do.call("rbind", lapply(1:n, function(i) {
                    if (treat[i] == 1){
                      data.table(id = i, time1 = 0, time2 = T[i], status = status[i], treat = 0,
                                 W = min(obsW[i], M), statusW = statusW[i], sex = sex[i])
                    } 
                    else{
                      data.table(id = i, time1 = c(0, D * W[i]), time2 = c(D * W[i], T[i]), status = c(0, status[i]),
                                 treat = c(0, 1), W = min(obsW[i], M), statusW = statusW[i], sex = sex[i]) 
                    }
                  }))
                  cox1 <- coxph(Surv(time1, time2, status) ~ treat + sex, data = tmp)
                  cc <- confint(cox1)
                  cox1_cov <- cc[1, 1] < beta & cc[1, 2] > beta
                  cox1_se <- sqrt(cox1$var[1, 1])
                  cox2 <- coxph(Surv(time1, time2, status) ~ treat + factor(W) + sex, data = tmp)
                  cc <- confint(cox2)
                  cox2_cov <- cc[1, 1] < beta & cc[1, 2] > beta
                  cox2_se <- sqrt(cox2$var[1, 1])
                  model <- em(Surv(time1, time2, status) ~ treat + factor(W) + sex, data = tmp, M = M, sex = TRUE)
                  ests1 <- coef(model$model)[1]
                  SE1 <- emSE(model, sex = TRUE)[1]
                  model1_cov <- (ests1 - 1.96 * SE1 < beta) & (ests1 + 1.96 * SE1 > beta)
                  cat("estimate", i, ": ", ests1, "\n")
                  cat("standard error", i, ": ", SE1, "\n")
                  model <- em(Surv(time1, time2, status) ~ treat + factor(W), data = tmp, M = M)
                  ests2 <- coef(model$model)[1]
                  SE2 <- emSE(model, M = M)[1]
                  model2_cov <- (ests2 - 1.96 * SE2 < beta) & (ests2 + 1.96 * SE2 > beta)
                  cat("estimate", i, ": ", ests2, "\n")
                  cat("standard error", i, ": ", SE2, "\n")
                  
                  nn <- floor((T - 1e-10) / D)
                  tmp2 <- do.call("rbind", lapply(1:n, function(i) {
                    sek <- seq(0, D * nn[i], D)
                    data.table(id = i, time1 = sek, time2 = c(sek[-1], T[i]), status = c(rep(0, length(sek) - 1), status[i]),
                               sex = sex[i], treat = c(rep(0, obsW[i]), rep(1, length(sek) - obsW[i])),
                               W = c(1:obsW[i], rep(obsW[i], length(sek) - obsW[i])))
                  }))
                  tmp2$W[tmp2$W >= M] <- M
                  cox3 <- coxph(Surv(time1, time2, status) ~ treat + sex + factor(W), data = tmp2)
                  cc <- confint(cox3)
                  cox3_cov <- cc[1, 1] < beta & cc[1, 2] > beta
                  cox3_se <- sqrt(cox3$var[1, 1])
                  cat("newCox", i, ": ", coef(cox3)[1], "\n")
                  c(ests1, SE1, model1_cov, ests2, SE2, model2_cov, coef(cox1)[1], cox1_se, cox1_cov,
                    coef(cox2)[1], cox2_se, cox2_cov, coef(cox3)[1], cox3_se, cox3_cov)
                }
endtime <- Sys.time()
difftime(endtime, starttime)
apply(ests, 2, mean)

mean(-ests[, 1])
exp(mean(-ests[, 1]))
mean(ests[, 2])
sd(ests[, 1])
mean(ests[, ncol(ests)])
exp(mean(ests[, ncol(ests)]))
mean(-ests[, 1] - 1.96 * ests[, 2] < log(1.5) & -ests[, 1] + 1.96 * ests[, 2] > log(1.5))
sum(is.nan(ests))
apply(ests, 2, mean)

# saveRDS(list(ests = ests, time = difftime(endtime, starttime)), "simresults_null.rds")
