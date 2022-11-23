## Estimation
em <- function(formula, data, M = 4, maxit = 1000, object = NULL, sex = FALSE){
  ## Nogle parametre til at afgøre konvergens
  eps1 <- 1e-10
  eps2 <- 0#1e-6
  n <- length(unique(data$id))                                        # Det er antaget, at der er en id variabel til at identificere individer
  if(is.null(object)){
    WW <- tapply(data$W, data$id, unique)                               # Det antages at antal behandlinger er kaldt W
    statusWW <- tapply(data$statusW, data$id, unique)                   # statusW fortæller os, hvorvidt vi har observeret faktisk antal behandlinger eller minimum
    probs <- -diff(c(1, survfit(Surv(WW, statusWW) ~ 1)$surv))          # Marginal fordeling af W skal estimeres ikke-parametrisk men bliver initialiseret med Kaplan-Meier 
    if(sex){
      probs <- -apply(rbind(1, matrix(survfit(Surv(W, statusW) ~ sex, data = data)$surv, ncol = 2, nrow = 4)),2,diff)
      probs <- cbind(probs[,2], probs[,1])
    }    
    model <- coxph(formula, data, model = TRUE, iter = 0)
    # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))      # Vi skal bruge baseline hazard til estimation af likelihood bidrag
  }
  else{
    probs <- object$probs
    model <- object$model
    # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))      # Vi skal bruge baseline hazard til estimation af likelihood bidrag
  }
  # if(!is.null(beta)){
  #   model <- coxph(formula, data, model = TRUE, iter = 0, init = beta)# Cox model fittes, hvor vi ignorerer at W ikke altid er det sande W men nogle gange blot er minimum
  # }
  # else{
  #   model <- coxph(formula, data, model = TRUE, iter = 0)             # Cox model fittes, hvor vi ignorerer at W ikke altid er det sande W men nogle gange blot er minimum
  # }
  crit <- TRUE                                                        # Et kriterie relateret til ændringen i den marginale fordeling af W
  bigdata <- rbind(do.call("rbind", lapply(1:(M-1), function(x){      # Vi laver et datasæt, hvor individer optræder flere gange med hver af deres potentielle værdier for W
    ko <- subset(data, statusW == 0 & W == x)                         # Vi fokuserer på de individer hvor vi ikke kender den sande værdi af W
    do.call("rbind", lapply(x:M, function(y) transform(ko, W = y)))   # Tilføj dem til vores datasæt med alle potentielle værdier af W
  })), subset(data, statusW == 1))                                    # Vi vil trods alt også godt have individer i vores datasæt hvor vi kender den sande værdi af W
  bigdata <- bigdata2 <- transform(bigdata, weight = 1)               # Vi vil lave en vægtet Cox model. Vægte er per default 1, men ændres om nogle linjers kode
  itt <- 0                                                            # Antal iterationer initialiseres
  while(crit & itt < maxit){
    itt <- itt + 1
    ## cat(itt, "\n")
    prPrev <- probs
    ## Likelihood bidrag starter med at være survivalfunktions bidraget. For obs med event skal senere ganges med hazard
    bigdata2 <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))
    ## Nedenunder ganger vi med hazard for obs med event
    # bigdata$L[bigdata$status == 1] <- bigdata$L[bigdata$status == 1] *
    # bas$baseline[sapply(bigdata$time2[bigdata$status == 1],
    # function(x) which.min(abs(x - bas$time)))] *
    # predict(model, bigdata[bigdata$status==1,], type = "risk")
    
    bigdata2$L[bigdata2$status == 1] <- bigdata2$L[bigdata2$status == 1] * 
      (predict(model, newdata = bigdata2[bigdata2$status == 1,], type = "expected") -
         predict(model, newdata = transform(bigdata2[bigdata2$status == 1,], time2 = time2-1e-7), type = "expected"))
    if(sex){
      help <- aggregate(x = bigdata2$L, by = list(bigdata2$id, bigdata2$W, bigdata2$sex), FUN = prod)
      help <- transform(help, x = x * probs[Group.2+Group.3*M])
      help2 <- aggregate(help[,4], by = list(help[,1]), FUN = sum)
      help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
      bigdata2$weight <- help$x[sapply(1:nrow(bigdata2),
                                       function(i) which(bigdata2$id[i] == help[,1] & bigdata2$W[i] == help[,2]))]  
    }
    else{
      ## Likelihood bidrag fra samme individ og med samme værdi af W ganges sammen for at få samlet likelihood bidrag fra hvert individ per værdi af W
      help <- aggregate(x = bigdata2$L, by = list(bigdata2$id, bigdata2$W), FUN = prod)
      ## Der ganges marginal fordeling af W på, da vægte per Bayes er givet ved L(W) * p(W) / sum(L(w) * p(w))
      help <- transform(help, x = x * probs[Group.2])
      ## Udregn sum(L(w)*p(w)) for hvert individ
      help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
      ## Udregn vægte
      help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
      ## Tilføj vægte til datasæt
      bigdata2$weight <- help$x[sapply(1:nrow(bigdata2),
                                       function(i) which(bigdata2$id[i] == help[,1] & bigdata2$W[i] == help[,2]))]  
    }
    bigdata2 <- bigdata2[weight > 0,]
    ## Fit vægtet Cox
    model <- coxph(formula, data = bigdata2, weights = weight, model = TRUE)
    ## Skaf ny baseline hazard
    # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))
    ## Estimer marginal fordeling af W ikke-parametrisk som gennemsnit af vægte
    if(sex){
      probs[,1] <- do.call("c", lapply(1:M, function(i){
        indx <- which(bigdata2$W == i & bigdata2$sex == 0)
        tmp1 <- sum(tapply(bigdata2$weight[indx], bigdata2$id[indx], unique))
        indx <- which(bigdata2$sex == 0)
        tmp2 <- sum(unlist(tapply(bigdata2$weight[indx], bigdata2$id[indx], unique)))
        tmp1 / tmp2
      })) 
      probs[,2] <- do.call("c", lapply(1:M, function(i){
        indx <- which(bigdata2$W == i & bigdata2$sex == 1)
        tmp1 <- sum(tapply(bigdata2$weight[indx], bigdata2$id[indx], unique))
        indx <- which(bigdata2$sex == 1)
        tmp2 <- sum(unlist(tapply(bigdata2$weight[indx], bigdata2$id[indx], unique)))
        tmp1 / tmp2
      }))
    }
    else{
      probs <- do.call("c", lapply(1:M, function(i){
        indx <- which(bigdata2$W == i)
        sum(tapply(bigdata2$weight[indx], bigdata2$id[indx], unique))
      })) / n 
    }
    ## Error er blot absolut værdien af ændringen i estimatet af marginal fordeling af W
    err <- abs(probs - prPrev)
    ## Hvis ændringen er tilpas stor fortsættes estimationen
    crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
    ## cat(probs, "\n")
  }
  cat("estimation iterations: ", itt, "\n")
  list(model = model, probs = probs, d = bigdata2)                                  # Returner model objekt samt estimat af marginal fordeling af W
}

# EM algorithm variance ------------------------------------------------------------------------------------------                   
emSE <- function(object, M = 4, maxit = 1000, sex = FALSE){
  eps1 <- 1e-10
  eps2 <- 0#1e-6
  bigdata <- object$d
  n <- length(unique(bigdata$id))
  # d <- 1 / n
  d <- 5.84e-6
  # d <- 1 / nrow(bigdata)
  probs <- object$probs
  model <- object$model
  formula <- formula(model)
  ## w <- model$weights
  ## X <- model.matrix(model)
  beta <- coef(model)
  ## time1 <- model$y[,1]
  ## time2 <- model$y[,2]
  ## status <- model$y[,3]
  npar <- length(beta)
  ## means <- m$means
  # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))
  information <- matrix(NA, nrow = npar, ncol = npar)
  for(i in 1:npar){
    beta2 <- beta
    beta2[i] <- beta2[i] + d
    crit <- TRUE                                                        # Et kriterie relateret til ændringen i den marginale fordeling af W
    itt <- 0
    while(crit & itt < maxit){
      itt <- itt + 1
      ## cat(itt, "\n")
      prPrev <- probs
      bigdata <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))
      ## Nedenunder ganger vi med hazard for obs med event
      bigdata$L[bigdata$status == 1] <- bigdata$L[bigdata$status == 1] *
        (predict(model, newdata = bigdata[bigdata$status == 1,], type = "expected") -
           predict(model, newdata = transform(bigdata[bigdata$status == 1,], time2 = time2-1e-7), type = "expected"))
      ## Likelihood bidrag fra samme individ og med samme værdi af W ganges sammen for at få samlet likelihood bidrag fra hvert individ per værdi af W
      if(sex){
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W, bigdata$sex), FUN = prod)
        help <- transform(help, x = x * probs[Group.2+Group.3*M])
        help2 <- aggregate(help[,4], by = list(help[,1]), FUN = sum)
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      else{
        ## Likelihood bidrag fra samme individ og med samme værdi af W ganges sammen for at få samlet likelihood bidrag fra hvert individ per værdi af W
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W), FUN = prod)
        ## Der ganges marginal fordeling af W på, da vægte per Bayes er givet ved L(W) * p(W) / sum(L(w) * p(w))
        help <- transform(help, x = x * probs[Group.2])
        ## Udregn sum(L(w)*p(w)) for hvert individ
        help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
        ## Udregn vægte
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        ## Tilføj vægte til datasæt
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      ## Fit vægtet Cox
      model <- coxph(formula, data = bigdata, weights = weight, control = coxph.control(iter.max=0),
                     init = beta2, model = TRUE)
      # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))
      ## Estimer marginal fordeling af W ikke-parametrisk som gennemsnit af vægte
      if(sex){
        probs[,1] <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i & bigdata$sex == 0)
          tmp1 <- sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
          indx <- which(bigdata$sex == 0)
          tmp2 <- sum(unlist(tapply(bigdata$weight[indx], bigdata$id[indx], unique)))
          tmp1 / tmp2
        })) 
        probs[,2] <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i & bigdata$sex == 1)
          tmp1 <- sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
          indx <- which(bigdata$sex == 1)
          tmp2 <- sum(unlist(tapply(bigdata$weight[indx], bigdata$id[indx], unique)))
          tmp1 / tmp2
        }))
        
      }
      else{
        probs <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i)
          sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
        })) / n
      }
      ## Error er blot absolut værdien af ændringen i estimatet af marginal fordeling af W
      err <- abs(probs - prPrev)
      ## Hvis ændringen er tilpas stor fortsættes estimationen
      crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
    }
    first <- model$first
    cat("Number of iterations: ", itt, "\n")
    beta2 <- beta
    beta2[i] <- beta2[i] - d
    crit <- TRUE                                                        # Et kriterie relateret til ændringen i den marginale fordeling af W
    itt <- 0
    while(crit & itt < maxit){
      itt <- itt + 1
      ## cat(itt, "\n")
      prPrev <- probs
      bigdata <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))
      ## Nedenunder ganger vi med hazard for obs med event
      bigdata$L[bigdata$status == 1] <- bigdata$L[bigdata$status == 1] *
        (predict(model, newdata = bigdata[bigdata$status == 1,], type = "expected") -
           predict(model, newdata = transform(bigdata[bigdata$status == 1,], time2 = time2-1e-7), type = "expected"))
      ## Likelihood bidrag fra samme individ og med samme værdi af W ganges sammen for at få samlet likelihood bidrag fra hvert individ per værdi af W
      if(sex){
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W, bigdata$sex), FUN = prod)
        help <- transform(help, x = x * probs[Group.2+Group.3*M])
        help2 <- aggregate(help[,4], by = list(help[,1]), FUN = sum)
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      else{
        ## Likelihood bidrag fra samme individ og med samme værdi af W ganges sammen for at få samlet likelihood bidrag fra hvert individ per værdi af W
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W), FUN = prod)
        ## Der ganges marginal fordeling af W på, da vægte per Bayes er givet ved L(W) * p(W) / sum(L(w) * p(w))
        help <- transform(help, x = x * probs[Group.2])
        ## Udregn sum(L(w)*p(w)) for hvert individ
        help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
        ## Udregn vægte
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        ## Tilføj vægte til datasæt
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      ## Fit vægtet Cox
      model <- coxph(formula, data = bigdata, weights = weight, control = coxph.control(iter.max=0),
                     init = beta2, model = TRUE)
      # bas <- transform(basehaz(model), baseline = diff(c(0,hazard)))
      ## Estimer marginal fordeling af W ikke-parametrisk som gennemsnit af vægte
      if(sex){
        probs[,1] <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i & bigdata$sex == 0)
          tmp1 <- sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
          indx <- which(bigdata$sex == 0)
          tmp2 <- sum(unlist(tapply(bigdata$weight[indx], bigdata$id[indx], unique)))
          tmp1 / tmp2
        })) 
        probs[,2] <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i & bigdata$sex == 1)
          tmp1 <- sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
          indx <- which(bigdata$sex == 1)
          tmp2 <- sum(unlist(tapply(bigdata$weight[indx], bigdata$id[indx], unique)))
          tmp1 / tmp2
        }))
        
      }
      else{
        probs <- do.call("c", lapply(1:M, function(i){
          indx <- which(bigdata$W == i)
          sum(tapply(bigdata$weight[indx], bigdata$id[indx], unique))
        })) / n
      }
      ## Error er blot absolut værdien af ændringen i estimatet af marginal fordeling af W
      err <- abs(probs - prPrev)
      ## Hvis ændringen er tilpas stor fortsættes estimationen
      crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
    }
    cat("Number of iterations: ", itt, "\n")
    second <- model$first
    ## information[i,i:npar] <- information[i:npar,i] <- 1/d * model$first[i:npar]
    information[i,] <- 1/(2*d) * (first - second)
  }
  sqrt(diag(solve(-information)))
  # -information
}

ate <- function(object, timepoint, D){
  probs <- object$probs
  m <- object$model
  d <- data.frame(treat = rep(0, length(probs)), time1 = 0, time2 = timepoint, W = 1:length(probs), status = 1)
  p1 <- sum((1 - predict(m, newdata = d, type = "survival")) * probs)
  d <- data.frame(treat = rep(c(0,1), length(probs)), 
                  time1 = rep(c(0, D), length(probs)), 
                  time2 = rep(c(D, timepoint), length(probs)), 
                  W = rep(1:length(probs), each = 2), status = rep(c(0,1), length(probs)))
  p0 <- aggregate(predict(m, newdata = d, type = "survival"), by = list(d$W), FUN = "prod")$x
  p0 <- sum((1 - p0) * probs)
  list(ate = p1 - p0, risk1 = p1, risk0 = p0)
}

stovlestrop <- function(object, data, timepoint, D, M = 4, nstraps = 100, maxit = 100){
  strapped_ests <- numeric(nstraps)
  for(i in 1:nstraps){
    cat("Iteration: ", i, "\n")
    # ids <- sample(unique(data$id), replace = TRUE)
    # d <- do.call(rbind, lapply(ids, function(j) data[id == j,]))
    d <- dsample(data, ~id)
    m <- em(formula(object$model), data = d, M = 4, maxit = 100)#, object = object)
    cat("Cox: ", m$model$coefficients, "\n")
    strapped_ests[i] <- ate(m, timepoint, D)$ate
    cat("est:", i, ":", strapped_ests[i], "\n")
  }
  strapped_ests
}

stovlestrop_par <- function(object, data, timepoint, D, M = 4, nstraps = 100, maxit = 100){
  strapped_ests <- foreach(i=1:nstraps, .combine="c", .export = c("em","ate","dsample"), .packages = c("data.table","survival","stats")) %dopar% {
    cat("Iteration: ", i, "\n")
    # ids <- sample(unique(data$id), replace = TRUE)
    # d <- do.call(rbind, lapply(ids, function(j) data[id == j,]))
    d <- dsample(data, ~id)
    m <- em(formula(object$model), data = d, M = M, maxit = 100)#, object = object)
    ate(m, timepoint, D)$ate
  }
  strapped_ests
}

cox_ate <- function(object, data, timepoint, D){
  firsts <- data[, match(unique(id), id)]
  dsub <- data[firsts,]
  dsub[, time2 := timepoint]
  p1 <- mean(1 - predict(object, newdata = transform(dsub, treat = 0), type = "survival"))
  dsub <- rbind(dsub, dsub)
  dsub <- dsub[order(id),]
  dsub[, time2 := rep(c(D, timepoint), nrow(dsub)/2)]
  dsub[, time1 := rep(c(0, D), nrow(dsub)/2)]
  dsub[, treat := rep(c(0, 1), nrow(dsub)/2)]
  p0 <- aggregate(predict(object, newdata = dsub, type = "survival"), by = list(dsub$id), FUN = "prod")$x
  p0 <- mean(1 - p0)
  list(ate = p1 - p0, risk1 = p1, risk0 = p0)
}

cox_stovlestrop <- function(formula, data, timepoint, D, nstraps = 100){
  strapped_ests <- numeric(nstraps)
  for(i in 1:nstraps){
    cat("Iteration: ", i, "\n")
    d <- dsample(data, ~id)
    # ids <- sample(unique(data$id), replace = TRUE)
    # d2 <- data[,,by=id]
    # d <- do.call(rbind, lapply(ids, function(j) data[id == j,]))
    # d <- data.frame(d)
    # setDT(d)
    m <- coxph(formula, data = d, model = TRUE) #, object = object)
    cat("Cox: ", m$model$coefficients, "\n")
    strapped_ests[i] <- cox_ate(m, data = d, timepoint = timepoint, D = D)$ate
    cat("est:", i, ":", strapped_ests[i], "\n")
  }
  strapped_ests
}



