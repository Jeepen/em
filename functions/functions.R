## Estimation
em <- function(formula, data, M = 4, maxit = 1000, object = NULL, sex = FALSE){
  eps1 <- 1e-10
  eps2 <- 0
  n <- length(unique(data$id))                                        # An id variable is assumed to exist
  if(is.null(object)){
    WW <- tapply(data$W, data$id, unique)                               # Number of treatments is assumed to be called W
    statusWW <- tapply(data$statusW, data$id, unique)                   # statusW tells us whether W is observed or not
    probs <- -diff(c(1, survfit(Surv(WW, statusWW) ~ 1)$surv))          # Marginal distribution of W is estimated non-parametrically. Initialized with Kaplan Meier 
    if(sex){
      probs <- -apply(rbind(1, matrix(survfit(Surv(W, statusW) ~ sex, data = data)$surv, ncol = 2, nrow = 4)),2,diff)
      probs <- cbind(probs[,2], probs[,1])
    }    
    model <- coxph(formula, data, model = TRUE, iter = 0)
  }
  else{
    probs <- object$probs
    model <- object$model
  }
  crit <- TRUE                                                        # Convergence criterion initialization
  bigdata <- rbind(do.call("rbind", lapply(1:(M-1), function(x){      # We create a dataset where subjects appear several times for each possible value of W
    ko <- subset(data, statusW == 0 & W == x)                         # We focus on subjects with unknown W
    do.call("rbind", lapply(x:M, function(y) transform(ko, W = y)))   # Add them to dataset with all possible values of W
  })), subset(data, statusW == 1))                                    # Include subjects where W is known again
  bigdata <- bigdata2 <- transform(bigdata, weight = 1)               # Initialize weights
  itt <- 0                                                            # Number of iterations
  while(crit & itt < maxit){
    itt <- itt + 1
    prPrev <- probs
    bigdata2 <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))

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
      help <- aggregate(x = bigdata2$L, by = list(bigdata2$id, bigdata2$W), FUN = prod)
      help <- transform(help, x = x * probs[Group.2])
      help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
      help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
      bigdata2$weight <- help$x[sapply(1:nrow(bigdata2),
                                       function(i) which(bigdata2$id[i] == help[,1] & bigdata2$W[i] == help[,2]))]  
    }
    bigdata2 <- bigdata2[weight > 0,]
    model <- coxph(formula, data = bigdata2, weights = weight, model = TRUE)
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
    err <- abs(probs - prPrev)
    crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
  }
  cat("estimation iterations: ", itt, "\n")
  list(model = model, probs = probs, d = bigdata2)                                  # Returns model object including estimate of marginal distribution of W
}

# EM algorithm variance ------------------------------------------------------------------------------------------                   
emSE <- function(object, M = 4, maxit = 1000, sex = FALSE){
  eps1 <- 1e-10
  eps2 <- 0
  bigdata <- object$d
  n <- length(unique(bigdata$id))
  d <- 5.84e-6
  probs <- object$probs
  model <- object$model
  formula <- formula(model)
  beta <- coef(model)
  npar <- length(beta)
  information <- matrix(NA, nrow = npar, ncol = npar)
  for(i in 1:npar){
    beta2 <- beta
    beta2[i] <- beta2[i] + d
    crit <- TRUE                                                        
    itt <- 0
    while(crit & itt < maxit){
      itt <- itt + 1
      prPrev <- probs
      bigdata <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))
      ## Multiply by hazard when event is observed
      bigdata$L[bigdata$status == 1] <- bigdata$L[bigdata$status == 1] *
        (predict(model, newdata = bigdata[bigdata$status == 1,], type = "expected") -
           predict(model, newdata = transform(bigdata[bigdata$status == 1,], time2 = time2-1e-7), type = "expected"))
      ## Likelihood contribution from same subject and with same value of W multiplied together to get total likelihood contribution per value of W
      if(sex){
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W, bigdata$sex), FUN = prod)
        help <- transform(help, x = x * probs[Group.2+Group.3*M])
        help2 <- aggregate(help[,4], by = list(help[,1]), FUN = sum)
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      else{
        ## Likelihood contribution from same subject and with same value of W multiplied together to get total likelihood contribution
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W), FUN = prod)
        ## Needed for weights in Bayes
        help <- transform(help, x = x * probs[Group.2])
        ## Calculate sum(L(w)*p(w)) for every subject
        help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
        ## Calculate weights
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        ## Add weights to dataset
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      ## Fit weighted Cox
      model <- coxph(formula, data = bigdata, weights = weight, control = coxph.control(iter.max=0),
                     init = beta2, model = TRUE)
      ## Estimate marginal distribution of W non-parametrically
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
      ## Error is absoulute change in probability estimates
      err <- abs(probs - prPrev)
      ## If change is sufficiently big, continue iteration
      crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
    }
    first <- model$first
    cat("Number of iterations: ", itt, "\n")
    beta2 <- beta
    beta2[i] <- beta2[i] - d
    crit <- TRUE                                                        # Criterion for iterations
    itt <- 0
    while(crit & itt < maxit){
      itt <- itt + 1
      prPrev <- probs
      bigdata <- transform(bigdata, L = predict(model, newdata = bigdata, type = "survival"))
      ## Multiply by hazard rate for subjects with the event
      bigdata$L[bigdata$status == 1] <- bigdata$L[bigdata$status == 1] *
        (predict(model, newdata = bigdata[bigdata$status == 1,], type = "expected") -
           predict(model, newdata = transform(bigdata[bigdata$status == 1,], time2 = time2-1e-7), type = "expected"))
      ## Likelihood contribution from same subject and with same value of W multiplied together
      if(sex){
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W, bigdata$sex), FUN = prod)
        help <- transform(help, x = x * probs[Group.2+Group.3*M])
        help2 <- aggregate(help[,4], by = list(help[,1]), FUN = sum)
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      else{
        ## Get likelihood conribution from same subject with same W value 
        help <- aggregate(x = bigdata$L, by = list(bigdata$id, bigdata$W), FUN = prod)
        ## Weights in Bayes are L(W) * p(W) / sum(L(w) * p(w))
        help <- transform(help, x = x * probs[Group.2])
        ## Calculate sum(L(w)*p(w)) for every subject
        help2 <- aggregate(help[,3], by = list(help[,1]), FUN = sum)
        ## Calculate weights
        help$x <- help$x / help2$x[sapply(help$Group.1, function(x) match(1, x == help2$Group.1))]
        ## Add weights to dataset
        bigdata$weight <- help$x[sapply(1:nrow(bigdata),
                                        function(i) which(bigdata$id[i] == help[,1] & bigdata$W[i] == help[,2]))]  
      }
      ## Fit weighted Cox
      model <- coxph(formula, data = bigdata, weights = weight, control = coxph.control(iter.max=0),
                     init = beta2, model = TRUE)
      ## Estimate marginal distribution of W
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
      ## Error is the absolute change in probability estimates
      err <- abs(probs - prPrev)
      ## If change is sufficiently big, continue iteration
      crit <- !all(err < (eps1 * (abs(prPrev) + eps2)))
    }
    cat("Number of iterations: ", itt, "\n")
    second <- model$first
    information[i,] <- 1/(2*d) * (first - second)
  }
  sqrt(diag(solve(-information)))
}

