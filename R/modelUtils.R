
multinomial <- function(link = "softmax") {
  structure(list(family = "multinomial"), class = "family")
}

fitModel <- function(formula, data, family = "gaussian",
                     weights = NULL, baseline = 1) {
  
  ## this was directly lifted from 'glm'
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## this was lifted from nnet. 
  class.ind <- function(cl) {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
    x[is.na(as.vector(unclass(cl))),] <- NA
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  if (missing(data)) {
    data <- environment(formula)
  }

  mr <- model.response(model.frame(formula, data))
  mm <- model.matrix(formula, data = data)
  
  if (is.factor(mr)) {
    Y <- class.ind(mr)
    if (is.character(baseline)) {
      toDrop <- which(baseline == colnames(Y))
    }
    else {
      toDrop <- baseline
    }

    if (length(toDrop) == 0 || toDrop > ncol(Y)) {
      stop("Bad baseline category.")
    }
    
    baseline <- colnames(Y)[toDrop]
    Y <- Y[, -toDrop, drop = FALSE]
  }
  else {
    Y <- mr
  }
  
  ## we want to match up the colnames from the binomial and
  ## multinomial regression.
  coefNames <-
    if(is.null(ncol(Y)) || ncol(Y) == 1)
      NULL
    else
      colnames(Y)
  
  if (is.null(dim(Y)))
    M <- 1
  else
    M <- ncol(Y)
  
  P <- ncol(mm)
  N <- nrow(mm)
  
  if (is.null(weights)) {
    user.weights <- rep(1, N)
  }
  else {
    user.weights <- as.double(weights)
  }
  
  betas.cur <- double(M*P)
  residuals <- rep(0, M*N)

  model.type <- switch (family$family,
                        "multinomial" = 3,
                        "binomial" = 2,
                        "gaussian" = 1)

  cout <- .C("fit_model_with_allocation",
             "model" = as.integer(model.type),
             "X" = as.double(mm),
             "Y" = as.double(Y),     
             "N" = as.integer(N),
             "P" = as.integer(P),
             "M" = as.integer(M),
             "weights" = as.double(user.weights),
             "estimates" = as.double(betas.cur),
             "loss" = as.double(1),
             "residuals" = as.double(residuals),
             "error" = as.integer(1),
             "tolerance" = as.double(1e-8),
             "max_iter" = as.integer(25))

  if(cout[[11]]!=0)stop("\nError code ",cout[[11]],".\n")

  ## reshape the estimates.
  estimates <- matrix(cout$estimates, nrow = M, byrow = TRUE)
  rownames(estimates) <- coefNames
  
  res <- list("family" = family,
              "X" = mm,
              "Y" = Y,     
              "N" = N,
              "P" = P,
              "M" = M,
              "weights" = user.weights,
              "estimates" = estimates,
              "loss" = cout$loss,
              "residuals" = cout$residuals,
              "error" = cout$error,
              "tolerance" = as.double(1e-8),
              "max_iter" = as.integer(25),
              "formula" = formula,
              "baseline" = baseline)

  ## wrap it up
  class(res) <- "modelUtils"
  return(res)
}

print.modelUtils <- function(x, ...) {
  print(x$family$family)
  print(x$formula)
  print(x$coefficients)
  print(x$loss)
}

##
## This might be a little strange that the type parameter is ignored
## and it always gives you predicted values on the same scale as the
## response.
## 
predict.modelUtils <- function(object, newdata = NULL, ...) {
  predictMulti <- function() {
    sm <- rowSums(etas <- exp(apply(object$estimates, 1, function(r) X %*% r)))
    pis <- etas/(1 + sm)
    nms <- colnames(pis)
    pis <- cbind(pis, 1 - rowSums(pis))
    colnames(pis) <- c(nms, object$baseline)
    pis
  }
  
  X <- if (missing(newdata) || is.null(newdata)) {
    object$X
  }
  else {
    model.matrix(object$formula, data = newdata)
  }

  switch(object$family$family,
         "gaussian" = X %*% t(object$estimates),
         "binomial" = plogis(X %*% t(object$estimates)),
         "multinomial" = predictMulti())
}

plot.modelUtils <- function(x, plot.compare = FALSE, ...) {
  plot(x$estimates)
}

summary.modelUtils <- function(object, ...) {
  class(object) <- "summary.modelUtils"
}

print.summary.modelUtils <- function(x, ...) {
  print(x$family$family)
  print(x$estimates)
}

coefficients.modelUtils <- function(object, ...) {
  return(object$estimates)
}

coef.modelUtils <- function(object, ...) {
  return(object$estimates)
}

##
## these aren't really anything sensible except
## for Gaussian family.
##
residuals.modelUtils <- function(object, ...) {
  residualsMulti <- function() {
     cbind(object$Y, 1 - rowSums(object$Y)) - yhat
   }

   yhat <- predict(object)
   switch(object$family$family,
          "gaussian" = (object$Y - yhat),
          "binomial" = (object$Y - yhat),
          "multinomial" = residualsMulti())

}



##
## what follows is an implementation in R of the C algorithm.
## they are about identical, but this is here mostly as a way
## to try different things quickly - it is absolutely the worst
## way you could write it in R. 
##
.fmi <- function(formula, data, init = NULL, baseline = 1,
                 weights = NULL, max.iter = 20) {
   class.ind <- function(cl) {
     n <- length(cl)
     x <- matrix(0, n, length(levels(cl)))
     x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
     dimnames(x) <- list(names(cl), levels(cl))
     x
   }
   if (missing(data)) {
     data <- environment(formula)
   }
   mr <- model.response(model.frame(formula, data))
   mm <- model.matrix(formula, data = data)
   
   ## get the matrix, drop the first.
   Y <- class.ind(mr)[, -baseline, drop = FALSE]
   
   M <- ncol(Y)
   P <- ncol(mm)
   N <- nrow(mm)

   if (is.null(init))
     betas.cur <- betas.next <- rep(0, M*P)
   else
     betas.cur <- betas.next <- init

   if (is.null(weights))
     user.weights <- rep(1, N)
   else
     user.weights <- weights

   
   betas.cur <- betas.next <- rep(0, M*P)

   X <- as.vector(mm)
   Y <- as.vector(Y)
   
   ## a place to put the vcovs
   bigW <- rep(NA,  N*(M^2))
   
   ## a place to put the mus
   mus <- rep(NA, M*N)

   ## a place to put the cross products
   XTWX <- rep(NA, (M*P)^2)

   ## a place to put the gradient. 
   g.theta <- rep(NA, M*P)

   ## a count of how many times you've gone through the loop.
   iter <- 0

   ## get a value for the likelihood.
   ll.current <- ll.previous <-
     .multinomial.log.likelihood(betas.next, N, P, M, X, Y, weights)
   
   while (TRUE) {
     for (n in 1:N) {
       denom <- 0
     
       for (m in 1:M) {
         eta <- 0 
         for (p in 1:P) {
           eta <- eta + (X[(p-1)*N + n] * betas.cur[p + P*(m-1)])
         }
         mus[(m-1)*N + n] <- exp(eta)
         denom <- denom + mus[(m-1)*N + n] 
       }
     
       for (m in 1:M) {
         mus[(m-1)*N + n] <- ((mus[(m-1)*N + n])/(1 + denom))
       }
     
       for (ii in 1:M) {
         for (jj in 1:M) {
           kk <- (n-1)*(M^2) + (ii-1)*M + jj
           bigW[kk] <- 
             if (ii == jj)
               mus[(ii-1)*N + n]*(1-mus[(ii-1)*N + n])*user.weights[n]
             else
               -mus[(ii-1)*N + n]*(mus[(jj-1)*N + n])*user.weights[n]
         }
       }
     }

     for (m in 1:M) {
       for (mm in 1:M) {
         for (p in 1:P) {
           for (pp in 1:P) {
             sm <- 0
             for (n in 1:N) {
               w <- (n - 1)*(M^2) + (m-1)*M + mm
               sm <- sm + (X[(p-1)*N + n] * bigW[w] * X[(pp-1)*N + n])
             }
             nn <- M*P*(p-1) + P*(mm-1) + (m-1)*(P^2)*M + pp
             XTWX[nn] <- sm
           }
         }
       }
     }
   
     for (m in 1:M) {
       for (p in 1:P) {
         sm <- 0
         for (n in 1:N) {
           sm <- sm + (X[(p-1)*N + n] * user.weights[n] *
                       (Y[(m-1)*N + n] - mus[(m-1)*N + n]))
         }
         g.theta[(m-1)*P + p] <- sm
       }
     }

     XTWX <- as.vector(solve(matrix(XTWX, ncol = M*P, nrow = M*P, byrow = F)))

     alpha <- 1

     while (TRUE) {
       for (m in 1:(M*P)) {
         sm <- 0
         for (p in 1:(M*P)) {
           sm <- sm + (XTWX[(m-1)*M*P + p] * g.theta[p])
         }
         betas.next[m] <- betas.cur[m] + alpha*sm
       }
       ll.trial <- .multinomial.log.likelihood(betas.next, N, P,
                                               M, X, Y, weights)
       
       if (ll.trial >= ll.current) {
         ll.current <- ll.trial
         break
       }
       else if (alpha < 1e-8) {
         stop("decreased alpha to lower limit.") # we can't do better!
       }
       else {
         alpha <- alpha/2
       }
     }

     if (abs(ll.current - ll.previous) < 1e-9) {
       break # we are done!
     }
     else if (iter >= max.iter) {
       break # we are done!
     }
     else {
       ll.previous <- ll.current
       betas.cur <- betas.next
       iter <- iter + 1
     }
   }

   return(list('coefficients' = betas.next,
               'likelihood' = ll.current,
               'iter' = iter))
 }


.multinomial.log.likelihood <- function(betas.trial, N, P, M, X, Y,
                                        weights) {
  ll <- 0

  for (n in 1:N) {
    denom <- sm <- 0
    tmp.mus <- rep(0, M)
    yy <- 0
    
    for (m in 1:M) {
      eta <- 0 
      for (p in 1:P) {
        eta <- eta + (X[(p-1)*N + n] * betas.trial[p + P*(m-1)])
      }
      tmp.mus[m] <- exp(eta)
      denom <- denom + tmp.mus[m]
      yy <- yy + Y[(m-1)*N + n]
    }
    
    for (m in 1:M) {
      tmp.mus[m] <- tmp.mus[m]/(1 + denom)
      sm <- sm + tmp.mus[m]
    }
    
    for (m in 1:M) {
      ll <- ll + Y[(m-1)*N + n] * log(tmp.mus[m])*weights[n]
    }

    ll <- ll + (1 - yy)*log(1 - sm)*weights[n]
  }
  return(ll)
}
