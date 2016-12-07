NN <-  1000
require(modelUtils)

X = cbind(rep(1,NN), rnorm(NN), runif(NN))
Y = X %*% c(2, 2, 4) + rnorm(NN, 0, 2)

xx <- fitModel(Y ~ X - 1, family = gaussian())
print(coefficients(xx) - coefficients(lm(Y ~ X - 1)))
print(paste("sum is (should be 0)", sum(predict(lm(Y ~ X - 1)) - predict(xx))))
print(paste("sum is (should be 0)", sum(residuals(lm(Y ~ X - 1)) - residuals(xx))))

P = plogis(X %*% c(.5, .2, -1))
P = runif(NN) < P

yy <- fitModel(P ~ X - 1, family = binomial())
print(coefficients(yy))
print(coefficients(yy) - coefficients(glm(P ~ X - 1, family = binomial())))
print(paste("sum is (should be around 0)", sum(predict(glm(P ~ X - 1, family = binomial()),
                                                       type = "response") - predict(yy))))
                                                         
sucfail <- t(sapply(P, function(pp) {
  N <- rpois(1, 10)
  s <- rbinom(1, prob = pp, size = N)
  f <- (N - s)
  c(s, f)
}))

glmTwoCol <- glm(sucfail ~ X - 1, family = binomial())
fModelTwoCol <- fitModel(sucfail ~ X - 1, family = binomial)


##
## sim 1
##
y0 <- runif(NN, -10, 10)
W <- y0/3 + rnorm(NN)

PY <- apply(cbind(W), 1, function(x) {
  bb <- exp(0*x) / ((exp(0*x) + exp(0.5*x) + exp(1.5*x)))
  b1 <- bb * exp(0.5*x)
  b2 <- bb * exp(1.5*x)
  c(bb, b1, b2)
})

Y <- apply(PY, 2, function(x) {
  ru <- runif(1)
  if(ru <= x[1])
    "A"
  else if (ru <= x[2] + x[1])
    "B"
  else
    "C"
})

dta.1 <- data.frame(W, "Y" = factor(Y))

zz <- fitModel(Y ~ ., data = dta.1, family = multinomial())

print(coefficients(zz))
predict(zz)

##
## sim 2
## 
W <- cbind(rep(1, NN), runif(NN, 0, 5), rnorm(NN, 5, 2))

PY <- apply(W, 1, function(x) {
  coefs <- matrix(c(0,0,0,
                    -9, 2, .8,
                    -6, 1, .7,
                    -4, .5, .6), nrow = 4, ncol = 3, byrow = T)

  denom <- sum(apply(coefs, 1, function(cc) {exp(x %*% cc)}))
  as.numeric(apply(coefs, 1, function(cc) {exp(x %*% cc)/denom}))
})

Y <- apply(PY, 2, function(x) {
  ru <- runif(1)
  if (ru <= x[1])
    0
  else if (ru <= x[2] + x[3])
    1
  else if (ru <= x[3] + x[4])
    2
  else
    3
})

dta.2 <- data.frame(W[,-1], "Y" = factor(Y))


ww <- fitModel(Y ~ ., data = dta.2, family = modelUtils:::multinomial(), baseline = 4)
pww <- predict(ww)
coefficients(ww)

pww[1:10,]

#require(VGAM)
#vg <- vglm(Y ~ ., data = dta.2, family = multinomial())
#pvg <- predict(vg, type = "response")


