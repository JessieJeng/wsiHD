## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
opt <- options()
options(continue="  ", width=70, prompt=" ")
on.exit(options(opt))
library(wsiHD, quietly=TRUE)
library(bigmemory, quietly = TRUE)
library(ff)

## ----eval=FALSE---------------------------------------------------------------
#  cSeq(pval_null, alpha = 0.1)

## ----eval=FALSE---------------------------------------------------------------
#  signalProb(pval, ..., c05, c1)

## ----eval=FALSE---------------------------------------------------------------
#  signalProb(pval, pval_null, ..., alpha = 0.1)

## ----eval=FALSE---------------------------------------------------------------
#  fnpOpt(pval, ..., beta, sHat)

## ----eval=FALSE---------------------------------------------------------------
#  fnpOpt(pval, pval_null, ..., beta, alpha = 0.1)

## -----------------------------------------------------------------------------
data(wsiData)

## -----------------------------------------------------------------------------
dim(wsiData)

## -----------------------------------------------------------------------------
summary(wsiData$q_RIBFLC)

## -----------------------------------------------------------------------------
summary(c(data.matrix(wsiData[,-1L])))

## -----------------------------------------------------------------------------
dm <- data.matrix(frame = wsiData)
pval <- apply(X = dm[,-1L], 
              MARGIN = 2L, 
              FUN = function(x,y) {
                summary(object = lm(formula = y~x))$coef[2L,4L]
              }, 
              y = dm[,1L])
p <- length(x = pval)

## -----------------------------------------------------------------------------
n <- 1000L
sig <- stats::cor(x = dm[,-1L])
zz <- MASS::mvrnorm(n = n, mu = rep(x = 0.0, times = p), Sigma = sig)
pval_null <- t(x = {1.0 - stats::pnorm(q = abs(x = zz))}*2.0)

## -----------------------------------------------------------------------------
pval_nullBM <- as.big.matrix(pval_null, type = "double")
pval_nullFF <- ff(vmode = "double", dim = c(p,n), pval_null)

## -----------------------------------------------------------------------------
cs <- cSeq(pval_null = pval_null, alpha = 0.2)

## -----------------------------------------------------------------------------
cs

## -----------------------------------------------------------------------------
cSeq(pval_null = pval_nullBM, alpha = 0.2)
cSeq(pval_null = pval_nullFF, alpha = 0.2)

## -----------------------------------------------------------------------------
piHat <- signalProp(pval = pval, c05 = cs$c05, c1 = cs$c1)

## -----------------------------------------------------------------------------
piHat

## -----------------------------------------------------------------------------
signalProp(pval = pval, pval_null = pval_null, alpha = 0.2)

## -----------------------------------------------------------------------------
signalProp(pval = pval, pval_null = pval_nullBM, alpha = 0.2)
signalProp(pval = pval, pval_null = pval_nullFF, alpha = 0.2)

## -----------------------------------------------------------------------------
sHat <- ceiling(x = p*piHat$piHat)
fnp <- fnpOpt(pval = pval, beta = 0.1, sHat = sHat)

## -----------------------------------------------------------------------------
fnp

## -----------------------------------------------------------------------------
fnpOpt(pval = pval, pval_null = pval_null, beta = 0.1, alpha = 0.2)

## -----------------------------------------------------------------------------
fnpOpt(pval = pval, pval_null = pval_nullBM, beta = 0.1, alpha = 0.2)
fnpOpt(pval = pval, pval_null = pval_nullFF, beta = 0.1, alpha = 0.2)

