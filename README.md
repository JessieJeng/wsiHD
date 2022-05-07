# wsiHD
Weak Signal Inference Under High Dimensionality
X. Jessie Jeng and Shannon Holloway
wsiHD is an R package that implements an analytic framework for weak signal estimation and inclusion under arbitrary covariance dependence. The framework comprises three steps, the estimation of the bounding sequences, the estimation of the signal proportion, and the estimation and control of the false negative proportion (FNP). An illustrative example based on the dataset provided with the package is also presented. 

# Examples
# Data
We use the dataset provided with the package, wsiData, to illustrate a typical analysis. This dataset is a publicly available genomic dataset (Buhlmann et al. (2014)) providing gene expression levels and the rate of riboflavin production with Bacillus subtilis for 71 individuals. The data comprises 4088 gene expression levels and the logarithm of the riboflavin production rate ($q_RRIBFLC). The data can be loaded in the usual way
data(wsiData)
dim(wsiData)
# [1] 71 4089
Consider the summary statistics of only the outcome
summary(wsiData$q_RIBFLC)
# Min. 1st Qu. Median Mean 3rd Qu. Max.
# -9.966 -7.688 -6.948 -7.159 -6.449 -5.673
We see that the range of the logarithm of the riboflavin production rate, y, is -9.9658 ≤ y ≤ -5.673. The
range of the gene expression levels, `, is 3.3228 ≤ ` ≤ 14.3896.
summary(c(data.matrix(wsiData[,-1L])))
# Min. 1st Qu. Median Mean 3rd Qu. Max.
# 3.323 6.946 7.665 7.669 8.352 14.390

# P-values
We first obtain the test statistics using marginal regression
dm <- data.matrix(frame = wsiData)
pval <- apply(X = dm[,-1L],
MARGIN = 2L,
FUN = function(x,y) {
summary(object = lm(formula = y~x))$coef[2L,4L]
},
y = dm[,1L])
p <- length(x = pval)

Next, we obtain 1000 sets of samples from the null distribution and their corresponding p-values
n <- 1000L
sig <- stats::cor(x = dm[,-1L])
zz <- MASS::mvrnorm(n = n, mu = rep(x = 0.0, times = p), Sigma = sig)
pval_null <- t(x = {1.0 - stats::pnorm(q = abs(x = zz))}*2.0)

where we have transposed the p-value matrix to put it into the expected input format.
Though our data is not of sufficient dimension to warrant the use of more memory efficient storage and access,
we will define variables of class “big.matrix” and “ff_matrix” for illustration purposes.
pval_nullBM <- as.big.matrix(pval_null, type = "double")
pval_nullFF <- ff(vmode = "double", dim = c(p,n), pval_null)
cSeq()
The first step of the framework is to estimate the bounding sequences. We will set α = 0.2.
Using the standard R “matrix” object, the call takes the form
cs <- cSeq(pval_null = pval_null, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
A message is generated indicating the number of samples (n) and variables (p). An S3 object of class “wsiHD”
comprising a list object is returned with element $c05 and $c1.
cs
## Bounding Sequences
## c05: 0.3165
## c1: 0.9809
For the “big.matrix” and “ff_matrix” objects,
cSeq(pval_null = pval_nullBM, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3166
## c1: 0.9809
cSeq(pval_null = pval_nullFF, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3164
## c1: 0.9809
yield similar return objects. Notice, however, that the estimates are not identical. Any differences are due to
the algorithm used to obtain the (1 − α)-th quantile. In base R, there are nine algorithms available to obtain
quantiles. We have opted to use the default algorithm in this implementation. However, this default algorithm
is not the same as that implemented by Armadillo (the underpinnings for the “big.matrix” implementation),
5
which uses type = 5, nor that of ff, which uses type = 1. Thus, estimates of the bounding sequences might
differ slightly for inputs of different classes but with equivalent p-values. For large p, any differences will be
very small.
signalProp()
The second step of the framework is to use the estimated bounding sequences to obtain the estimated signal
proportion.
piHat <- signalProp(pval = pval, c05 = cs$c05, c1 = cs$c1)
An S3 object of class “wsiHD” comprising a list object is returned containing πˆ ($piHat), πˆ0.5 ($piHat05),
and πˆ1 ($piHat1).
piHat
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0791
## piHat1: 0.0978
These results can be equivalently obtained using a slightly different input structure. Namely,
signalProp(pval = pval, pval_null = pval_null, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3165
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0791
## piHat1: 0.0978
Here we see that the bounding sequences were estimated internally and are returned through the value object.
Again, we see that using the alternative input classes leads to slightly different results due to the underlying
quantile algorithms.
signalProp(pval = pval, pval_null = pval_nullBM, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3166
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0791
## piHat1: 0.0978
signalProp(pval = pval, pval_null = pval_nullFF, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3164
6
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0792
## piHat1: 0.0978
fnpOpt()
The final step of the framework is to use the estimated signal proportion to obtain the estimated number of
signal variables and then estimate the FNP. We will take β = 0.1 in this example.
sHat <- ceiling(x = p*piHat$piHat)
fnp <- fnpOpt(pval = pval, beta = 0.1, sHat = sHat)
An S3 object of class “wsiHD” comprising a list object is returned containing $ind, the rank satisfying the
threshold condition; $pvalue, the maximum p-value for the variables that satisfy the threshold condition; and
$FNP, the estimated FNP for all variables. Note that the class of the $FNP element will depend on the class
input pval.
fnp
## Number of Signals: 398
## max p-value: 0.01049627
## Summary of FNP:
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.00000 0.00000 0.00000 0.05199 0.00000 0.99750
Though the full p-dimensional vector of F NP \ is returned, the print function only displays the summary
statistics.
Similar in spirit to the signalProp() function, these results can be equivalently obtained using a slightly
different input structure. Namely,
fnpOpt(pval = pval, pval_null = pval_null, beta = 0.1, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3165
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0791
## piHat1: 0.0978
## Number of Signals: 398
## max p-value: 0.01049627
## Summary of FNP:
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.00000 0.00000 0.00000 0.05199 0.00000 0.99750
Here we see that the bounding sequences were estimated internally as well as the estimated signal proportions.
These are returned through the value object.
And using the alternative input classes
7
fnpOpt(pval = pval, pval_null = pval_nullBM, beta = 0.1, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3166
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0791
## piHat1: 0.0978
## Number of Signals: 398
## max p-value: 0.01049627
## Summary of FNP:
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.00000 0.00000 0.00000 0.05199 0.00000 0.99750
fnpOpt(pval = pval, pval_null = pval_nullFF, beta = 0.1, alpha = 0.2)
## estimating bounding sequence using 1000 samples, each containing 4088 variables
## Bounding Sequences
## c05: 0.3164
## c1: 0.9809
## Signal Proportions
## piHat: 0.0978
## piHat05: 0.0792
## piHat1: 0.0978
## Number of Signals: 398
## max p-value: 0.01049627
## Summary of FNP:
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.00000 0.00000 0.00000 0.05199 0.00000 0.99750
