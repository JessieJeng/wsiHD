# wsiHD
Weak Signal Inference Under High Dimensionality
X. Jessie Jeng and Shannon Holloway
wsiHD is an R package that implements an analytic framework for weak signal estimation and inclusion under arbitrary covariance dependence. The framework comprises three steps, the estimation of the bounding sequences, the estimation of the signal proportion, and the estimation and control of the false negative proportion (FNP). The implementation of each step is described in the following section. An illustrative example based on the dataset provided with the package is also presented.

This function estimates the bounding sequences using the empirical distributions

$$
V_{0.5,a} = \max_{1 \le j \le p} 
\frac{|j/p - p_{nul, (j)}|}
{\sqrt{p_{nul, (j)}}}
~~~~~~ \mathrm{and} ~~~~~~
V_{1,a} = \max_{1 \le j \le p} 
\frac{|j/p - p_{nul, (j)}|}
{p_{nul, (j)} },
$$
