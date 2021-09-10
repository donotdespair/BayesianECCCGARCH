# BayesianECCCGARCH
## Bayesian Estimation and Inference for the ECCC-GARCH Model in R

by Tomasz Woźniak

> A random walk Metropolis-Hastings algorithm for the Bayesian estimation of the Vector Autoregressive models with the conditional volatility process being the Extended Constant Conditional Correlation GARCH(1,1) model is provided, as well as an appropriate estimator for the marginal data density. 
>
> Keywords: R, Multivariate GARCH Models, Metropolis-Hastings Sampler, Marginal Data Density

The codes are available under the GNU General Public License v3.0. To refer to the codes in publications, please, cite one of the following papers:

Woźniak, Tomasz (2015) Testing Causality Between Two Vectors in Multivariate GARCH Models, *International Journal of Forecasting*, 31(3), pp. 876--894, DOI: [10.1016/j.ijforecast.2015.01.005](http://doi.org/10.1016/j.ijforecast.2015.01.005).

Wo\'zniak, Tomasz (2018) Granger-Causal Analysis of GARCH Models: a Bayesian Approach, *Econometric Reviews*, 37(4), pp. 325-346, DOI: [10.1080/07474938.2015.1092839](http://doi.org/10.1080/07474938.2015.1092839).

The project contains the following files:

- [`BayesianECCCGARCH.pdf`](https://github.com/donotdespair/BayesianECCCGARCH/blob/master/BayesianECCCGARCH.pdf) a vignette presenting the model and R codes
- `BayesianECCCGARCH.R` containing the utility functions
- `BayesianECCCGARCH-example.R` presenting an application of the functions in a simple example
- `priorConstant-M0-shrinkage.R` reproduction of the results from the paper for one of the models
- `priorConstant-M1-shrinkage.R` reproduction of the results from the paper for another model
- `dataECB.RData` data used in the paper
