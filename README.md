# BayesianECCCGARCH
### Bayesian Estimation and Inference for the ECCC-GARCH Model in R

by Tomasz Woźniak

> **Summary.** A random walk Metropolis-Hastings algorithm for the Bayesian estimation of the Vector Autoregressive models with the conditional volatility process being the Extended Constant Conditional Correlation GARCH(1,1) model is provided, as well as an appropriate estimator for the marginal data density. 
>
> **Keywords.** R, Multivariate GARCH Models, Metropolis-Hastings Sampler, Marginal Data Density

## Citations

To refer to the code in publications, please, cite one of the following papers:

> Woźniak, Tomasz (2015) Testing Causality Between Two Vectors in Multivariate GARCH Models, *International Journal of Forecasting*, 31(3), pp. 876--894, DOI: [10.1016/j.ijforecast.2015.01.005](http://doi.org/10.1016/j.ijforecast.2015.01.005).

> Woźniak, Tomasz (2018) Granger-Causal Analysis of GARCH Models: a Bayesian Approach, *Econometric Reviews*, 37(4), pp. 325-346, DOI: [10.1080/07474938.2015.1092839](http://doi.org/10.1080/07474938.2015.1092839).

## Project contents

The project contains the following files:

- [`BayesianECCCGARCH.pdf`](https://gitlab.com/tomaszwozniak/BayesianECCCGARCH/-/blob/master/BayesianECCCGARCH.pdf) a vignette presenting the model and R code
- `BayesianECCCGARCH.R` containing the utility functions
- `BayesianECCCGARCH-example.R` presenting an application of the functions in a simple example
- `priorConstant-M0-shrinkage.R` reproduction of the results from the paper for one of the models
- `priorConstant-M1-shrinkage.R` reproduction of the results from the paper for another model
- `dataECB.RData` data used in the paper

## Downloading the code

To download the code simply click on the download icon on the top of this page

![](gl-download.png)

and select the format of the compressed file to be downloaded.

## Forking and contributing 

You can also choose to fork the project and clone your copy from the repository. Just follow the steps (requires a [GitLab account](https://gitlab.com)):

1. **On this page:** fork the project by clicking the icon on the top of this page ([more info](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html#creating-a-fork))

   ![](gl-fork.png)

2. **On you computer:** clone the repository you have just forked to create its local copy that you can work with.

3. **Optional:** if you find a bug or if you improve the code, please feel free to submit a pull/merge request. ([more info](https://docs.gitlab.com/ee/topics/gitlab_flow.html#mergepull-requests-with-gitlab-flow))

