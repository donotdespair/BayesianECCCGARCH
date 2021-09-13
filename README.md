# BayesianECCCGARCH
## Bayesian Estimation and Inference for the ECCC-GARCH Model in R

by Tomasz Woźniak

> **Summary.** A random walk Metropolis-Hastings algorithm for the Bayesian estimation of the Vector Autoregressive models with the conditional volatility process being the Extended Constant Conditional Correlation GARCH(1,1) model is provided, as well as an appropriate estimator for the marginal data density. 
>
> **Keywords.** R, Multivariate GARCH Models, Metropolis-Hastings Sampler, Marginal Data Density

### Citations

To refer to the codes in publications, please, cite one of the following papers:

> Woźniak, Tomasz (2015) Testing Causality Between Two Vectors in Multivariate GARCH Models, *International Journal of Forecasting*, 31(3), pp. 876--894, DOI: [10.1016/j.ijforecast.2015.01.005](http://doi.org/10.1016/j.ijforecast.2015.01.005).

> Woźniak, Tomasz (2018) Granger-Causal Analysis of GARCH Models: a Bayesian Approach, *Econometric Reviews*, 37(4), pp. 325-346, DOI: [10.1080/07474938.2015.1092839](http://doi.org/10.1080/07474938.2015.1092839).

### Project contents

The project contains the following files:

- [`BayesianECCCGARCH.pdf`](https://github.com/donotdespair/BayesianECCCGARCH/blob/master/BayesianECCCGARCH.pdf) a vignette presenting the model and R codes
- `BayesianECCCGARCH.R` containing the utility functions
- `BayesianECCCGARCH-example.R` presenting an application of the functions in a simple example
- `priorConstant-M0-shrinkage.R` reproduction of the results from the paper for one of the models
- `priorConstant-M1-shrinkage.R` reproduction of the results from the paper for another model
- `dataECB.RData` data used in the paper

## Downloading the codes

To download the codes follow the steps (requires a [GitHub account](https://github.com) and [GitHub Desktop](https://desktop.github.com/)):

1. **On this page:** fork the project by clicking the icon on the top of this page ([more info](https://guides.github.com/activities/forking/))

   ![](https://github-images.s3.amazonaws.com/help/bootcamp/Bootcamp-Fork.png)

2. **In GitHub Desktop:** go to the application menu to select: `File -> Clone Repository...`, choose the forked repository from the list of available repositories from GitHub.com, select your local folder to store the repository, and click `Clone`. ([even more info](https://docs.github.com/en/get-started/quickstart/fork-a-repo))

3. **Optional:** if you find a bug or if you improve the code, please feel free to submit a pull request. (more info [here](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) and [here](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request))

