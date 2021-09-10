# File name: 		BayesianECCCGARCH.R
# By: 				Tomasz Wozniak, email: tomasz.wozniak@unimelb.edu.au
# Purpose:      	a set of functions for the Bayesian Extiamtion of VAR-ECCC-GARCH models
# This Version:   21.01.2015


# Codes for trivariate models
                
#   zero.restrictions(param, restrictions=NULL, N, lag)     # creates a vector of parameters including non-restricted parameters and those set to zero  
#	PAR.MAT(param, N, lag, ccc=TRUE, diag=FALSE, s=1)		# from vector to matrices of parameters
#	INV.PAR.MAT(pa,N,lag,ccc=TRUE,diag=FALSE,s=1)
#	ConradKaranasos(pa)										# Checks non-positive conditions for a, A, B matrices for bivariate GARCH process as in Conrad,Karanasos (2010)
#	CHECK(param,N,lag,ccc=TRUE,diag=FALSE,s=1,negative=FALSE,alb=-7)# Checks conditions for parameters
# 	RESID(param,data,lag,ccc=TRUE,diag=FALSE)				# Computes residuals from VAR model
#   VAR.CHANGE(param,T,N,lag,sc,ccc=TRUE,diag=FALSE)
#	VOL(resid,pa)											# Computes series of conditional variances
#	LIKE(param,data,lag,ccc=TRUE,diag=FALSE,sc=NULL)		# Likelihood value for EDCC
#	LIKELIHOOD(param,data,lag,ccc=TRUE,diag=FALSE,sc=NULL,log=TRUE)		# Vectorised version of LIKE
#  deschamps(ni, delta, lambda)        					# a formula for the prior from Deschamps (2006)
#  shrink.garch(pa,restrictions,N,lag, hyper.parameters = c(100,.1) )      # Computes the logarithm of a prior distribution given by:
#	prior(param,priorr,bound,N,lag,ccc=TRUE,diag=FALSE,s=1)				# Logarithm of prior N(0,diag(100))
#	vprior(mpar,priorr,bound,N,lag,ccc=TRUE,diag=FALSE,s=1,log=TRUE)	# vectorised version of prior
#	MH(S,data,par0,sigma0,C=1,lag,ccc=TRUE,diag=FALSE,negative=FALSE,pii,vars=c(30,100),ni=0,sc=NULL)	# Metropolis-Hastings algorithm
#  ml.cj2001(mcmc,kernel,rej="NULL",data,restrictions,sigma0,den,lag=1)       # computes the marginal density of data from MH algorithm output for the ECCC-GARCH model using the method described in Chib, Jeliazkov (2001,JASA)

#  Some utility functions
#   short.mcmc(mcmc,s=100) 
#	vech2m(vec)
#	Q2corr(x,N)
#	lt_sub(ECOV,N,ni)
#	cov_sub(HHCORR,N)


library(coda)
library(fMultivar)
library(fMultivar)
library(mnormt)
library(tmvtnorm)

zero.restrictions = function(param, restrictions=NULL, N, lag){
    # creates a vector of parameters including non-restricted parameters and those set to zero
    # param         - vector of unrestricted parameters
    # restrictions  - vector indicating which parameters are set to zero        
    # N             - number of variables in the system
    # lag           - order of VAR process

    if (!is.vector(restrictions)) {
        z.param = param         # Unrestricted model case
    } else {
        k = (lag + 2.5)*(N^2) + 1.5*N + 1
        if ((length(param) + length(restrictions)) != k) {
            print("ERROR: Inappropriate No. of unrestricted and restricted parameters")
        } else {
            z.param = rep(NA, k)
            index.z = 0
            index.1 = 0

            for (i in 1:length(restrictions)){
                if (index.z == (restrictions[i] - 1)){
                    z.param[restrictions[i]] = 0
                } else {
                    z.param[(index.z + 1):(restrictions[i] - 1)] = param[(index.1 + 1):(restrictions[i] - i)]
                    z.param[restrictions[i]] = 0
                }
                index.z = restrictions[i]
                index.1 = restrictions[i] - i
            }

            if (restrictions[length(restrictions)] != k) {
                z.param[(restrictions[length(restrictions)] + 1):k] = param[(index.1 + 1):length(param)]
            }
        }   
    }

    return(z.param)
}

PAR.MAT = function(para, restrictions=NULL, N, lag){
    # Order: v0,...,vLAG,a,A,B,ni,ccc
    # para          - vector of unrestricted parameters
    # restrictions  - vector indicating which parameters are set to zero

  param = zero.restrictions(para, restrictions, N, lag)

  nv = matrix(0,11,2)
  nv[1,] = c(1,N)
  for (i in 2:11) {nv[i,] = c(nv[i-1,2]+1,nv[i-1,2]+N^2) }
  var = vector("list", lag)
  for (i in 1:lag) {
	var[[i]] = matrix(param[nv[(i+1),1]:nv[(i+1),2]],ncol=N)
  }

  a = param[(N+lag*(N^2)+1):(N+lag*(N^2)+N)]
 
	A = matrix(param[(N+lag*(N^2)+N+1):(N+lag*(N^2)+N+N^2)],ncol=N)
	B = matrix(param[(N+lag*(N^2)+N+N^2+1):(N+lag*(N^2)+N+2*N^2)],ncol=N)
	cc = diag(N)
	
	param.m = new.env()
	if (lag>=0){param.m$v0 = matrix(param[nv[1,1]:nv[1,2]],ncol=1)}
	if (lag>0) {param.m$var = var}
	param.m$a = a
	param.m$A = A
	param.m$B = B
	param.m$ni = param[(2*N+(2+lag)*(N^2)+1)]
	
	cc[lower.tri(cc)==TRUE] = param[(2*N+(2+lag)*(N^2)+2):(2*N+(2+lag)*(N^2)+1+(N-1)*N/2)]
	cc = cc + t(cc)
	diag(cc) = rep(1,N)
	param.m$ccc = cc
  
    param.m <- as.list(param.m)

  return(param.m)
}

'INV.PAR.MAT' = function(pa,restrictions,N,lag){
    cvar = (1:(lag*(N^2)))*0
    for (i in 1:lag) {cvar[((i-1)*N^2+1):(i*N^2)] = vec(pa$var[[i]])}
    para = c(as.vector(pa$v0),cvar,pa$a,vec(pa$A),vec(pa$B),pa$ni,pa$ccc[lower.tri(pa$ccc)==TRUE])
    if (!is.vector(restrictions)) {
    	param = para
    } else {
    	param = para[-restrictions]
    }
    return(param)
}

'CHECK' = function(param,restrictions,N,lag,alb=0.0000001){
  # Check for parameter conditions
  # param - parameters in PAR.MAT(par)
  # N   - time series dimention
  # lag - var process order
  # ccc - CCC indicator

  	pa = PAR.MAT(para = param, restrictions=restrictions, N=N, lag=lag)
  	q = TRUE

  # Positivity cnd.:
  	if (min(pa$ni) < 2) q = FALSE
	if (min(as.vector(pa$a))< alb) q= FALSE

  # Positivity cnd. for GARCH:
	if (!is.vector(restrictions)) {
		if (min(c(pa$A,pa$B)) < alb) q=FALSE	
	} else {
		restrictions.AB = restrictions[(restrictions > (2*N+lag*(N^2)))&(restrictions <= (2*N + (2+lag)*(N^2)))] - (2*N+lag*(N^2))
		AB = c(pa$A,pa$B)[-(restrictions.AB)]
		if (min(AB) < alb) q=FALSE
	}
	
  # Correlations:
  	if (min(pa$ccc[lower.tri(pa$ccc)==TRUE])<=(-1+alb) | max(pa$ccc[lower.tri(pa$ccc)==TRUE])>=(1-alb)) q = FALSE

  # Stationarity
    # Stationarity of VAR:
  	if (lag!=0) {
    	VV <- matrix(0,lag*N,lag*N)
    	for (i in 1:lag) { VV[(1:N),(i-1)*N+(1:N)] = pa$var[[i]] }
    	if (lag>1) {VV[(N+1):(lag*N),1:((lag-1)*N)] = diag(N*(lag-1))}
    	if (max(abs(eigen(VV)$values)) >= 1) q = FALSE 
  	}

  # Stationarity of GARCH:
  	if (max(abs(eigen(pa$A+pa$B)$values)) >= 1) q = FALSE

  	return(q)
}


'RESID' = function(pa,data,lag) {
  # Computes residuals from VAR model
  # param	- vector of parameters of the model
  # Y	- a matrix TxN of data ln(xt/xt-1)
  # N   - number of time series
  # lag - order of VAR process
 
  N = ncol(data)
  T = nrow(data)-lag
  Y = data[(lag+1):nrow(data),]
  Z = matrix(0,T,1+lag*N)	# matrix with a constant and lagged data as in Lutkepohl (2005)
  Z[,1] = 1
  if (lag!=0){
	for (i in 0:(lag-1)) {Z[,(1+i*N+1):(1+i*N+N)] = data[(lag-i):(nrow(data)-(i+1)),]}
  }
  
  PV = matrix(0,(1+lag*N),N)
  PV[1,] = pa$v0
  if (lag!=0) {
	for (i in 1:lag) { PV[(1+(i-1)*N+1):(1+(i-1)*N+N),] = t(pa$var[[i]]) }
  }
  
  resid = Y-Z%*%PV
  return(resid)
}

'VOL' = function(resid,pa) {
  # Computes series of conditional variances
  # resid	- a matrix TxN of RESIDUALS
  # pa - a parameters

	H = resid*0
	H = H + t(apply(0*H,1,"+",pa$a))	# add intercepts
	H = H + (rbind(rep(0,ncol(resid)),resid[1:(nrow(resid)-1),])^2)%*%t(pa$A)	# add lagged squared residuals times A matrix
	
	H[1,] <- H[1,] + (apply(resid,2,sd)^2)%*%t(pa$B)	# add to first cond. variance h0 - var(resid)

	for (i in 2:nrow(resid)) {
		H[i,] <- H[i,] + H[i-1,]%*%t(pa$B)
  	}
  	
	return(H)
}

'LIKE' = function(param,data,restrictions,lag){

  pa = PAR.MAT(para = param, restrictions=restrictions, N=ncol(data), lag=lag)  

  resid = RESID(pa=pa, data=data, lag=lag)
  HT = VOL(resid,pa)
  CC = t(apply(matrix(0,nrow(resid),length(vech(pa$ccc))),1,"+",vech(pa$ccc)))
  COV <- t(apply(cbind(HT,CC),1,cov_sub,N=ncol(data)))

  lt <- try(apply(cbind(resid,COV),1,lt_sub,N=ncol(data),ni=pa$ni))
  l = try(sum(lt))
  
  if (is.numeric(l)) {l=l} else {l=NA}
  return(l)
}

'LIKELIHOOD' <- function(param,data,restrictions,lag,log=TRUE){
  if (is.vector(param)) param <- matrix(param,nrow=1)
  q <- apply(param,1,LIKE,data=data,restrictions=restrictions,lag=lag)
  if (!log) q <- exp(q)
  as.numeric(q)
}

'deschamps' <- function(ni, delta=2, lambda){
    # Computes the logarithm of the Deschamps(2006) prior for the degrees of freedom parameter
  if (ni <= delta) {d <- -Inf} else {d <- log(lambda*exp(-lambda*(ni-delta)))}
  return(d)
}

'shrink.garch' = function(pa,restrictions,N,lag, hyper.parameters = c(100,.1) ){
    # Computes the logarithm of a prior distribution given by:
    # multivariate normal distribution with mean zero and diagonal covariance matrix with variances equal 100 for VAR parameters and constant terms of the GARCH and 0.5 for A and B matrices of the GARCH process
  # hyper.parameters    - hyper-parameters of the prior distribution, 
  # hyper.parameters[1] - variance of the prior for all parameters but A, B matrices
  # hyper.parameters[2] - variance of the prior for A and B matrices  

    # define matrices for prior mean and variance
    pa.prior.mean = pa.prior.var = pa

    # define the number of parameters taken to the prior
	k = (4*N+lag*(N^2))
	k = k - length(restrictions[restrictions<=k])

    # set prior mean and variances for VAR and GARCH parameters only
    pa.prior.mean$v0    = matrix(rep(0, N),ncol=1)
    pa.prior.var$v0     = matrix(rep(hyper.parameters[1], N),ncol=1)

    for (i in 1:lag) pa.prior.mean$var[[i]]   = matrix(rep(0, N^2),ncol=N)
    for (i in 1:lag) pa.prior.var$var[[i]]   = matrix(rep(hyper.parameters[1], N^2),ncol=N)

    pa.prior.mean$a     = rep(0, N)
    pa.prior.var$a      = rep(hyper.parameters[1], N)

    pa.prior.mean$A     = matrix(rep(0, N^2),ncol=N)
    pa.prior.var$A      = matrix(rep(hyper.parameters[2], N^2),ncol=N)

    pa.prior.mean$B     = matrix(rep(0, N^2),ncol=N)
    pa.prior.var$B      = matrix(rep(hyper.parameters[2], N^2),ncol=N)

    # construct the vectors of prior means and variances
    par.prior.mean      = INV.PAR.MAT(pa.prior.mean,restrictions,N=N,lag=lag)[1:k]
    par.prior.var       = INV.PAR.MAT(pa.prior.var,restrictions,N=N,lag=lag)[1:k]

    # compute the value of the prior distribution
    prior.ordinate  =  dmnorm(INV.PAR.MAT(pa,restrictions,N=N,lag=lag)[1:k], mean=par.prior.mean  , varcov=diag(par.prior.var) ,log=T)

    return(prior.ordinate)                                    
}


'prior' <- function(param, restrictions, N, lag, hyper.parameters = c(100,.1)){
  # Computes ln of the prior
  # param - vector of parameters
  # prior - output from: 'Litterman'
  # bound - output from: 'priorBound'
  # hyper.parameters    - hyper-parameters of the prior distribution, 
  # hyper.parameters[1] - variance of the prior for all parameters but A, B matrices
  # hyper.parameters[2] - variance of the prior for A and B matrices  

    pa = PAR.MAT(para = param, restrictions=restrictions, N=N, lag=lag)  
    prior.ni = deschamps(pa$ni, delta=2, lambda=.04) # lambda = .04 gives 32.6% of chances that ni>30 : Integrate[0.04/E^(0.04 (-2 + x)), {x, 30, 100000}] (Wolfram alpha)
    prior.shrink    = shrink.garch(pa, restrictions=restrictions, N=N,lag=lag,hyper.parameters=hyper.parameters)

  return(prior.ni + prior.shrink)
}

'vprior' <- function(mpar, restrictions, N, lag, hyper.parameters = c(100,.1),log=TRUE){
  # vectorized version of prior
  # mpar - a matrix of parameters (in rows the vectors of model's parameters)
  # hyper.parameters    - hyper-parameters of the prior distribution, 
  # hyper.parameters[1] - variance of the prior for all parameters but A, B matrices
  # hyper.parameters[2] - variance of the prior for A and B matrices 

  if (is.vector(mpar)) mpar <- matrix(mpar,nrow=1)
  q <- apply(mpar,1,prior, restrictions=restrictions, N=N, lag=lag, hyper.parameters = hyper.parameters)
  if (!log) q <- exp(q)
  return(as.numeric(q))
}

'MH' = function(S,data,par0,sigma0,C=1,restrictions=NULL, lag=1, hyper.parameters = c(100,0.1), print.iterations=100){
  # Metropolis-Hasting algorithm for MGARCH
  # S                   - size of simulation from posterior
  # data                - data: vector TxN
  # par0                - starting values, as param
  # sigma0              - vcm of candidate
  # C                   - a scaling constant
  # lag		            - lag of VAR
  # restrictions        - vector indicating which parameters are set to zero  
  # hyper.parameters    - hyper-parameters of the prior distribution, 
  # hyper.parameters[1] - variance of the prior for all parameters but A, B matrices
  # hyper.parameters[2] - variance of the prior for A and B matrices

  t0 = proc.time()
  
  sigma0 = C*sigma0
  N = ncol(data)
  
  THETA = matrix(0,S,length(par0))
  KERNEL = rep(0,S)
  LIKELI = rep(0,S)
  PRIOR = rep(0,S)
  rej = rep(0,S)

  # First iteration: check the starting values (if does not work draw new one)
  for (i in 1:S){
    rejj = 0
    if (i==1){
      th.star = par0
      q = CHECK(par0,restrictions=restrictions,N=N,lag=lag)
	  while (!q) {
		th.star = rmst(n=1,xi=par0,Omega=sigma0, alpha=rep(0,length(par0)),nu=5)
		q = CHECK(th.star,restrictions=restrictions,N=N,lag=lag)
        rejj = rejj+1
	  }

	  THETA[1,] = th.star
	  l = LIKELI[1]  = LIKELIHOOD(th.star,data,restrictions,lag,log=TRUE)
	  PRIOR[1] = vprior(th.star, restrictions, N, lag, hyper.parameters =  hyper.parameters, log=TRUE)
	  KERNEL[1] = LIKELI[1] + PRIOR[1]
    } else {
      th.star = rmst(n=1,xi=THETA[i-1,],Omega=sigma0,alpha=rep(0,length(par0)),nu=5)
      q = CHECK(th.star,restrictions=restrictions,N=N,lag=lag)
	  qq = TRUE
  	  while(qq){
		while (!q) {
		  th.star = rmst(n=1,xi=THETA[i-1,],Omega=sigma0,alpha=rep(0,length(par0)),nu=5)
		  q = CHECK(th.star,restrictions=restrictions,N=N,lag=lag)
          rejj = rejj+1
		}
                  
		l = LIKELIHOOD(th.star,data,restrictions,lag,log=TRUE)
		p = vprior(th.star, restrictions, N, lag, hyper.parameters =  hyper.parameters,log=TRUE)
		k = l+p
		qq = is.na(k)
	  }

      if (runif(1)<=min(exp(k-KERNEL[i-1]),1)){
        THETA[i,] = th.star
        KERNEL[i] = k
        LIKELI[i]  = l
        PRIOR[i] = p
      }
      else {
        THETA[i,] = THETA[i-1,]
        KERNEL[i] = KERNEL[i-1]
        LIKELI[i]  = LIKELI[i-1]
        PRIOR[i] = PRIOR[i-1]
      }
      rej[i] = rejj

	  if (i%%50==0){
		temp = new.env()
		temp$THETA = THETA
		temp$KERNEL = KERNEL
		temp$LIKELI = LIKELI
		temp$PRIOR = PRIOR
		temp$i = i
		temp = as.list(temp)
		save(temp,file="m1_tmp.RData")
	  }
    }
    # Print iteration results
    if((i %% print.iterations)==0) cat(" ",i)
  }
  
  t1 = proc.time()
	results = new.env()
	results$THETA = THETA
	results$KERNEL = KERNEL
    results$LIKELI = LIKELI
    results$PRIOR = PRIOR
	results$time = (t0-t1)/3600
    results$rejections = rej
	return(as.list(results))
}






'ml.cj2001' = function(mcmc,kernel,rej=NULL,data,restrictions=NULL,sigma0,den,lag=1){
   ###########################################################################
   # The function computes the marginal density of data from MH algorithm output for 
   # ECCC-GARCH model using the method described in Chib, Jeliazkov (2001,JASA)
   ###########################################################################
   # Inputs:   mcmc - a matrix of draws from the posterior distribution an output of the estimation using MH algorithm - function MH()
   #		kernel - a vector of kernel = log(prior) + log(likelihood) values: an output of the estimation using MH algorithm - function MH()
   #		data - matrix of data
   #		sigma0 - covariance matrix of the candidate generating distribution
   #		den - number of simulation draws for the denumerator
   ###########################################################################
   # Output:	ln(p(y|M)) - logarithm of the marginal density of data
   ###########################################################################
   # Requires:	source("01fnct.R")
   ###########################################################################
   
   # some settings
   N = ncol(data)
   
   # Posterior mean and characteristics
   theta.star = apply(mcmc,2,mean)									# Posterior mean
   theta.star.likelihood = LIKELIHOOD(param = theta.star,data = data, restrictions = restrictions, lag = lag, log = TRUE)	# Likelihood function
   theta.star.prior = vprior(mpar = theta.star, restrictions = restrictions, N = N, lag = lag,log=TRUE)			# Prior distribution
   theta.star.kernel = theta.star.likelihood + theta.star.prior						# Evaluate kernel
   
   # Numerator of formula (9)
   alfa = function(s){return(min(c(s,1)))}                							# a function min{a,1}
   
   to.alpha = matrix(exp(theta.star.kernel - kernel),ncol=1)						# Compute alfa function
   alpha.matrix = apply(to.alpha,1,alfa)
   
   to.q = dmst(mcmc, xi = theta.star, Omega = sigma0,alpha=rep(0,length(theta.star)), nu = 5)			# Evaluate the candidate generating distribution
   
   if (rej==NULL) {
      rej.no   = 0
   } else {
      rej.no   = sum(rej)
   }
   draws.no    = nrow(mcmc)
   numerator = sum(alpha.matrix * to.q)/(rej.no + draws.no)									# Value of the numerator of (9)
   
   # Denumerator of formula (9)
   theta.j = rmst(den,xi=theta.star,Omega=sigma0,alpha=rep(0,length(theta.star)),nu=5)				# Draw from q
   check.indicator = apply(theta.j,1,CHECK,restrictions=restrictions,N=N,lag=lag,alb=0.0000001)
   theta.j = theta.j[check.indicator,]
   
   theta.j.likelihood = LIKELIHOOD(param = theta.j,data = data, restrictions = restrictions, lag = lag, log = TRUE)	# Likelihood function
   theta.j.prior = vprior(mpar = theta.j, restrictions = restrictions, N = N, lag = lag,log=TRUE)				# Prior distribution
   wh =  which(!is.na(theta.j.likelihood)==TRUE)
   # 	cat("Fraction of successful draws: ",length(wh)/den)
   theta.j.kernel = theta.j.likelihood[wh] + theta.j.prior[wh]							# Evaluate kernel
   
   to.alpha.j = matrix(exp(theta.j.kernel - theta.star.kernel),ncol=1)					# Compute alfa function
   alpha.matrix.j = apply(to.alpha.j,1,alfa)
   
   denumerator = sum(alpha.matrix.j)/den									# Value od the denumerator of (9)
   
   # Output:
   ml.log = theta.star.likelihood + theta.star.prior - log(numerator/denumerator)				# Value of (10) Chib, Jeliazkov (2001)   
   # TO DO: Numerical Standard error (Section 2.4 of Chib, Jeliazkov (2001))
   
   return(ml.log)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'short.mcmc' <- function(mcmc,s=100){
  q <- nrow(mcmc)
  short <- mcmc[(1:(q%/%s))*s,]
  as.matrix(short)
}

'vech2m' = function(vec){
	# vech() inverted
	# vec - vector of given parameters
        N = .5*(-1+sqrt(1+8*length(vec)))
	d = matrix(0,N,2)
	d[1,] = c(1,N)
	if (N>1){
		for (i in 2:N){
			d[i,] = c((d[i-1,2]+1),(d[i-1,2]+(N-i+1)))
		}
	}
	final = matrix(0,N,N)
	for (i in 1:(N)){
		final[i,(i):(N)] = vec[d[i,1]:d[i,2]]
		final[(i):(N),i] = vec[d[i,1]:d[i,2]]
	}
	return(final)	
}

'Q2corr' <- function(x,N){vech(cov2cor(vech2m(x)))}

'lt_sub' <- function(ECOV,N,ni){
	q <- try(log(gamma((ni+N)/2)) -log(gamma(ni/2)) -(N/2)*log((ni-2)*pi) -0.5*log(det(vech2m(ECOV[(N+1):length(ECOV)]))) - ((ni+N)/2)*log(1 + (1/(ni-2))*t(ECOV[1:N])%*%solve(vech2m(ECOV[(N+1):length(ECOV)]))%*%ECOV[1:N]))
	if (is.numeric(q)) {q=q} else {q=NA}
	return(q)
}

'cov_sub' <- function(HHCORR,N){
	cov <- vech(diag(sqrt(HHCORR[1:N]))%*%vech2m(HHCORR[(N+1):length(HHCORR)])%*%diag(sqrt(HHCORR[1:N])))
	return(as.vector(cov))
}