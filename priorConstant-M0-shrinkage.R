# File name: 		00 priorConstant-M0-shrinkage.R
# By: 				Tomasz Wozniak, email: tomasz.wozniak@unimelb.edu.au
# Purpose:      	Computations of the normalizing constraint of the prior distribution, when the normalization is due to the truncation
# First Version:	29.10.2014

####################################################################
# packages and functions
####################################################################
rm(list=ls())

library(truncnorm)
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

'CHECK' = function(param,restrictions,N,lag,alb=0.0000001){
  # Check for parameter conditions
  # param - parameters in PAR.MAT(par)
  # N   - time series dimention
  # lag - var process order
  # ccc - CCC indicator

  	pa = PAR.MAT(para = param, restrictions=restrictions, N=N, lag=lag)
  	q = TRUE

  # Positivity cnd.:
# 	if (min(as.vector(pa$a))< alb) q= FALSE
# 
#   # Positivity cnd. for GARCH:
# 	restrictions.AB = restrictions[(restrictions > (2*N+lag*(N^2)))&(restrictions <= (2*N + (2+lag)*(N^2)))] - (2*N+lag*(N^2))
# 	if (length(restrictions.AB)==0) {
# 		if (min(c(pa$A,pa$B)) < alb) q=FALSE	
# 	} else {
# 		AB = c(pa$A,pa$B)[-(restrictions.AB)]
# 		if (min(AB) < alb) q=FALSE
# 	}
	
  # Stationarity of GARCH:
  	if (max(abs(eigen(pa$A+pa$B)$values)) >= 1) q = FALSE

  	return(q)
}

####################################################################
# compute normalizing constraints for models:
#     + M_0
#     + shinkage
#     + 16 constants all together
####################################################################

S     = 1e7
N     = 2
lag   = 1

ub    = 1.8
hyper = sqrt(0.1)

restrictions         = c(1:6, 17:18)
no.not.restricted    = 8 - length(restrictions[(restrictions > (2*N+lag*(N^2)))&(restrictions <= (2*N + (2+lag)*(N^2)))] - (2*N+lag*(N^2)))

t0 = proc.time()
parameters           = matrix(c(rtruncnorm(n=S*2, a=0, b=Inf, mean = 0, sd = 10), rtruncnorm(n=S*no.not.restricted, a=0, b=ub, mean = 0, sd = hyper) ),ncol=2+no.not.restricted)
t1 = proc.time()
check                = apply(parameters, 1, CHECK, restrictions=restrictions, N=N, lag=lag, alb=0.0000001)
t2 = proc.time()

Pr.simulated         = sum(check)/S

save(Pr.simulated, file= "m00-s.RData")