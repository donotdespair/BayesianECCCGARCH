# File name:      go.R
# By:             Tomasz Wozniak, email: tomasz.wozniak@unimelb.edu.au
# Purpose:        Estimation of the restricted VAR(1)-ECCC-GARCH(1,1) model
# Data:           dataECB.RData
# This Version:  18.03.2016

rm(list=ls())
library(ccgarch); library(coda)
source("BayesianECCCGARCH.R")

load("dataECB.RData")
Y =  gu.after

restrictions = c(10,14)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimation of VAR(1)-ECCC-GARCH(1,1) for  bivariate system
# Restrictions M1: A_21 = B_21 = 0
# Daily exchange rates: GBP/EUR, USD/EUR 
# period AFTER the crisis
# N = 2, T = 1379
# Matrix of data: gu.after  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(35532); s = 100000; C = 0.2
load("aft-M1-dif-00.RData"); qq = as.mcmc(m$THETA[50000:150000,]); sigma0 = cov(qq)
load("aft-M1-dif-00.RData"); qq = as.mcmc(m$THETA); param0 = qq[nrow(qq),];
m = MH(S=s,data=Y,par0=param0,sigma0=sigma0,C=C,restrictions=restrictions,lag=1,hyper.parameters=c(100,.1),print.iterations=1000)
save(m,file="aft-M1-shr-00.RData")

set.seed(35532); s = 200000; C = 0.2
load("aft-M1-dif-00.RData"); qq = as.mcmc(m$THETA[50000:150000,]); sigma0 = cov(qq)
load("aft-M1-shr-00.RData"); qq = as.mcmc(m$THETA); param0 = qq[nrow(qq),];
m = MH(S=s,data=Y,par0=param0,sigma0=sigma0,C=C,restrictions=restrictions,lag=1,hyper.parameters=c(100,.1),print.iterations=1000)
save(m,file="aft-M1-shr-01.RData")

set.seed(234567); s = 200000; C = 0.2
load("aft-M1-dif-00.RData"); qq = as.mcmc(m$THETA[50000:150000,]); sigma0 = cov(qq)
load("aft-M1-shr-01.RData"); qq = as.mcmc(m$THETA); param0 = qq[nrow(qq),];
m = MH(S=s,data=Y,par0=param0,sigma0=sigma0,C=C,restrictions=restrictions,lag=1,hyper.parameters=c(100,.1),print.iterations=1000)
save(m,file="aft-M1-shr-02.RData")

####################################
# MDD
####################################
load("aft-M1-dif-00.RData"); qq = as.mcmc(m$THETA[50000:150000,]); sigma0 = 0.2*cov(qq)

load("aft-M1-shr-00.RData"); mcmc = m$THETA; cernel = m$KERNEL; likeli = m$LIKELI; rej = m$rejections
load("aft-M1-shr-01.RData"); mcmc = rbind(mcmc,m$THETA); cernel = c(cernel,m$KERNEL); likeli = c(likeli,m$LIKELI); rej = c(rej,m$rejections)
load("aft-M1-shr-02.RData"); mcmc = rbind(mcmc,m$THETA); cernel = c(cernel,m$KERNEL); likeli = c(likeli,m$LIKELI); rej = c(rej,m$rejections)

ss    = 10
mcmc  = short.mcmc(mcmc,s=ss); cernel = as.vector(short.mcmc(matrix(cernel,ncol=1),s=ss)); likeli = as.vector(short.mcmc(matrix(likeli,ncol=1),s=ss)); rej = as.vector(short.mcmc(matrix(rej,ncol=1),s=ss))

ml.cj = ml.cj2001(mcmc=mcmc, kernel=cernel, rej=rej, data=Y,restrictions=restrictions,sigma0=sigma0,den=5000,lag=1)
save(ml.cj, file = "mdd-aft-M1-shr.RData")