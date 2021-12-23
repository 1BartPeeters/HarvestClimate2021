

## NOISE-only MODELS ----

###########################################.
# ~ Ricker vs. Beverton-Holt (BH) ----
# ... additive vs. multiplicative noise ----
# 
tmb_RickBH <- '
// 
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(N);
DATA_VECTOR(obs_r);
DATA_SCALAR(penalty);
DATA_INTEGER(model); // 0: Ricker, 1: BH
DATA_INTEGER(error); // 0: Additive, 1: Multiplicative

PARAMETER_VECTOR(z);

// parameters:
PARAMETER_VECTOR(par);
Type a = par[0];
Type b = exp(par[1]);
Type c = exp(par[2]);

Type nll = 0.0; // initialize negative log likelihood
nll -= dnorm(z, 0, 1, true).sum();

Type K = 0.0;

vector<Type> pred_r(N.size());
pred_r.fill(0.0);

if(model == 0){ // Ricker
  if(error == 0){ // Additive
    pred_r = a - b*N + c*z;
  } 
  if(error == 1){ // Multiplicative
    pred_r = a - b*N * exp(-c*z);
  } 
  K = a/b;
}
if (model == 1){ // BH
  if(error == 0){ // Additive
    pred_r = a - log(1+b*N) + c*z;
  } 
  if(error == 1){ // Multiplicative
    pred_r = a - log(1+b*N) * exp(-c*z);
  } 
  K = (exp(a)-1)/b;
}

ADREPORT(a);
ADREPORT(b);
ADREPORT(c);
ADREPORT(K);

// force pred.r == obs.r !
nll -= dnorm(pred_r, obs_r, penalty, true).sum();

//Rcout << eta(0)  << eta(1) <<  sq(0)  << sq(1) << " " << nSS << std::endl;

return nll;
}'
write(tmb_RickBH, file = "TMB_models/tmb_RickBH.cpp")
#
library(TMB)
compile("TMB_models/tmb_RickBH.cpp")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ... BOTH noise ----
# 
tmb_RickBH_bothNoise <- '
// 
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(N);
DATA_VECTOR(obs_r);
DATA_SCALAR(penalty);
DATA_INTEGER(model); // 0: Ricker, 1: BH

PARAMETER_VECTOR(z1);
PARAMETER_VECTOR(z2);

// parameters:
PARAMETER_VECTOR(par);
Type a = par[0];
Type b = exp(par[1]);
Type c = exp(par[2]);
Type d = exp(par[3]);
Type rho = 2/(1 + exp(-par[4])) - 1;

Type nll = 0.0; // initialize negative log likelihood
nll -= dnorm(z1, 0, 1, true).sum();
nll -= dnorm(z2, 0, 1, true).sum();

Type K = 0.0;

vector<Type> z2trans = z1*rho + z2*sqrt(1-rho*rho);
vector<Type> pred_r(N.size());
pred_r.fill(0.0);

if(model == 0){ // Ricker
  pred_r = a - b*N * exp(-c*z1) + d*z2trans;
  K = a/b;
}
if (model == 1){ // BH
  pred_r = a - log(1+b*N)* exp(-c*z1) + d*z2trans;
  K = (exp(a)-1)/b;
}

ADREPORT(a);
ADREPORT(b);
ADREPORT(c);
ADREPORT(K);

// force pred.r == obs.r !
nll -= dnorm(pred_r, obs_r, penalty, true).sum();

//Rcout << eta(0)  << eta(1) <<  sq(0)  << sq(1) << " " << nSS << std::endl;

return nll;
}'
write(tmb_RickBH_bothNoise, file = "TMB_models/tmb_RickBH_bothNoise.cpp")
#
compile("TMB_models/tmb_RickBH_bothNoise.cpp")

library(TMB)
dyn.load(dynlib("TMB_models/tmb_RickBH"))
dyn.load(dynlib("TMB_models/tmb_RickBH_bothNoise"))


AICc.bp <- function(lnL, k, n=length(obs.r)) {
  AIC <- 2*k + 2*lnL # technically -ln(L), but L should be max likelihood, thus already multiplied by -1 in the optim function
  AICc <- AIC + (2*k^2 + 2*k)/(n-k-1)
  return(AICc)
}


# ~ empirical data, example
dat <- read.table("data_TMB/ibex.txt", header=T)

dat$r <- c(diff(log(dat$N)), NA)
dat2 <- dat[complete.cases(dat),]
N <- dat2$N
obs.r <- dat2$r
#
parameters <- list(z=rep(0, length(obs.r)), par=c(0.5, log(0.01), -1))
parameters2 <- list(z1=rep(0, length(obs.r)), z2=rep(0, length(obs.r)), par=c(1, log(0.01), -1, 0, 0))
PENALTY <- 0.001
# input
# model: 0 = Ricker, 1 = BH
# error: 0 = additive environmental noise, 1 = multiplicative noise
# for model with both types of noise, error is not included in input list
input0 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 0, error = 0)
input1 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 0, error = 1)
input2 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 0)
input3 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 1, error = 0)
input4 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 1, error = 1)
input5 <- list(N = N, obs_r = obs.r, penalty=PENALTY, model = 1)
#
obj0 <- MakeADFun(input0, parameters, random="z",DLL = "tmb_RickBH", method = "BFGS",silent = FALSE)
obj1 <- MakeADFun(input1, parameters, random="z",DLL = "tmb_RickBH", method = "BFGS",silent = FALSE)
obj2 <- MakeADFun(input2, parameters2,random=c("z1","z2"),DLL = "tmb_RickBH_bothNoise", method = "BFGS",silent = FALSE)
obj3 <- MakeADFun(input3, parameters, random="z",DLL = "tmb_RickBH", method = "BFGS",silent = FALSE)
obj4 <- MakeADFun(input4, parameters, random="z",DLL = "tmb_RickBH", method = "BFGS",silent = FALSE)
obj5 <- MakeADFun(input5, parameters2,random=c("z1","z2"),DLL = "tmb_RickBH_bothNoise", method = "BFGS",silent = FALSE)
#
obj0$hessian <- obj1$hessian <- obj2$hessian <- obj3$hessian<- obj4$hessian <- obj5$hessian <- FALSE
obj0$control <- obj1$control <- obj2$control <- obj3$control<- obj4$control <- obj5$control <- list("maxit"=10000)
opt0 <- do.call("optim", obj0)
opt1 <- do.call("optim", obj1)
opt2 <- do.call("optim", obj2)
opt3 <- do.call("optim", obj3)
opt4 <- do.call("optim", obj4)
opt5 <- do.call("optim", obj5)
opt0 <- nlminb(start=obj0$par, objective=obj0$fn, gradient=obj0$gr, control= list(iter.max=1000, eval.max= 1500))
opt1 <- nlminb(start=obj1$par, objective=obj1$fn, gradient=obj1$gr, control= list(iter.max=1000, eval.max= 1500))
opt2 <- nlminb(start=obj2$par, objective=obj2$fn, gradient=obj2$gr, control= list(iter.max=1000, eval.max= 1500))
opt3 <- nlminb(start=obj3$par, objective=obj3$fn, gradient=obj3$gr, control= list(iter.max=1000, eval.max= 1500))
opt4 <- nlminb(start=obj4$par, objective=obj4$fn, gradient=obj4$gr, control= list(iter.max=1000, eval.max= 1500))
opt5 <- nlminb(start=obj5$par, objective=obj5$fn, gradient=obj5$gr, control= list(iter.max=1000, eval.max= 1500))
AICc0 <- AICc.bp(opt0$objective, k = length(opt0$par))
AICc1 <- AICc.bp(opt1$objective, k = length(opt1$par))
AICc2 <- AICc.bp(opt2$objective, k = length(opt2$par))
AICc3 <- AICc.bp(opt3$objective, k = length(opt3$par))
AICc4 <- AICc.bp(opt4$objective, k = length(opt4$par))
AICc5 <- AICc.bp(opt5$objective, k = length(opt5$par))
#
aicc <- c(AICc0,AICc1,AICc2,AICc3, AICc4, AICc5)
round(aicc,2)
round(aicc-min(aicc),2)

# estimates
tab0 <-rbind(summary(sdreport(obj0), "report", p.value=T))
tab1 <-rbind(summary(sdreport(obj1), "report", p.value=T))
tab2 <-rbind(summary(sdreport(obj2), "report", p.value=T))
tab3 <-rbind(summary(sdreport(obj3), "report", p.value=T))
tab4 <-rbind(summary(sdreport(obj4), "report", p.value=T))
tab5 <-rbind(summary(sdreport(obj5), "report", p.value=T))
tab0
tab1
tab2
tab3
tab4
tab5

plot(obs.r~N, xlim=c(0, max(N)*1.5), ylim=c(-1,1))
Nvec <- 0:10000
lines(Nvec, tab1[1,1] - tab1[2,1]*Nvec)
lines(Nvec, tab4[1,1] - log(1+tab4[2,1]*Nvec),col=2)
# additive
lines(Nvec, tab0[1,1] - tab0[2,1]*Nvec)
lines(Nvec, tab3[1,1] - log(1+tab3[2,1]*Nvec),col=2)
# both
lines(Nvec, tab2[1,1] - tab2[2,1]*Nvec,col=3)
lines(Nvec, tab5[1,1] - log(1+tab5[2,1]*Nvec),col=4)


#_________________________________________________________________________________________####
# covariate model ############
tmb_RickBH_cov2 <- '
// 
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(N);
DATA_VECTOR(obs_r);
DATA_VECTOR(cov1);
DATA_SCALAR(penalty);
DATA_INTEGER(model); // 0: Ricker, 1: BH
DATA_INTEGER(error); // 0: Additive, 1: Multiplicative
DATA_INTEGER(variable); // 0: Additive, 1: Multiplicative

PARAMETER_VECTOR(z);

// parameters:
PARAMETER_VECTOR(par);
Type a = par[0];
Type b = exp(par[1]);
Type c = exp(par[2]);
Type gamma1 = par[3];

Type nll = 0.0; // initialize negative log likelihood
nll -= dnorm(z, 0, 1, true).sum();

Type K = 0.0;

vector<Type> pred_r(N.size());
pred_r.fill(0.0);

if(model == 0){ // Ricker
  if(error == 0){ // Additive
    if(variable == 0){ // Additive
      pred_r = a - b*N + c*z + gamma1*cov1;
    }
    if(variable == 1){ // Multiplicative
      pred_r = a - b*N * exp(gamma1*cov1) + c*z;
    }
  } 
  if(error == 1){ // Multiplicative
    if(variable == 0){ // Additive
      pred_r = a - b*N * exp(-c*z) + gamma1*cov1;
    }
    if(variable == 1){ // Multiplicative
      pred_r =  a - b*N*exp(-c*z)*exp(gamma1*cov1);
    }
  } 
  K = a/b;
}

if (model == 1){ // BH
  if(error == 0){ // Additive
    if(variable == 0){ // Additive
      pred_r = a - log(1+b*N) + c*z + gamma1*cov1;
    }
    if(variable == 1){ // Multiplicative
      pred_r = a - log(1+b*N) * exp(gamma1*cov1) + c*z;
      
    }
  } 
  if(error == 1){ // Multiplicative
    if(variable == 0){ // Additive
      pred_r = a - log(1+b*N) * exp(-c*z) + gamma1*cov1;
    }
    if(variable == 1){ // Multiplicative
      pred_r = a - log(1+b*N) * exp(gamma1*cov1) * exp(-c*z);
    }
  } 
  K = (exp(a)-1)/b;
}

ADREPORT(a);
ADREPORT(b);
ADREPORT(c);
ADREPORT(gamma1);
ADREPORT(K);

// force pred.r == obs.r !
nll -= dnorm(pred_r, obs_r, penalty, true).sum();

//Rcout << eta(0)  << eta(1) <<  sq(0)  << sq(1) << " " << nSS << std::endl;

return nll;
}'
write(tmb_RickBH_cov2, file = "TMB_models/tmb_RickBH_cov2.cpp")
#
library(TMB)
compile("TMB_models/tmb_RickBH_cov2.cpp")

########################################################-

library(TMB)
dyn.load(dynlib("TMB_models/tmb_RickBH_cov2"))


AICc.bp <- function(lnL, k, n=length(obs.r)) {
  AIC <- 2*k + 2*lnL # technically -ln(L), but L should be max likelihood, thus already multiplied by -1 in the optim function
  AICc <- AIC + (2*k^2 + 2*k)/(n-k-1)
  return(AICc)
}
adjR2 <- function(obs, pred, n.par){
  R2 <- cor(obs,pred,use="complete.obs")^2
  n <- length(obs)
  1-(1-R2)*(n - 1)/(n-n.par-1)
}

# Empirical data models ----
dat <- read.table("data_TMB/ibex.txt", header=T)
dat <- read.table("data_TMB/SvreinIPM.txt", header=T); dat$ROS <- log(dat$ROS+1)
dat <- read.table("data_TMB/RedDeer.txt", header=T); dat$snow <- c(dat$snow[-1],NA)
dat <- read.table("data_TMB/SoaySheep.txt", header=T)
dat <- read.table("data_TMB/moskus.txt", header=T); dat$snowdepth_may <- log(dat$snowdepth_may+1)
dat <- read.table("data_TMB/MuleDeer.txt", header=T); dat$snowH20content_April <- c(-log(dat$snowH20content_April)[-1],NA) # inverse log-scale!

dat$r <- c(diff(log(dat$N)), NA)
dat2 <- dat[complete.cases(dat),]
N <- dat2$N
obs.r <- dat2$r
cov1 <- as.vector(scale(dat2[,!(names(dat2) %in% c("N", "year","r"))]))

# RICKER
#model = 0
# Bev Holt
#model = 1
#
PENALTY <- 0.0001 
input0 <- list(N = N, obs_r = obs.r, cov1 = cov1, penalty=PENALTY, model = model, error = 0, variable = 0)
input1 <- list(N = N, obs_r = obs.r, cov1 = cov1, penalty=PENALTY, model = model, error = 0, variable = 1)
input3 <- list(N = N, obs_r = obs.r, cov1 = cov1, penalty=PENALTY, model = model, error = 1, variable = 0)
input4 <- list(N = N, obs_r = obs.r, cov1 = cov1, penalty=PENALTY, model = model, error = 1, variable = 1)
#
parameters <- list(z=rep(0, length(obs.r)), par=c(1, log(.01), -1, 0))
# FINE-TUNE SV.REIN obj5
#parameters2 <- list(z=rep(0, length(obs.r)), par=c(0.296, 0.00098, 0.784, 0.495, 1.02))
#
obj0 <- MakeADFun(input0, parameters,random="z",DLL = "tmb_RickBH_cov2", method = "BFGS",silent = FALSE)
obj1 <- MakeADFun(input1, parameters,random="z",DLL = "tmb_RickBH_cov2", method = "BFGS",silent = FALSE)
obj3 <- MakeADFun(input3, parameters,random="z",DLL = "tmb_RickBH_cov2", method = "BFGS",silent = FALSE)
obj4 <- MakeADFun(input4, parameters,random="z",DLL = "tmb_RickBH_cov2", method = "BFGS",silent = FALSE)
#
obj0$hessian <- obj1$hessian <- obj3$hessian <-obj4$hessian  <- FALSE
obj0$control <- obj1$control <- obj3$control <-obj4$control  <- list("maxit"=10000)
opt0 <- do.call("optim", obj0)
opt1 <- do.call("optim", obj1)
opt3 <- do.call("optim", obj3)
opt4 <- do.call("optim", obj4)
opt0 <- nlminb(start=obj0$par, objective=obj0$fn, gradient=obj0$gr, control= list(iter.max=1000, eval.max= 1500))
opt1 <- nlminb(start=obj1$par, objective=obj1$fn, gradient=obj1$gr, control= list(iter.max=1000, eval.max= 1500))
opt3 <- nlminb(start=obj3$par, objective=obj3$fn, gradient=obj3$gr, control= list(iter.max=1000, eval.max= 1500))
opt4 <- nlminb(start=obj4$par, objective=obj4$fn, gradient=obj4$gr, control= list(iter.max=1000, eval.max= 1500))
AICc0 <- AICc.bp(opt0$objective, k = length(opt0$par))
AICc1 <- AICc.bp(opt1$objective, k = length(opt1$par))
AICc3 <- AICc.bp(opt3$objective, k = length(opt3$par))
AICc4 <- AICc.bp(opt4$objective, k = length(opt4$par))
#
aicc <- c(AICc0,AICc1,AICc3, AICc4)
round(aicc,2)
round(aicc-min(aicc),2)
# Ricker: add, mult, both;  BH: add, mult, both
tab0 <-rbind(summary(sdreport(obj0), "report", p.value=T))
tab1 <-rbind(summary(sdreport(obj1), "report", p.value=T))
tab3 <-rbind(summary(sdreport(obj3), "report", p.value=T))
tab4 <-rbind(summary(sdreport(obj4), "report", p.value=T))
tab0
tab1
tab3
tab4

# R2 - Ricker
# add cov
tab <- tab0
ested <- tab["a",1] - tab["b",1]*N + tab["gamma1",1]*cov1
(R2 <- cor(ested,obs.r,use="complete.obs")^2)
adjR2(ested, obs.r, length(opt0$par))
# mult cov
tab <- tab1
tab <- tab4
ested <- tab["a",1] - tab["b",1]*N* exp(tab["gamma1",1]*cov1)
(R2 <- cor(ested,obs.r,use="complete.obs")^2)
adjR2(ested, obs.r, length(opt0$par))

# R2 - BevHolt
# add cov
tab <- tab0
ested <- tab["a",1] - log(1 +tab["b",1]*N) + tab["gamma1",1]*cov1
(R2 <- cor(ested,obs.r,use="complete.obs")^2)
adjR2(ested, obs.r, length(opt0$par))
# mult cov
tab <- tab4
ested <- tab["a",1] - log(1 +tab["b",1]*N)* exp(tab["gamma1",1]*cov1)
(R2 <- cor(ested,obs.r,use="complete.obs")^2)
adjR2(ested, obs.r, length(opt0$par))
