

# functions ----

# harvest function
# for constant (cte) harvest, set constant yield as a proportion of K
h <- function(N, p, K, type = "prop"){
  if(type=="prop") X <- N-p*N
  if(type=="cte") X <- ifelse(N - K*p <=0, 0, N - K*p)
  X
}
# growth rate
R.fn <- function(N, K, a, c=0, d=0, z1=0, z2=0, model){
  if(model == "Ri"){ # Ricker 
    b <- a / K
    R <- a + c*z1 - b*N*exp(-d*z2)
  }
  if(model == "BH"){ # Beverton Holt
    b <- (exp(a)-1) / K
    R <- a + c*z1 - log(1 + b*N)*exp(-d*z2)
  }
  R
} 
# pop size next time step
N.fn <- function(N, K=100, a=1, c=0, d=0, z1=0, z2=0, model="Ri", p=0, ...){
  hN <- h(N, p, K, ...) # harvest
  hN*exp(R.fn(N=hN, K, a, c, d, z1, z2, model))
}

# time series with nt time steps and correlated additive/multiplicative noise
ts.fn <- function(K=100, a=1, c=0, d=0, rho=1, model="Ri", p=0, nt=100,...){
  require(mvtnorm)
  z12 <- rmvnorm(nt, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
  N <- K
  for(t in 1:nt) N[t+1] <- N.fn(N[t], K, a, c, d, z1=z12[t,1], z2=z12[t,2], model, p, ...)
  N
}

# examples
# harvest effect on growth
plot(log(N.fn(1:200,  K=100, a=1, model="Ri", p=0)) - log(1:200),type="l")
lines(R.fn(N=1:200, K=100, a=1, model="Ri"), col=2, lty=2)
lines(log(N.fn(1:200,  K=100, a=1, model="Ri", p=0.2)) - log(1:200),col=2)


# SIMULATIONS ####

# Ricker vs. Beverton Holt ####
# comparing variance in log N for Ricker vs. BH ~ Figure 3 main paper 
# to run simulation on multiple cores
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
library(foreach)

#### hyperparameters
K=100
pp = seq(0,.3,le=4)  # harvest proportions for low beta0
pp2 = seq(0,0.5,le=4) # for medium beta0
pp3 = seq(0,0.7,le=4)  # for high beta0
aa=c(.5,1,1.5)
dd=c(0.39715, 0.2159, 0.14668)
# about the choice of parameters aa and dd:
# a = beta0, low to high max growth rate.
# d, i.e. effect of multiplicative variance, chosen so that variance in R = 0.05, same as additive variance
a = .5; d = 0.39715
var(R.fn(N=100, K=100, a=a, c=0, d=d, z1=0, z2=rnorm(10000000,0,1), model="BH"))
var(R.fn(N=100, K=100, a=a, c=0, d=d, z1=0, z2=rnorm(10000000,0,1), model="Ri"))
#a = 1; d = 0.2159; p = 0.25
#a = 1.5; d = 0.14668; p = 0.4

# hyper
hyperparam <- expand.grid(p=pp,a=aa, model=c("Ri","BH"))
hyperparam$d <- ifelse(hyperparam$a == aa[1], dd[1],
                       ifelse(hyperparam$a == aa[2], dd[2], dd[3]))
hyperparam[hyperparam$a==aa[2],]$p <- round(pp2,2)
hyperparam[hyperparam$a==aa[3],]$p <- round(pp3,2)
hyperparam$c <- sqrt(0.05)
hyperparam$K <- K
hyperparam$type <- "prop"
hyperparam1 <- hyperparam2 <- hyperparam
hyperparam1$c <- 0; hyperparam2$d <- 0
hyperparam1$noise <- "mult"; hyperparam2$noise <- "add"
hyperparam <- rbind(hyperparam1,hyperparam2)
str(hyperparam)
# simulations
nt=10000
bean.sim <- foreach(
  p=hyperparam$p, 
  a=hyperparam$a,
  d=hyperparam$d,
  c=hyperparam$c,
  K=hyperparam$K,
  model=hyperparam$model,
  type=hyperparam$type
) %dopar% {
  ts.fn(K=K, a=a, c=c, d=d, model=model, nt=nt, p=p)
}
str(bean.sim)
str(hyperparam)
nt = length(bean.sim[[1]])

# EXAMPLE BEANPLOT
library(beanplot)
hyperparam$noise <- as.factor(hyperparam$noise)
hyperparam$p.noise <- as.factor(paste(hyperparam$p,hyperparam$noise))
nt = length(bean.sim[[1]])
par.a1.ri <- hyperparam$p.noise[which(hyperparam$a == aa[1] & hyperparam$model == "Ri")]
sim.a1.ri <- unlist(bean.sim[which(hyperparam$a == aa[1] & hyperparam$model == "Ri")]) 
sub.a1.ri <- data.frame(sim = sim.a1.ri, mod = rep(par.a1.ri, each=nt))
sub.a1.ri$mod <- factor(sub.a1.ri$mod)
head(sub.a1.ri)

range(sub.a1.ri$sim)
ylim <- quantile(sub.a1.ri$sim, c(0.001,0.9999))
beanplot(sim ~ mod, data = sub.a1.ri, ll = 0,side="both",log="y",
        ylab = "N",xlab="Harvest proportion",
         maxwidth=1, what=c(0,1,1,1), beanlinewd=1, ylim=ylim,
         border = NA, col = list("black", c("grey", "white")))
abline(h=100,lty=3,col="grey")
legend("bottomleft", fill = c("black", "grey"),bty="n",
       legend = c("Additive", "Multiplicative"),title="Stochasticity")





# Quasi-extinction probabilities ####

# quasi-extinction probability for a range of beta0 and variance in r
# multiplicative, additive and mixed env. noise models (Ricker only)

h <- function(N, p) N - p*N
N.rma <-  function(a,b,c,d, N, p, z1, z2) h(N,p)*exp(a - h(N,p)*exp(b + c*z1) + d*z2)
VAR.rma <- function(a,b,c,d,rho, N=K) N^2 * exp(2*b+c^2)*(exp(c^2)-1) + d^2 + N*d*c*rho*exp(b+0.5*c^2)
# find parameters for a given variance in R for N = K
fc <- function(c, d, varR)  K^2 * exp(2*b+c^2)*(exp(c^2)-1) + d^2 + K*d*c*rho*exp(b+0.5*c^2) - varR

# simulate timeseries stopping when N below threshold
fn.ext.rma <- function(a,c,d,rho,P,ns=1000,nt=100,K=50, tres = K/5 ){
  require(mvtnorm)
  b = log(a/K)
  S <- c()
  
  N.sim <- matrix(NA, nr=ns, nc=nt)
  N.sim[,1] <- K
  
  for(n in 1:ns){
    t = 1
    z12 <- rmvnorm(nt, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
    while(N.sim[n,t] > tres & t < nt){
      N.sim[n,t+1] = N.rma(a,b,c,d,  N.sim[n,t], p=P, z1=z12[t,1], z2=z12[t,2])
      t = t+1
    }
  }
  return(N.sim)
}

#..............................................................................
### Multiplicative Variance
K=50
d= 0
rho=0

aa <- seq(0.2,2,by=0.1)
varR = seq(0.01,0.1,le=length(aa))
bb <- log(aa/K)
tres=K/5 # quasi-extinction threshold below 20% of K
ns= 10 # number of simulations per parameter set
nt= 1000 # number of time steps
arSim0 <- arSim1 <-arSim2 <- opt.cc <- array(NA, dim=c(length(aa),length(varR)))

pp1 = 0.1
pp2 = 0.2

for(i in 1:length(aa)){
  a=aa[i]
  b=bb[i]
  for(j in 1:length(varR)){
    c <- opt.cc[i,j] <- uniroot(fc, d=d,varR=varR[j], lower=0.001, upper=10)$root
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=a,c=c, d=d, rho=rho, P=0,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim0[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=d, rho=rho, P=pp1,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim1[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=d, rho=rho, P=pp2,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim2[i,j] <- N.ext
  }
}
params <- list(aa=aa,bb=bb,K=K, rho=rho, d=d,varR=varR, c=opt.c, tres=tres, nt=nt, ns=ns)
simlist <- list(params=params, arSim0=arSim0, arSim1=arSim1, arSim2=arSim2,opt.c=opt.c)

#..............................................................................
### Additive Variance
K=50
c=0
rho=0

aa <- seq(0.2,2,by=0.1)
varR = seq(0.01,0.1,le=length(aa))
bb <- log(aa/K)
tres=K/5
ns= 10
nt= 1000
arSim0 <- arSim1 <-arSim2 <- opt.dd <- array(NA, dim=c(length(aa),length(varR)))

pp1 = 0.1
pp2 = 0.2

for(i in 1:length(aa)){
  a=aa[i]
  b=bb[i]
  for(j in 1:length(varR)){
    d <- opt.dd[i,j] <- uniroot(fc, c=c,varR=varR[j], lower=0.0001, upper=10)$root
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=a,c=c, d=d, rho=rho, P=0,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim0[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=d, rho=rho, P=pp1,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim1[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=d, rho=rho, P=pp2,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim2[i,j] <- N.ext
  }
}

arSim2 <- simlist$arSim2
params <- list(aa=aa,bb=bb,K=K, rho=rho, d=opt.dd,varR=varR, c=c, tres=tres, nt=nt, ns=ns)
simlist <- list(params=params, arSim0=arSim0, arSim1=arSim1, arSim2=arSim2)

#..............................................................................
### Mixed Variance

fc <- function(c, d, varR)  K^2 * exp(2*b+c^2)*(exp(c^2)-1) + d^2 + K*d*c*rho*exp(b+0.5*c^2) - varR
# For d^2 = 40% of total Variance and rho = 0.5, ~ 40% of Var is due to multipl.var and 20% due to the correlation

K=50
rho=0.25
aa <- seq(0.2,2,by=0.1)
varR = seq(0.01,0.1,le=length(aa)) # total variance in R
dd = sqrt(varR*.3) # so that additive variance contributes to 30% of the total variance

bb <- log(aa/K)
tres=K/5
ns= 10
nt= 1000
arSim0 <- arSim1 <-arSim2 <- opt.cc <- array(NA, dim=c(length(aa),length(varR)))

pp1 = 0.1
pp2 = 0.2

for(i in 1:length(aa)){
  a=aa[i]
  b=bb[i]
  for(j in 1:length(varR)){
    c <- opt.cc[i,j] <- uniroot(fc, d=dd[j],varR=varR[j], lower=0.001, upper=10)$root
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=a,c=c, d=dd[j], rho=rho, P=0,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim0[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=dd[j], rho=rho, P=pp1,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim1[i,j] <- N.ext
    
    set.seed(1)
    sim.tmp <- fn.ext.rma(a=aa[i],c=c,d=dd[j], rho=rho, P=pp2,ns=ns,nt=nt,tres=tres)
    N.ext <- sum(is.na(sim.tmp[,nt]))
    arSim2[i,j] <- N.ext
  }
}

params <- list(description="d^2 composed 30% of total variance. With rho=0.25, this gave approx. 60%                
               variance in the multiplicative part, and 10% from the covariance part in: K^2 * exp(2*b+c^2)*(exp(c^2)-1) + d^2 + K*d*c*rho*exp(b+0.5*c^2)",
               aa=aa,bb=bb,K=K, rho=rho, d=dd,varR=varR, c=opt.cc, tres=tres, nt=nt, ns=ns)
simlist <- list(params=params, arSim0=arSim0, arSim1=arSim1, arSim2=arSim2, opt.cc=opt.cc)
