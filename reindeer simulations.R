
# top ----
library(arm)
library(scales)
library(mvtnorm)

### LOAD DATA ###
# Full age structure during 1994-2014
load("reindeer data/AgeStruc.rdata") 
# AgeStr.BP = age structure (0-13+ yr olds)  estimated based on cohort analyses of the final IPM with 6 age classes,
#      the last number is the total population size.

# ROS simulations
load("reindeer data/ROS.proj.rdata") 
# ~ 5 scenarios (very low, low, medium, high, very high frequencies of ROS events)
#     for Peeters et al. 2022 Ecol Lett, we used scenarios 1, 3, 5
# ~ for each scenario, there are 6000 simulated trajectories of 100 time steps

### model parameters from Hansen et al. 2019 Nat. Commun. for 9090 posterior models from the IPM
load(file="reindeer data/coeffic.Rdata") # survival
load(file="reindeer data/coefficF.Rdata") # fecundity
colnames(coeffic) <- c("int","A2","A3","A4","A5","A6", "N","year","ROS","win","A2:ROS","A3:ROS","A4:ROS","A5:ROS","A6:ROS","k")
colnames(coefficF) <- c("int","A3","A4","A5","A6", "N","year","ROS","win","A3:ROS","A4:ROS","A5:ROS","A6:ROS","k")
# residual covariance matrices
load(file="reindeer data/SigmaInt.Rdata")

# scale N -- estimated effect of population density on survival and fecundity is for standardized population sizes
meanN <- mean(AgeStr.BP[-21,15])
sdN <- sd(AgeStr.BP[-21,15])
# mean winter length during 1994-2014 (study period for the)
mean.winlen = 240.6


# Simulations using mean parameter estimates from survival and fecundity models
coefS <- colMeans(coeffic)
coefF <- colMeans(coefficF)
meanSigmaInt <- apply(SigmaInt, c(1,2),mean)
input.data.mean <- list(coefS2=coefS, coefF2 = coefF, SigmaInt = meanSigmaInt,
                             ROS.proj = ROS.proj[c(1,3,5),,], winlen = mean.winlen,
                             age.str = AgeStr.BP[21,], meanN=meanN, sdN=sdN)
### Function to simulate reindeer trajectories with harvesting
### simSV
simSV <- function(h.par, h.type="prop", nsim=10, ti=100,input.data.=input.data.mean,seed=50){
  require(arm)
  require(mvtnorm)
  
  ### h.ty should be a vector describing harvest types:
  # "n" = no hunt, "c" = constant, "p" = proportion, "t" = threshold, "pt" = proportional threshold
  ### h.par should be matrix with 2 columns: 
  # ~ first column (n) indicates number of animals (for constant or threshold)
  # ~ second column (p) indicates proportion of N harvested
  ### input.data = list with model output for survival (reg.S) and fecundity (reg.F), 
  # ~ ROS.proj = ROS projections : array with dimensions m=5, s=6000, t=100 (m: scenarios very low - very high frequency (1-5), s = 6000 simulations, t= 100 timesteps)
  # ~ age.str = age structure column 1:14 for first year (2014) and total N in column 15
  
  # -- prepare arrays #
  fxS <- input.data.$coefS2
  fxF <- input.data.$coefF2
  
  nscen = dim(input.data.$ROS.proj)[1] # nbr of ROS scenarios
  nhunt = length(h.par)
  
  fec.mat <- surv.mat <- c()
  N.sim <-  H.sim <- Posthunt <- array(NA,c(nscen,nhunt, nsim, 14,ti+1))
  Ntot.sim<- Htot.sim <- Posthunt.tot <- array(NA, c(nscen, nhunt, nsim, ti+1))
  r.sim <- lambda.sim <- array(NA,c(nscen,nhunt, nsim,ti))
  
  #-- model parameters #
  # N.sim first year = 2014
  for(i in 1:14){
    N.sim[,,,i,1] <- input.data.$age.str[i] # only need last row, start from 2014
  }
  Ntot.sim[,,,1] <- sum(input.data.$age.str[1:14])
  
  # + SIMULATE  -
  nr.sim=0; ntot=prod(nscen,nhunt,nsim) # for tracking progress
  
  # use righ order of loops: first posterior sample, then ROS scenario, then subsample ROS, then loop ROS simulation, and then harvest
  
  for(m in 1:nscen){ # nber of ROS scenarios
    for(s in 1:nsim)  { # nber of ROS simulations
      
      res.sim <- rmvnorm(ti, sigma=meanSigmaInt)
      
      for(h in 1:nhunt){ # use h for changing singular harvest parameters
        for (t in 1:ti){
          
          # 1) post-hunt N
          # assume same age-structure in harvest offtake as in life population
          age.dist <- N.sim[m,h,s,,t]/sum(N.sim[m,h,s,,t])
          
          # Proportional harvest
          if(h.type=="prop"){
            n.shot <- sum(N.sim[m,h,s,,t])*h.par[h]
            H.sim[m,h,s,,t] <- round(n.shot*age.dist)
          }
          if(h.type=="cte"){ # cte = constant
            H.sim[m,h,s,,t] <- round(h.par[h]*age.dist)
          }
          
          
          # total yearly harvest
          Htot.sim[m,h,s,t] <- sum(H.sim[m,h,s,,t])
          # 0 if posthunt results in < 0 (e.g. constant harvest)
          Posthunt[m,h,s,,t] <- N.sim[m,h,s,,t] - H.sim[m,h,s,,t]
          Posthunt[m,h,s,,t] <-  ifelse(Posthunt[m,h,s,,t] <0, 0, Posthunt[m,h,s,,t])
          Posthunt.tot[m,h,s,t] <- sum(Posthunt[m,h,s,,t])
          
          # in case population extinct:
          if(sum(Posthunt[m,h,s,,t],na.rm=T) == 0){
            N.sim[m,h,s,,t+1] <- 0
            Ntot.sim[m,h,s,t+1] <- 0
            Posthunt[m,h,s,,t+1] <- 0
            Posthunt.tot[m,h,s,t+1] <- 0
            next # don't need to calculate any vital rates anymore
          }
          
          Ntotal <- Posthunt.tot[m,h,s,t]
          Ntot.sc <- (Ntotal - meanN)/sdN
          
          ROS_s <- input.data.$ROS.proj[m,s,t]*exp(unname(coefS["k"])*Ntot.sc)
          ROS_f <- input.data.$ROS.proj[m,s,t]*exp(unname(coefF["k"])*Ntot.sc)
          
          ### 2) Vital rates
          # Age-specific survival
          surv.mat[1]     <-invlogit(fxS["int"] + fxS["ROS"]*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,1])
          surv.mat[2]     <-invlogit(fxS["A2"] + fxS["int"] + (fxS["ROS"] + fxS["A2:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,2])
          surv.mat[3]     <-invlogit(fxS["A3"] + fxS["int"] + (fxS["ROS"] + fxS["A3:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,3])
          surv.mat[4:9]   <-invlogit(fxS["A4"] + fxS["int"] + (fxS["ROS"] + fxS["A4:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,4])
          surv.mat[10:12] <-invlogit(fxS["A5"] + fxS["int"] + (fxS["ROS"] + fxS["A5:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,5])
          surv.mat[13]    <-invlogit(fxS["A6"] + fxS["int"] + (fxS["ROS"] + fxS["A6:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,6])
          
          
          # Age-specific fecundity
          fec.mat[1]     <-NA
          fec.mat[2]     <-invlogit(fxF["int"] + fxF["ROS"]*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,7])
          fec.mat[3]     <-invlogit(fxF["A3"] + fxF["int"] + (fxF["ROS"] + fxF["A3:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,8])
          fec.mat[4:9]   <-invlogit(fxF["A4"] + fxF["int"] + (fxF["ROS"] + fxF["A4:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,9])
          fec.mat[10:12] <-invlogit(fxF["A5"] + fxF["int"] + (fxF["ROS"] + fxF["A5:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,10])
          fec.mat[13]    <-invlogit(fxF["A6"] + fxF["int"] + (fxF["ROS"] + fxF["A6:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,11])
          
          
          ### 3) Population development
          # survival: survivors from hunt have to go through winter
          for(a in 2:13){
            N.sim[m,h,s,a,t+1] <- rbinom(1,Posthunt[m,h,s,a-1,t], surv.mat[a-1])
          } 
          N.sim[m,h,s,14,t+1] <- rbinom(1,Posthunt[m,h,s,13,t] + Posthunt[m,h,s,14,t], surv.mat[13])
          # calves: only next year's survivors can give birth
          calves <- 0
          for(a in 2:13){
            calves[a] <- rbinom(1, N.sim[m,h,s,a+1,t+1], fec.mat[a]) 
          }
          N.sim[m,h,s,1,t+1] <- rbinom(1,sum(calves), 0.5) # Sex ratio = 0.5
          # total pop
          Ntot.sim[m,h,s,t+1] <- sum(N.sim[m,h,s,,t+1]) # pre-hunt pop size
          
          # growth rate
          lambda.sim[m,h,s,t] <- ifelse(Ntot.sim[m,h,s,t] == 0, NA,
                                        Ntot.sim[m,h,s,t+1]/Ntot.sim[m,h,s,t]) # if Nt = 0, growth rate = NA
          r.sim[m,h,s,t] <- ifelse(lambda.sim[m,h,s,t]==0, NA,log(lambda.sim[m,h,s,t]))
          
        }#t
        
        # progress
        nr.sim = nr.sim+1
        if(nr.sim==1 | floor(nr.sim/ntot*10)>floor((nr.sim-1)/ntot*10)){
          cat(round(nr.sim/ntot*100), "%  ", format(Sys.time(), "%H:%M:%S"),"\n",sep="")
        } # progress
      }#h
    }#s
  }#m
  # ~ gather in list
  #sim3 <- list(dimensions=c(scen=nscen, h=nhunt, posterior=npost, sim=nsim, a= 14, t=ti+1), h.par=h.par, h.ty=h.ty, Ntot.sim=Ntot.sim, 
  #             N.sim=N.sim, H.sim=H.sim, Htot.sim=Htot.sim, Posthunt.tot = Posthunt.tot, Posthunt = Posthunt, lambda.sim = lambda.sim, r.sim=r.sim)
  sim3 <- list(dimensions=c(scen=nscen, h=nhunt, sim=nsim, a= 14, t=ti+1), 
               h.par=h.par,  N.sim=N.sim, Ntot.sim=Ntot.sim, r.sim=r.sim)
  return(sim3)
}

# example
h.par.traj=c(0,0.15)
opt.traj <- simSV(nsim=3, ti=100, h.par=h.par.traj, input.data.=input.data.mean)
opt.traj$dim
ts.plot(t(opt.traj$Ntot.sim[,1,1,]),col=c(1,2,4)) # 1 simulation of each ROS scenario, no harvest
ts.plot(t(opt.traj$Ntot.sim[1,,1,]),col=c(1,2)) # 1 simulation of each harvest proportion, very low ROS frequency


# Simulations using samples from the 9090 posterior models of survival and fecundity 
input.data.posterior <- list(coefS2=coeffic, coefF2 = coefficF, SigmaInt = SigmaInt,
                             ROS.proj = ROS.proj, winlen = mean.winlen,
                             age.str = AgeStr.BP[21,], meanN=meanN, sdN=sdN)
### simSV_posterior
simSV_posterior <- function(h.par, h.type="prop", scen=c(1,3,5), nsim=10, npost=3, ti=100, input.data.=input.data.posterior, seed=50){
  require(arm)
  require(mvtnorm)
  
  # h.type = either "prop" for proportional harvesting or "cte" for constant harvesting
  # h.par =  vector of at least size 2; 
  # ~ for "prop" -- proportion of annually harvested reindeer 
  # ~ for "cte" -- total number of reindeer annually harvested
  # scen = which ROS scenarios to select: 1-5 are very low, low, medium, high, and very high frequencies of ROS events
  
  # npost = number of posterior models to use 
  # nsim = number of simulated ROS trajectories, randomly drawn for each npost and ROS scenario
  # ti = number of timesteps in each trajectory
  
  # -- sample from posterior estimates
  set.seed(seed)
  npost.sample = sample(1:9090,npost)
  fixefS <- input.data.$coefS2[npost.sample,]
  fixefF <- input.data.$coefF2[npost.sample,]
  postSigmaInt <- input.data.$SigmaInt[,,npost.sample]
  # -- prepare arrays #
  nscen = length(scen) # no. of ROS scenarios
  nhunt = length(h.par) # no. of harvest scenarios
  fec.mat <- surv.mat <- c()
  N.sim <-  H.sim <- Posthunt <- array(NA,c(npost,nscen,nhunt,nsim, 14,ti+1))
  Ntot.sim<- Htot.sim <- Posthunt.tot <- array(NA, c(npost,nscen, nhunt, nsim, ti+1))
  r.sim <- lambda.sim <- array(NA,c(npost,nscen,nhunt, nsim,ti))
  
  #-- initial age structure and pop size #
  # N.sim first year = 2014
  for(i in 1:14){
    N.sim[,,,,i,1] <- input.data.$age.str[i] # only need last row, start from 2014
  }
  Ntot.sim[,,,,1] <- sum(input.data.$age.str[1:14])
  
  # + SIMULATE  -
  nr.sim=0; ntot=prod(nscen,npost,nhunt,nsim) # for tracking progress
  
  # use right order of loops: first posterior sample, then ROS scenario, then subsample ROS, then loop ROS simulation, and then harvest
  for(p in 1:npost){
    fxS <- fixefS[p,]
    fxF <- fixefF[p,]
    set.seed(seed)
    res.sim <- rmvnorm(ti, sigma=postSigmaInt[,,p]) # residuals from multivariate normal distribution 
    
    # sample ROS simulations anew for each posterior sample
    set.seed(seed + p)
    ROS.sel <- input.data.$ROS.proj[scen,sample(1:6000, nsim),]
    
    for(m in 1:nscen){ # no. of ROS scenarios
      for(s in 1:nsim)  { # no. of simulated trajectories per ROS scenario and posterior model
        
        for(h in 1:nhunt){ # use h for changing singular harvest parameters
          for (t in 1:ti){
            
            # 1) post-hunt N
            # assume same age-structure in harvest offtake as in life population
            age.dist <- N.sim[p,m,h,s,,t]/sum(N.sim[p,m,h,s,,t])
            
            # Proportional harvest
            if(h.type=="prop"){
              n.shot <- sum(N.sim[p,m,h,s,,t])*h.par[h]
              H.sim[p,m,h,s,,t] <- round(n.shot*age.dist)
            }
            if(h.type=="cte"){ # cte = constant
              H.sim[p,m,h,s,,t] <- round(h.par[h]*age.dist)
            }
            
            # total yearly harvest
            Htot.sim[p,m,h,s,t] <- sum(H.sim[p,m,h,s,,t])
            # 0 if posthunt results in < 0 (e.g. constant harvest)
            Posthunt[p,m,h,s,,t] <- N.sim[p,m,h,s,,t] - H.sim[p,m,h,s,,t]
            Posthunt[p,m,h,s,,t] <-  ifelse(Posthunt[p,m,h,s,,t] <0, 0, Posthunt[p,m,h,s,,t])
            Posthunt.tot[p,m,h,s,t] <- sum(Posthunt[p,m,h,s,,t])
            
            # in case population extinct:
            if(sum(Posthunt[p,m,h,s,,t],na.rm=T) == 0){
              N.sim[p,m,h,s,,t+1] <- 0
              Ntot.sim[p,m,h,s,t+1] <- 0
              Posthunt[p,m,h,s,,t+1] <- 0
              Posthunt.tot[p,m,h,s,t+1] <- 0
              next # don't need to calculate any vital rates anymore
            }
            
            Ntotal <- Posthunt.tot[p,m,h,s,t]
            Ntot.sc <- (Ntotal - meanN)/sdN
            
            ROS_s <- ROS.sel[m,s,t]*exp(unname(fxS["k"])*Ntot.sc)
            ROS_f <- ROS.sel[m,s,t]*exp(unname(fxF["k"])*Ntot.sc)
            
            ### 2) Vital rates
            # Age-specific survival
            surv.mat[1]     <-invlogit(fxS["int"] + fxS["ROS"]*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,1])
            surv.mat[2]     <-invlogit(fxS["A2"] + fxS["int"] + (fxS["ROS"] + fxS["A2:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,2])
            surv.mat[3]     <-invlogit(fxS["A3"] + fxS["int"] + (fxS["ROS"] + fxS["A3:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,3])
            surv.mat[4:9]   <-invlogit(fxS["A4"] + fxS["int"] + (fxS["ROS"] + fxS["A4:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,4])
            surv.mat[10:12] <-invlogit(fxS["A5"] + fxS["int"] + (fxS["ROS"] + fxS["A5:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,5])
            surv.mat[13]    <-invlogit(fxS["A6"] + fxS["int"] + (fxS["ROS"] + fxS["A6:ROS"])*ROS_s + fxS["N"]*Ntot.sc + fxS["year"]*2014 + fxS["win"]*input.data.$winlen + res.sim[t,6])
            
            
            # Age-specific fecundity
            fec.mat[1]     <-NA
            fec.mat[2]     <-invlogit(fxF["int"] + fxF["ROS"]*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,7])
            fec.mat[3]     <-invlogit(fxF["A3"] + fxF["int"] + (fxF["ROS"] + fxF["A3:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,8])
            fec.mat[4:9]   <-invlogit(fxF["A4"] + fxF["int"] + (fxF["ROS"] + fxF["A4:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,9])
            fec.mat[10:12] <-invlogit(fxF["A5"] + fxF["int"] + (fxF["ROS"] + fxF["A5:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,10])
            fec.mat[13]    <-invlogit(fxF["A6"] + fxF["int"] + (fxF["ROS"] + fxF["A6:ROS"])*ROS_f + fxF["N"]*Ntot.sc + fxF["year"]*2014 + fxF["win"]*input.data.$winlen + res.sim[t,11])
            
            
            ### 3) Population development
            # survival: survivors from hunt have to go through winter
            for(a in 2:13){
              N.sim[p,m,h,s,a,t+1] <- rbinom(1,Posthunt[p,m,h,s,a-1,t], surv.mat[a-1])
            } 
            N.sim[p,m,h,s,14,t+1] <- rbinom(1,Posthunt[p,m,h,s,13,t] + Posthunt[p,m,h,s,14,t], surv.mat[13])
            # calves: only next year's survivors can give birth
            calves <- 0
            for(a in 2:13){
              calves[a] <- rbinom(1, N.sim[p,m,h,s,a+1,t+1], fec.mat[a]) 
            }
            N.sim[p,m,h,s,1,t+1] <- rbinom(1,sum(calves), 0.5) # Sex ratio = 0.5
            # total pop
            Ntot.sim[p,m,h,s,t+1] <- sum(N.sim[p,m,h,s,,t+1]) # pre-hunt pop size
            
            # growth rate
            # if Nt = 0, growth rate = NA
            lambda.sim[p,m,h,s,t] <- ifelse(Ntot.sim[p,m,h,s,t] == 0, NA,
                                            Ntot.sim[p,m,h,s,t+1]/Ntot.sim[p,m,h,s,t]) 
            r.sim[p,m,h,s,t] <- ifelse(lambda.sim[p,m,h,s,t]==0, NA,
                                       log(lambda.sim[p,m,h,s,t]))
          }#t
          
          # progress
          nr.sim = nr.sim+1
          if(nr.sim==1 | floor(nr.sim/ntot*20)>floor((nr.sim-1)/ntot*20)){
            cat(round(nr.sim/ntot*100), "%  ", format(Sys.time(), "%H:%M:%S"),"\n",sep="")
          }
        }#h
      }#s
    }#m
  }#p
  
  # ~ gather in list
  #sim3 <- list(dimensions=c(scen=nscen, h=nhunt, posterior=npost, sim=nsim, a= 14, t=ti+1), h.par=h.par, h.ty=h.ty, Ntot.sim=Ntot.sim, 
  #             N.sim=N.sim, H.sim=H.sim, Htot.sim=Htot.sim, Posthunt.tot = Posthunt.tot, Posthunt = Posthunt, lambda.sim = lambda.sim, r.sim=r.sim)
  sim3 <- list(dimensions=c(post=npost, scen=nscen, h=nhunt, sim=nsim, a= 14, t=ti+1), 
               h.par=h.par,h.type=h.type, 
               N.sim=N.sim, Ntot.sim=Ntot.sim, r.sim=r.sim)
  return(sim3)
}

### Example proportional harvesting
h.par <- seq(0,0.3,by=0.1) # have to use a vector of at least size 2
# simulate
test.sim <- simSV_posterior(h.par=h.par, npost=5, nsim=10, ti=100)
#
test.sim$dimensions #
# post = no. posterior models used
# scen = no. ROS scenarios
# h = no. of harvest scenarios (proportions/constant harvest rates)
# sim = no. of simulated trajectories for each post, scen and h
# a = no. of age classes
# t = number of time steps in each trajectory, incl. the initial population structure

str(test.sim)
# h.par = harvest parameters (proportions or no. of harvested reindeer per year)
# h.type = harvest type, "prop" or "cte"
# N.sim = array of annual pop sizes in each age class (14), for each simulation
# Ntot.sim = array of total annual pop size pre-harvesting
# Posthunt.tot = array of total annual population size after harvesting
# r.sim = array of annual logistic growth rates (i.e. difference in log population sizes between years)

### Example constant harvesting
#h.par.cte <- c(0,20,40,60)
#test.sim2 <- simSV_posterior(h.par=h.par.cte,h.type="cte", npost=5, nsim=10, ti=100)

# check quantiles
apply(apply(test.sim$Ntot.sim, c(2,3,1,4),mean),c(1,2),quantile,c(0.025,0.975))

### Calculate population properties ###
# variation in r
var.r <- apply(apply(test.sim$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)
var.r
# rows = ROS scenarios, columns = harvest scenarios

# mean N
N_mean <- apply(test.sim$Ntot.sim, c(2,3),mean)
N_mean

# to calculate probability of population crash
# crash was defined as population reduction by half, i.e. lambda < 0.5, or r < log(0.5)
crash <- apply(apply(test.sim$r.sim, c(2,3,1,4), function(x) any(x< log(0.5))),c(1,2),function(x) sum(x)/prod(test.sim$dim[c("post","sim")]))
crash <- ifelse(is.na(crash),0,crash)
crash

# to calculate probability of quasi-extinction
round(AgeStr.BP[21,15]*.20) # 20% of initital pop size in 2014 as threshold for quasi-extinction
quasi.N <- 350
test.sim$dim
Pquasi <- apply(apply(test.sim$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),c(1,2),function(x)sum(x)/prod(test.sim$dim[c("post","sim")]))
Pquasi # rows = Ros scenarios, columns = harvest proportions

