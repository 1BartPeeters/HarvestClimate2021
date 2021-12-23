
# top ----
library(arm)
library(scales)
library(mvtnorm)

### LOAD DATA ###
### climate parameters and age structure 2014
load("reindeer data/simSV_inputData.rdata") 
# AgeStr.BP = age structure (0-13+ yr olds) from 1994-2014, estimated based on cohort analyses of the final IPM with 6 age classes,
#      the last number is the total population size.
# mean.winlen = mean winter length during the study period for the IPM (1994-2014)
# ROS.proj¨= ROS simulations. 
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

### input.data = list with 
# - model output for survival (coeffic) and fecundity (coefficF), 
# ~ ROS.proj = ROS projections : array with dimensions m=5, s=6000, t=100 (m: scenarios very low to very high frequency (1-5), s = 6000 simulations, t= 100 timesteps)
# - winlen = winter length, which was kept constant for all simulations at the mean winter length during 1994-2014,
# ~ age.str = age structure column 1:14 for first year (2014) and total N in column 15, used to initiate the simulations

# for simulations using mean survival and fecundity model estimates
coefS <- colMeans(coeffic)
coefF <- colMeans(coefficF)
meanSigmaInt <- apply(SigmaInt, c(1,2),mean)
input.data.mean <- list(coefS2=coefS, coefF2 = coefF, SigmaInt = meanSigmaInt,
                             ROS.proj = ROS.proj[c(1,3,5),,], winlen = mean.winlen,
                             age.str = AgeStr.BP[21,], meanN=meanN, sdN=sdN)

# for simulations using posterior models
input.data.posterior <- list(coefS2=coeffic, coefF2 = coefficF, SigmaInt = SigmaInt,
                             ROS.proj = ROS.proj, winlen = mean.winlen,
                             age.str = AgeStr.BP[21,], meanN=meanN, sdN=sdN)

### Functions to simulate reindeer trajectories with harvesting
# ~ simSV uses the mean estimates of the posterior survival and fecundity models
# ~ simSV_posterior uses npost random draws from the 9090 posterior models for each nsim randomly drawn ROS trajectories

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

# Example of how to calculate properties: quantiles of total pop size
apply(apply(test.sim$Ntot.sim, c(2,3,1,4),mean),c(1,2),quantile,c(0.025,0.975))

### Example constant harvesting
h.par.cte <- c(0,20,40,60)
test.sim2 <- simSV_posterior(h.par=h.par.cte,h.type="cte", npost=5, nsim=10, ti=100)
# check quantiles
apply(apply(test.sim2$Ntot.sim, c(2,3,1,4),mean),c(1,2),quantile,c(0.025,0.975))


### Calculate population properties ###
# variation in r
var.r <- apply(apply(test.sim$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)

# mean N
N_mean <- apply(test.sim$Ntot.sim, c(2,3),mean)

# to calculate probability of population crash
crash <- apply(apply(test.sim$lambda.sim, c(2,3,1,4), function(x) any(x< 0.5)),c(1,2),function(x) sum(x)/prod(test.sim$dim[c("post","sim")]))
crash <- ifelse(is.na(crash),0,crash)

# to calculate probability of quasi-extinction
round(AgeStr.BP[21,15]*.20) # 20% of initital pop size in 2014 as threshold for quasi-extinction
quasi.N <- 350
test.sim$dim
Pquasi <- apply(apply(test.sim$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),c(1,2),function(x)sum(x)/prod(test.sim$dim[c("post","sim")]))
Pquasi # rows = Ros scenarios, columns = harvest proportions


#________________________________________________________________________________________________________####
##FULL SIMULATION ####

## NOT RUN ##
#nROS=10 # nr ROS simulations per posterior model
#nPost=1000 # nr posterior models

### proportional harvest 
## 1st simulation
#h.par <- seq(0,0.2,by=0.01)
#FULL.SIM <- simSV_posterior(h.par=h.par, npost=nPost, nsim=nROS, ti=100)
#str(FULL.SIM)
#FULL.SIM$dim <- FULL.SIM$dim[c("post","scen","h","sim","a","t")]
#save(FULL.SIM, file="reindeer data/simulations/FULL_SIM_posterior_20211007.Rdata")

## 2nd simulation
#h.par.p2 <- seq(0.21,0.30,by=0.01)
#FULL.SIM.prop2 <- simSV_posterior(h.par=h.par.p2, npost=nPost, nsim=nROS, ti=100)
#FULL.SIM.prop2$dim <- FULL.SIM.prop2$dim[c("post","scen","h","sim","a","t")]
#save(FULL.SIM.prop2, file="reindeer data/simulations/FULL_SIM_2_posterior_20211007.Rdata")

### constant harvest
## 1st simulation
# harvest proportions
#h.par <- seq(0,200,by=10)
#h.par <- seq(210,300,by=10)
#FULL.SIM.cte <- simSV_posterior(h.par=h.par,h.type="cte", npost=nPost, nsim=nROS, ti=100)
#FULL.SIM.cte$dim <- FULL.SIM.cte2$dim[c("post","scen","h","sim","a","t")]
#save(FULL.SIM.cte, file="reindeer data/simulations/FULL_SIM_CONSTANT_1_posterior_20211007.Rdata")

## 2nd simulation
#h.par <- seq(210,300,by=10)
#FULL.SIM.cte2 <- simSV_posterior(h.par=h.par,h.type="cte", npost=nPost, nsim=nROS, ti=100)
#FULL.SIM.cte2$dim <- FULL.SIM.cte2$dim[c("post","scen","h","sim","a","t")]
#save(FULL.SIM.cte2, file="reindeer data/simulations/FULL_SIM_CONSTANT_2_posterior_20211007.Rdata")

#________________________________________________________________________________________________________####
## Suppl. figure Constant vs. prop harvest ####
# load simulation data
memory.limit(size=500000)
#gc() # free unused memory
# proportional harvest simulations
load("reindeer data/simulations/FULL_SIM_posterior_20191213.Rdata") 
load("reindeer data/simulations/FULL_SIM_2_posterior_20211007.Rdata")
# constant harvest simulations
load("reindeer data/simulations/FULL_SIM_CONSTANT_1_posterior_20211007.Rdata")
load("reindeer data/simulations/FULL_SIM_CONSTANT_2_posterior_20211007.Rdata")

# parameters
h.par.p1 <- FULL.SIM$h.par
h.par.p2 <- FULL.SIM.prop2$h.par
h.par.p <- c(h.par.p1,h.par.p2)
h.par.c1 <- FULL.SIM.cte$h.par
h.par.c2 <- FULL.SIM.cte2$h.par
h.par.c <- c(h.par.c1,h.par.c2)

# var in r
Rv.p1 <- apply(apply(FULL.SIM$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)
Rv.p2 <- apply(apply(FULL.SIM.prop2$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)
Rv.p <- cbind(Rv.p1,Rv.p2)
Rv.c1 <- apply(apply(FULL.SIM.cte$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)
Rv.c2 <- apply(apply(FULL.SIM.cte2$r.sim, c(2,3,1,4), var),c(1,2),mean,na.rm=T)
Rv.c <- cbind(Rv.c1,Rv.c2)
# mean in N
Nm.p1 <- apply(FULL.SIM$Ntot.sim, c(2,3),mean)
Nm.p2 <- apply(FULL.SIM.prop2$Ntot.sim, c(2,3),mean)
Nm.p <- cbind(Nm.p1,Nm.p2)
Nm.c1 <- apply(FULL.SIM.cte$Ntot.sim, c(2,3),mean)
Nm.c2 <- apply(FULL.SIM.cte2$Ntot.sim, c(2,3),mean)
Nm.c <- cbind(Nm.c1,Nm.c2)
# var in log N
Xv.p1 <- apply(ifelse(FULL.SIM$Ntot.sim<350,NA,log(FULL.SIM$Ntot.sim)), c(2,3),var,na.rm=T)
Xv.p2 <- apply(ifelse(FULL.SIM.prop2$Ntot.sim<350,NA,log(FULL.SIM.prop2$Ntot.sim)), c(2,3),var,na.rm=T)
Xv.p <- cbind(Xv.p1,Xv.p2)
Xv.c1 <- apply(ifelse(FULL.SIM.cte$Ntot.sim<350,NA,log(FULL.SIM.cte$Ntot.sim)), c(2,3),var,na.rm=T)
Xv.c2 <- apply(ifelse(FULL.SIM.cte2$Ntot.sim<350,NA,log(FULL.SIM.cte2$Ntot.sim)), c(2,3),var,na.rm=T)
Xv.c <- cbind(Xv.c1,Xv.c2)

# quasi.ext
quasi.N <- round(AgeStr.BP[21,15]/5) # 20% of initital pop size in 2014 as threshold for quasi-extinction
quasi.N <- 350
# same as for generic model: K/5
#quasi.N <- 350
#quasi.N <- 100 # same as Hansen et al. 2019 Nat.Com. 
Q.p1 <- apply(apply(FULL.SIM$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),
              c(1,2),function(x)sum(x)/prod(FULL.SIM$dim[c("post","sim")]))
Q.p2 <- apply(apply(FULL.SIM.prop2$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),
              c(1,2),function(x)sum(x)/prod(FULL.SIM.prop2$dim[c("post","sim")]))
Q.p <- cbind(Q.p1,Q.p2)
Q.c1 <- apply(apply(FULL.SIM.cte$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),
              c(1,2),function(x)sum(x)/prod(FULL.SIM.cte$dim[c("post","sim")]))
Q.c2 <- apply(apply(FULL.SIM.cte2$Ntot.sim, c(2,3,1,4), function(x) any(x< quasi.N)),
              c(1,2),function(x)sum(x)/prod(FULL.SIM.cte2$dim[c("post","sim")]))
Q.c <- cbind(Q.c1,Q.c2)


# plot
#png(paste0("C:/R/04_Harvest-climate/FINAL script/Figures/Figure S4 reindeer propVScte.png"), 
#    res=300, units="in",w=6,h=10)

par(mfrow=c(4,2), mar=c(4.5,5,3,1), oma=c(0,0,.5,0),cex.lab=1.5,cex.axis=1.2,las=0,bty="l",tck=-0.02,
    mgp=c(3,0.5,0),las=1, yaxs="r")
lwd = 1; lty= c(3,2,1); col=1
#par(mfrow=c(1,2),mar=c(4.5,5,1.5,1),mgp=c(2.75,0.5,0),tck=-0.02,las=1,bty="l",cex.lab=1.2)
# var r
matplot(h.par.p, t(Rv.p), type="l",col=col, xlab="Harvest proportion",ylab=expression("Var("~italic(r)~")"),
        lwd=lwd,lty=lty)
legend("topright",bty="n", legend=c("Low","Medium","High"), title="Frequency of ROS events",col=col,
       lwd=lwd,lty=lty,cex=1.3)
mtext("(a)",3,line=0.5,adj=-.15,cex=1,font=2)
mtext("Proportional harvesting",line=2,font=2,cex=1)

matplot(h.par.c, t(Rv.c), type="l",col=col, xlab="Harvest constant",ylab=expression("Var("~italic(r)~")"),
        lwd=lwd,lty=lty)
mtext("(b)",3,line=0.5,adj=-.15,cex=1,font=2)
mtext("Constant harvesting",line=2,font=2,cex=1)

# N 
matplot(h.par.p, t(Nm.p), type="l", ylim=c(0,2000), lty=lty,xlab="Harvest proportion",
        ylab=expression("Mean("~italic(N)~")"),lwd=lwd,col=col)

mtext("(c)",3,line=0.5,adj=-.15,cex=1,font=2)

matplot(h.par.c, t(Nm.c), type="l", ylim=c(0,2000), lty=lty,xlab="Harvest constant",
        ylab=expression("Mean("~italic(N)~")"),lwd=lwd,col=col)
mtext("(d)",3,line=0.5,adj=-.15,cex=1,font=2)
# var N
matplot(h.par.p, t(Xv.p), type="l",ylim=c(range(Xv.c,Xv.p)), lty=lty,xlab="Harvest proportion",
        ylab=expression("Var( ln(" ~italic(N)~") )"),lwd=lwd,col=col)
mtext("(e)",3,line=0.5,adj=-.15,cex=1,font=2)

matplot(h.par.c, t(Xv.c), type="l",ylim=c(range(Xv.c,Xv.p)), lty=lty,xlab="Harvest constant",
        ylab=expression("Var( ln(" ~italic(N)~") )"),lwd=lwd,col=col)
mtext("(f)",3,line=0.5,adj=-.15,cex=1,font=2)


# quasi ext
matplot(h.par.p, t(Q.p), type="l",col=col, lty=lty,ylab="P(quasi-extinction)",lwd=lwd,xlab="Harvest proportion")
mtext("(g)",3,line=0.5,adj=-.15,cex=1,font=2)

matplot(h.par.c, t(Q.c), type="l",col=col, lty=lty,ylab="P(quasi-extinction)",lwd=lwd,xlab="Harvest constant")
mtext("(h)",3,line=0.5,adj=-.15,cex=1,font=2)

#dev.off()

