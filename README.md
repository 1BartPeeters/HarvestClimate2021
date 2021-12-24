# HarvestClimate2021
This repository is Archived and read-only. It contains the code and data for our paper "Harvesting can stabilize population fluctuations and buffer the impacts of extreme climatic events" [Peeters et al. 2022 Ecology Letters](https://www.authorea.com/users/395263/articles/545846-harvesting-can-stabilize-population-fluctuations-and-buffer-the-impacts-of-extreme-climatic-events?commit=5bdcca62ca0e5021fc1fc5c16c430f9b2d9129b4). The r-scripts provide the code to reproduce the simulations and model selection of population growth rate models decribed in the methods of our paper.

**theoretical model.R** = functions for logistic population growth rates and population trajectories with Ricker and Beverton-Holt types of density regulation, and additive (i.e. density-independent) vs. multiplicative (i.e. density-dependent) environmental noise. It also contains the code to reproduce the beanplots of variance in ln(N) and simulate quasi-extinction risks for varying harvest proportions and strenghts of environmnental drivers.

**TMB models.R** = code to build models using Template Model Builder (TMB) to test Ricker vs. Beverton-Holt logistic growth rates with additive vs. multiplicative noise, and climate covariates, on empirical data for six ungulates (provided in the folder "data_TMB"). See Table S1 in our paper for references to the published raw data. 

**reindeer simulations.R** = code to simulate effects of proportional and constant harvesting on Svalbard reindeer population trajectories for varying frequencies of rain-on-snow (ROS). The input parameters are estimated effects of ROS and N on survival and fecundity of six age classes of female Svalbard reindeer, derived from an Integrated Population Model (see Lee et al. 2015 Oikos, Hansen et al. 2019 Nature Communications). 
The folder "reindeer data" contains:
- parameter estimates for 9090 posterior models of age-specific survival (*coeffic.Rdata*) and fecundity (*coefficF.Rdata*), age-specific effects of interactions between ROS and density. 
- *SigmaInt.Rdata*: covariance matrices between age-specific survival and fecundity for 9090 posterior models
- *ROS.proj.Rdata* = array with simulations of ROS trajectories: 6000 simulations of 100 time steps for each of 5 scenarios (very low to very high frequencies of ROS; see Hansen et al. 2019 Nature Communications).
- *AgeStruc.Rdata* = age structures of 0 to 13+ year old female reindeer during 1994-2014. These were estimated based on cohort analyses as the IPM estimated annual population sizes for six age classes.

DOI: 10.5281/zenodo.5803068
