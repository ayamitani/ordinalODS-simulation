#----------------------------------------------------------------------------------------------------------------------------
#   Simulation program for ordinal ODS
#      
#   This program runs the simulation study presented in
#   "Applying survey-weighted ordinal regression models for unbiased inference in outcome-dependent samples with ordinal outcomes"
#   Authors: Aya Mitani, Osvaldo Espin-Garcia, Daniel Fernandez, Victoria Landsman
#   
#   Code written by: Aya Mitani
#
#----------------------------------------------------------------------------------------------------------------------------

### Load packages

library(MASS)
library(tidyverse)
library(VGAM)
library(nnet)
library(here)
library(xtable)
library(parallel)
library(svyVGAM)

### First read in the arguments listed at the command line
args <- (commandArgs(TRUE))
print(args)

### args is now a list of character vectors
### First check to see if arguments are passed.
### Then cycle through each element of the list and evaluate the expressions.

if(length(args) == 0){
  print("No arguments supplied.")
  
  ## supply default values
  baseseed <- 934
  numsim <- 2 # number of simulations to run in each core
  N <- 10000 # sample size in population
  n <- 80 # sample n from each outcome k = 1,...,K
  K <- 5 # number of outcome categories
  # proportion of patients in each category in the full cohort
  baseprob1 <- 0.8
  baseprob2 <- 0.05
  baseprob3 <- 0.05
  baseprob4 <- 0.05
  baseprob5 <- 0.05
  beta1 <- 1
  beta2 <- -0.5
}else{
  for (i in (1:length(args))) eval(parse(text=args[[i]]))
}



#------------------------------------------------------------------------
# Some functions used in simulation
#------------------------------------------------------------------------

outclm <- function(datain, name){
  clmmodel <- vglm(ordered(y) ~ x1 + x2, data = datain, family = cumulative(parallel = T, reverse = F))
  clm_est <- coef(clmmodel)
  names(clm_est) <- paste0(name, "_clm_est_", clmparams)
  clm_var <- diag(vcovvlm(clmmodel))
  names(clm_var) <- paste0(name, "_clm_var_", clmparams)
  clm_se <- sqrt(clm_var)
  names(clm_se) <- paste0(name, "_clm_se_", clmparams)
  outclm <- list(clm_est, clm_var, clm_se)
}

outacat <- function(datain, name){
  acatmodel <- vglm(ordered(y) ~ x1 + x2, data = datain, family = acat(parallel = TRUE, reverse = FALSE))
  acat_est <- coef(acatmodel)
  names(acat_est) <- paste0(name, "_acat_est_", acatparams)
  acat_var <- diag(vcovvlm(acatmodel))
  names(acat_var) <- paste0(name, "_acat_var_", acatparams)
  acat_se <- sqrt(acat_var)
  names(acat_se) <- paste0(name, "_acat_se_", acatparams)
  outacat <- list(acat_est, acat_var, acat_se)
}

outcratio <- function(datain, name){
  cratiomodel <- vglm(ordered(y) ~ x1 + x2, data = datain, family = cratio(parallel = TRUE, reverse = FALSE))
  cratio_est <- coef(cratiomodel)
  names(cratio_est) <- paste0(name, "_cratio_est_", cratioparams)
  cratio_var <- diag(vcovvlm(cratiomodel))
  names(cratio_var) <- paste0(name, "_cratio_var_", cratioparams)
  cratio_se <- sqrt(cratio_var)
  names(cratio_se) <- paste0(name, "_cratio_se_", cratioparams)
  outcratio <- list(cratio_est, cratio_var, cratio_se)
}

outmlm <- function(datain, name){
  mlm <- vglm(y ~ x1 + x2, data = datain, family = multinomial())
  mlm_est <- coef(mlm)
  names(mlm_est) <- paste0(name, "_mlm_est_", mlmparams)
  mlm_var <- diag(vcovvlm(mlm))
  names(mlm_var) <- paste0(name, "_mlm_var_", mlmparams)
  mlm_se <- sqrt(mlm_var)
  names(mlm_se) <- paste0(name, "_mlm_se_", mlmparams)
  outmlm <- list(mlm_est, mlm_var, mlm_se)
}


outwclm <- function(datain, name, svydes){
  clmmodel <- svy_vglm(ordered(y) ~ x1 + x2, design = svydes, family = cumulative(parallel = T, reverse = F))
  clm_est <- coef(clmmodel)
  names(clm_est) <- paste0(name, "_wclm_est_", clmparams)
  clm_var <- diag(vcov(clmmodel))
  names(clm_var) <- paste0(name, "_wclm_var_", clmparams)
  clm_se <- sqrt(clm_var)
  names(clm_se) <- paste0(name, "_wclm_se_", clmparams)
  outclm <- list(clm_est, clm_var, clm_se)
}

outwacat <- function(datain, name, svydes){
  acatmodel <- svy_vglm(ordered(y) ~ x1 + x2, design = svydes, family = acat(parallel = TRUE, reverse = FALSE))
  acat_est <- coef(acatmodel)
  names(acat_est) <- paste0(name, "_wacat_est_", acatparams)
  acat_var <- diag(vcov(acatmodel))
  names(acat_var) <- paste0(name, "_wacat_var_", acatparams)
  acat_se <- sqrt(acat_var)
  names(acat_se) <- paste0(name, "_wacat_se_", acatparams)
  outacat <- list(acat_est, acat_var, acat_se)
}

outwcratio <- function(datain, name, svydes){
  cratiomodel <- svy_vglm(ordered(y) ~ x1 + x2, design = svydes, family = cratio(parallel = TRUE, reverse = FALSE))
  cratio_est <- coef(cratiomodel)
  names(cratio_est) <- paste0(name, "_wcratio_est_", cratioparams)
  cratio_var <- diag(vcov(cratiomodel))
  names(cratio_var) <- paste0(name, "_wcratio_var_", cratioparams)
  cratio_se <- sqrt(cratio_var)
  names(cratio_se) <- paste0(name, "_wcratio_se_", cratioparams)
  outcratio <- list(cratio_est, cratio_var, cratio_se)
}

outwmlm <- function(datain, name, svydes){
  mlm <- svy_vglm(y ~ x1 + x2, design = svydes, family = multinomial())
  mlm_est <- coef(mlm)
  names(mlm_est) <- paste0(name, "_wmlm_est_", mlmparams)
  mlm_var <- diag(vcov(mlm))
  names(mlm_var) <- paste0(name, "_wmlm_var_", mlmparams)
  mlm_se <- sqrt(mlm_var)
  names(mlm_se) <- paste0(name, "_wmlm_se_", mlmparams)
  outmlm <- list(mlm_est, mlm_var, mlm_se)
}


relbiasf <- function(popbetas, sampbetas, name, parms){
  relbias <- (sampbetas - popbetas)/popbetas * 100
  names(relbias) <- paste0(name, "_relbias_", parms)
  relbias
}

absrelbiasf <- function(popbetas, sampbetas, name, parms){
  absrelbias <- abs(sampbetas - popbetas)/popbetas * 100
  names(absrelbias) <- paste0(name, "_absrelbias_", parms)
  absrelbias
}

sqerrf <- function(popbetas, sampbetas, name, parms){
  sqerr <- (sampbetas - popbetas) ** 2
  names(sqerr) <- paste0(name, "_sqerr_", parms)
  sqerr
}

###--------- end of functions --------###


k <- 1:K
K <- length(k)
K1 <- K - 1
K2 <- K - 2
neach <- n
prevs <- baseprobs <- c(baseprob1, baseprob2, baseprob3, baseprob4, baseprob5)
covs <- covs <- c("x1", "x2")
p <- length(c(beta1, beta2))


### parameters


clmparams <- acatparams <- cratioparams <- c(rep(paste0("a", 1:K1)), "beta1", "beta2")
mlmparams <- c(rep(paste0("a", 1:K1)), rep(paste0("beta1", 1:K1)), rep(paste0("beta2", 1:K1)))
smparams <- c(rep(paste0("phi", 1:K2)), rep(paste0("a", 1:K1)), "beta1", "beta2")

trueprob <- baseprobs
names(trueprob) <- paste0("trueprob", 1:K)
truebeta <- c(beta1, beta2)
names(truebeta) <- paste0("truebeta", 1:p)

ncores <- detectCores()
simtotal <- ncores * numsim



#######################################     Begin simulation     #########################################



ods_sim <- function(x){
  
  from <- seq(1, simtotal, by = numsim)[x]
  to <- seq(numsim, simtotal, by = numsim)[x]
  
  resultmat <- matrix(NA, numsim, 434)
  
  for (s in 1:numsim){
    
    sim <- (from + s - 1)
    
    simseed <- baseseed * sim + sim ^ 2
    set.seed(simseed)
    
    ###-------------------------------------
    ### generate full cohort
    ###-------------------------------------
    
    x1 <- rnorm(N, mean = 0, sd = 1)
    x2 <- rbinom(N, size = 1, prob = 0.5)
    z2 <- rnorm(N, mean = 28.5, sd = 4)
    z3 <- runif(N, min = 45, max = 75)
    
    xb <- x1 * beta1 + x2 * beta2 
    truebeta <- c(beta1, beta2)
    names(truebeta) <- paste0("truebeta", 1:p)
    

    alpha1 <- qnorm(prevs[1]) 
    alpha2 <- qnorm(prevs[1] + prevs[2]) 
    alpha3 <- qnorm(prevs[1] + prevs[2] + prevs[3]) 
    alpha4 <- qnorm(prevs[1] + prevs[2] + prevs[3] + prevs[4])
    
    
    trueavec <- c(-5000, alpha1, alpha2, alpha3, alpha4, 5000)
    names(trueavec) <- c("truea0", "truea1", "truea2", "truea3", "truea4", "truea5")
    
    yvectrue <- NULL
    problist <- list()
    
    for (k in 1:K){ 
      problist[[k]] <- pnorm(trueavec[k+1] - xb) - pnorm(trueavec[k] - xb)
    }
    
    probmatrix <- do.call("rbind", problist)
    
    ylist <- vector("list", N)
    ydummylist <- vector("list", N)
    
    for(i in 1:N){
      ydummy <- rmultinom(n = 1, size = 1, probmatrix[,i])
      for (q in 1:K)  
      { if (ydummy[q] == 1) { y = q } }
      ydummylist[[i]] <- ydummy
      ylist[[i]] <- y
    }  
    y <- as.factor(unlist(ylist))
    
    
    propy <- prop.table(table(y))
    names(propy) <- c("prop_y1", "prop_y2", "prop_y3", "prop_y4", "prop_y5")
    
    id <- c(1:N)
    popdata <- as.data.frame(cbind(id,y,x1,x2)) %>% 
      mutate(case = as.numeric(y > 1)) 
    

    tryCatch(clm_pop <- outclm(datain = popdata, name = "pop"))
    clm_popest <- clm_pop[[1]]
    clm_popvar <- clm_pop[[2]]
    clm_popse <- clm_pop[[3]]
    
    tryCatch(acat_pop <- outacat(datain = popdata, name = "pop"))
    acat_popest <- acat_pop[[1]]
    acat_popvar <- acat_pop[[2]]
    acat_popse <- acat_pop[[3]]
    
    tryCatch(cratio_pop <- outcratio(datain = popdata, name = "pop"))
    cratio_popest <- cratio_pop[[1]]
    cratio_popvar <- cratio_pop[[2]]
    cratio_popse <- cratio_pop[[3]]

    y_pop <- popdata$y
    x1_pop <- popdata$x1
    x2_pop <- popdata$x2
    tryCatch(sm <- rrvglm(-as.numeric(y_pop) ~ x1_pop + x2_pop, multinomial))
    sumsm <- summary(sm)
    sm_popest0 <- sumsm@coef3[,1]
    names(sm_popest0) <- paste0("pop", "_sm_est_", smparams)
    sm_popse0 <- sumsm@coef3[,2]
    names(sm_popse0) <- paste0("pop", "_sm_se_", smparams)
    sm_popvar0 <- sm_popse0 ^ 2
    names(sm_popvar0) <- paste0("pop", "_sm_var_", smparams)
    vcovsm <- vcov(sm)[c(1:3,8,9),c(1:3,8,9)]
    logodds1 <- c(1,sm_popest0[1:3]) * sm_popest0[8]
    var10 <- sm_popvar0[8]
    var11 <- sm_popest0[8]^2 * vcovsm[1,1] + sm_popest0[1]^2 * vcovsm[4,4] + 2 * sm_popest0[1] * sm_popest0[8] * vcovsm[1,4]
    var12 <- sm_popest0[8]^2 * vcovsm[2,2] + sm_popest0[2]^2 * vcovsm[4,4] + 2 * sm_popest0[2] * sm_popest0[8] * vcovsm[2,4]
    var13 <- sm_popest0[8]^2 * vcovsm[3,3] + sm_popest0[3]^2 * vcovsm[4,4] + 2 * sm_popest0[3] * sm_popest0[8] * vcovsm[3,4]
    logodds2 <- c(1,sm_popest0[1:3]) * sm_popest0[9]
    var20 <- sm_popvar0[9]
    var21 <- sm_popest0[9]^2 * vcovsm[1,1] + sm_popest0[1]^2 * vcovsm[5,5] + 2 * sm_popest0[1] * sm_popest0[9] * vcovsm[1,5]
    var22 <- sm_popest0[9]^2 * vcovsm[2,2] + sm_popest0[2]^2 * vcovsm[5,5] + 2 * sm_popest0[2] * sm_popest0[9] * vcovsm[2,5]
    var23 <- sm_popest0[9]^2 * vcovsm[3,3] + sm_popest0[3]^2 * vcovsm[5,5] + 2 * sm_popest0[3] * sm_popest0[9] * vcovsm[3,5]
    varlogodds1 <- c(var10,var11,var12,var13)
    varlogodds2 <- c(var20,var21,var22,var23)
    sm_popest2 <- c(logodds1, logodds2)
    sm_popvar2 <- c(varlogodds1, varlogodds2)
    sm_popse2 <- sqrt(sm_popvar0)
    sm_popor <- exp(c(logodds1, logodds2))
    names(sm_popest2) <- c(paste0("pop", "_sm_est_phi", c(0:(k-2)), "beta1"), paste0("pop", "_sm_est_phi", c(0:(k-2)), "beta2"))
    names(sm_popvar2) <- c(paste0("pop", "_sm_var_phi", c(0:(k-2)), "beta1"), paste0("pop", "_sm_var_phi", c(0:(k-2)), "beta2"))
    names(sm_popse2) <- c(paste0("pop", "_sm_se_phi", c(0:(k-2)), "beta1"), paste0("pop", "_sm_se_phi", c(0:(k-2)), "beta2"))
    sm_popest <- sm_popest0
    sm_popvar <- sm_popvar0
    sm_popse <- sm_popse0
    names(sm_popest) <- paste0("pop", "_sm_est_", smparams)
    names(sm_popvar) <- paste0("pop", "_sm_var_", smparams)
    names(sm_popse) <- paste0("pop", "_sm_se_", smparams)
    
    
    ###-------------------------------------
    ### simple random sampling
    ###-------------------------------------
    
    srs_sample0 <- popdata %>%
      dplyr::sample_n(neach * K) %>%
      mutate(srs = 1) %>%
      dplyr::select(id, srs)
    srs_sample <- left_join(popdata, srs_sample0, by = "id") %>%
      mutate(srs = replace_na(srs, 0)) %>%
      group_by(y) %>%
      mutate(sampwt = mean(srs),
             invsampwt = 1/sampwt) %>%
      ungroup() %>%
      filter(srs == 1) %>%
      as.data.frame()
    
    
    tryCatch(clm_srs <- outclm(datain = srs_sample, name = "srs"))
    clm_srsest <- clm_srs[[1]]
    clm_srsvar <- clm_srs[[2]]
    clm_srsse <- clm_srs[[3]]
    
    tryCatch(acat_srs <- outacat(datain = srs_sample, name = "srs"))
    acat_srsest <- acat_srs[[1]]
    acat_srsvar <- acat_srs[[2]]
    acat_srsse <- acat_srs[[3]]
    
    tryCatch(cratio_srs <- outcratio(datain = srs_sample, name = "srs"))
    cratio_srsest <- cratio_srs[[1]]
    cratio_srsvar <- cratio_srs[[2]]
    cratio_srsse <- cratio_srs[[3]]

    y_srs <- srs_sample$y
    x1_srs <- srs_sample$x1
    x2_srs <- srs_sample$x2
    tryCatch(sm <- rrvglm(-as.numeric(y_srs) ~ x1_srs + x2_srs, multinomial))
    sumsm <- summary(sm)
    sm_srsest0 <- sumsm@coef3[,1]
    names(sm_srsest0) <- paste0("srs", "_sm_est_", smparams)
    sm_srsse0 <- sumsm@coef3[,2]
    names(sm_srsse0) <- paste0("srs", "_sm_se_", smparams)
    sm_srsvar0 <- sm_srsse0 ^ 2
    names(sm_srsvar0) <- paste0("srs", "_sm_var_", smparams)
    vcovsm <- vcov(sm)[c(1:3,8,9),c(1:3,8,9)]
    logodds1 <- c(1,sm_srsest0[1:3]) * sm_srsest0[8]
    var10 <- sm_srsvar0[8]
    var11 <- sm_srsest0[8]^2 * vcovsm[1,1] + sm_srsest0[1]^2 * vcovsm[4,4] + 2 * sm_srsest0[1] * sm_srsest0[8] * vcovsm[1,4]
    var12 <- sm_srsest0[8]^2 * vcovsm[2,2] + sm_srsest0[2]^2 * vcovsm[4,4] + 2 * sm_srsest0[2] * sm_srsest0[8] * vcovsm[2,4]
    var13 <- sm_srsest0[8]^2 * vcovsm[3,3] + sm_srsest0[3]^2 * vcovsm[4,4] + 2 * sm_srsest0[3] * sm_srsest0[8] * vcovsm[3,4]
    logodds2 <- c(1,sm_srsest0[1:3]) * sm_srsest0[9]
    var20 <- sm_srsvar0[9]
    var21 <- sm_srsest0[9]^2 * vcovsm[1,1] + sm_srsest0[1]^2 * vcovsm[5,5] + 2 * sm_srsest0[1] * sm_srsest0[9] * vcovsm[1,5]
    var22 <- sm_srsest0[9]^2 * vcovsm[2,2] + sm_srsest0[2]^2 * vcovsm[5,5] + 2 * sm_srsest0[2] * sm_srsest0[9] * vcovsm[2,5]
    var23 <- sm_srsest0[9]^2 * vcovsm[3,3] + sm_srsest0[3]^2 * vcovsm[5,5] + 2 * sm_srsest0[3] * sm_srsest0[9] * vcovsm[3,5]
    varlogodds1 <- c(var10,var11,var12,var13)
    varlogodds2 <- c(var20,var21,var22,var23)
    sm_srsest2 <- c(logodds1, logodds2)
    sm_srsvar2 <- c(varlogodds1, varlogodds2)
    sm_srsse2 <- sqrt(sm_srsvar0)
    sm_srsor <- exp(c(logodds1, logodds2))
    names(sm_srsest2) <- c(paste0("srs", "_sm_est_phi", c(0:(k-2)), "beta1"), paste0("srs", "_sm_est_phi", c(0:(k-2)), "beta2"))
    names(sm_srsvar2) <- c(paste0("srs", "_sm_var_phi", c(0:(k-2)), "beta1"), paste0("srs", "_sm_var_phi", c(0:(k-2)), "beta2"))
    names(sm_srsse2) <- c(paste0("srs", "_sm_se_phi", c(0:(k-2)), "beta1"), paste0("srs", "_sm_se_phi", c(0:(k-2)), "beta2"))
    sm_srsest <- sm_srsest0
    sm_srsvar <- sm_srsvar0
    sm_srsse <- sm_srsse0
    names(sm_srsest) <- paste0("srs", "_sm_est_", smparams)
    names(sm_srsvar) <- paste0("srs", "_sm_var_", smparams)
    names(sm_srsse) <- paste0("srs", "_sm_se_", smparams)
    
    
    ###-------------------------------------
    ### Outcome-dependent sampling
    ###-------------------------------------
    
    ### unweighted models
    
    ods_sample0 <- popdata %>%
      dplyr::group_by(y) %>%
      dplyr::mutate(num_rows=n()) %>%
      dplyr::sample_n(neach) %>%
      ungroup %>%    
      dplyr::mutate(ods = 1) %>%
      dplyr::select(id, ods)
    ods_sample <- left_join(popdata, ods_sample0, by = "id") %>%
      dplyr::mutate(ods = replace_na(ods, 0)) %>%
      dplyr::group_by(y) %>%
      dplyr::mutate(sampwt = mean(ods),
             invsampwt = 1/sampwt) %>%
      ungroup() %>%
      dplyr::filter(ods == 1) %>%
      as.data.frame()
    
    tryCatch(clm_ods <- outclm(datain = ods_sample, name = "ods"))
    clm_odsest <- clm_ods[[1]]
    clm_odsvar <- clm_ods[[2]]
    clm_odsse <- clm_ods[[3]]
    
    tryCatch(acat_ods <- outacat(datain = ods_sample, name = "ods"))
    acat_odsest <- acat_ods[[1]]
    acat_odsvar <- acat_ods[[2]]
    acat_odsse <- acat_ods[[3]]
    
    tryCatch(cratio_ods <- outcratio(datain = ods_sample, name = "ods"))
    cratio_odsest <- cratio_ods[[1]]
    cratio_odsvar <- cratio_ods[[2]]
    cratio_odsse <- cratio_ods[[3]]

    
    y_ods <- ods_sample$y
    x1_ods <- ods_sample$x1
    x2_ods <- ods_sample$x2
    w_ods <- ods_sample$invsampwt
    tryCatch(sm <- rrvglm(-as.numeric(y_ods) ~ x1_ods + x2_ods, multinomial))
    sumsm <- summary(sm)
    sm_odsest0 <- sumsm@coef3[,1]
    sm_odsse0 <- sumsm@coef3[,2]
    sm_odsvar0 <- sm_odsse0 ^ 2
    vcovsm <- vcov(sm)[c(1:3,8,9),c(1:3,8,9)]
    logodds1 <- c(1,sm_odsest0[1:3]) * sm_odsest0[8]
    var10 <- sm_odsvar0[8]
    var11 <- sm_odsest0[8]^2 * vcovsm[1,1] + sm_odsest0[1]^2 * vcovsm[4,4] + 2 * sm_odsest0[1] * sm_odsest0[8] * vcovsm[1,4]
    var12 <- sm_odsest0[8]^2 * vcovsm[2,2] + sm_odsest0[2]^2 * vcovsm[4,4] + 2 * sm_odsest0[2] * sm_odsest0[8] * vcovsm[2,4]
    var13 <- sm_odsest0[8]^2 * vcovsm[3,3] + sm_odsest0[3]^2 * vcovsm[4,4] + 2 * sm_odsest0[3] * sm_odsest0[8] * vcovsm[3,4]
    logodds2 <- c(1,sm_odsest0[1:3]) * sm_odsest0[9]
    var20 <- sm_odsvar0[9]
    var21 <- sm_odsest0[9]^2 * vcovsm[1,1] + sm_odsest0[1]^2 * vcovsm[5,5] + 2 * sm_odsest0[1] * sm_odsest0[9] * vcovsm[1,5]
    var22 <- sm_odsest0[9]^2 * vcovsm[2,2] + sm_odsest0[2]^2 * vcovsm[5,5] + 2 * sm_odsest0[2] * sm_odsest0[9] * vcovsm[2,5]
    var23 <- sm_odsest0[9]^2 * vcovsm[3,3] + sm_odsest0[3]^2 * vcovsm[5,5] + 2 * sm_odsest0[3] * sm_odsest0[9] * vcovsm[3,5]
    varlogodds1 <- c(var10,var11,var12,var13)
    varlogodds2 <- c(var20,var21,var22,var23)
    sm_odsest2 <- c(logodds1, logodds2)
    sm_odsvar2 <- c(varlogodds1, varlogodds2)
    sm_odsse2 <- sqrt(sm_odsvar2)
    sm_odsor <- exp(c(logodds1, logodds2))
    names(sm_odsest2) <- c(paste0("ods", "_sm_est_phi", c(0:(k-2)), "beta1"), paste0("ods", "_sm_est_phi", c(0:(k-2)), "beta2"))
    names(sm_odsvar2) <- c(paste0("ods", "_sm_var_phi", c(0:(k-2)), "beta1"), paste0("ods", "_sm_var_phi", c(0:(k-2)), "beta2"))
    names(sm_odsse2) <- c(paste0("ods", "_sm_se_phi", c(0:(k-2)), "beta1"), paste0("ods", "_sm_se_phi", c(0:(k-2)), "beta2"))
    sm_odsest <- sm_odsest0
    sm_odsvar <- sm_odsvar0
    sm_odsse <- sm_odsse0
    names(sm_odsest) <- paste0("ods", "_sm_est_", smparams)
    names(sm_odsvar) <- paste0("ods", "_sm_var_", smparams)
    names(sm_odsse) <- paste0("ods", "_sm_se_", smparams)
    
    
    
    ### weighted models
    
    odsdes <- svydesign(id = ~1,
                        strata = ~ y, 
                        weights = ~ invsampwt, 
                        data = ods_sample)
    
    tryCatch(wclm_ods <- outwclm(datain = ods_sample, name = "ods", svydes = odsdes))
    wclm_odsest <- wclm_ods[[1]]
    wclm_odsvar <- wclm_ods[[2]]
    wclm_odsse <- wclm_ods[[3]]
    
    tryCatch(wacat_ods <- outwacat(datain = ods_sample, name = "ods", svydes = odsdes))
    wacat_odsest <- wacat_ods[[1]]
    wacat_odsvar <- wacat_ods[[2]]
    wacat_odsse <- wacat_ods[[3]]
    
    tryCatch(wcratio_ods <- outwcratio(datain = ods_sample, name = "ods", svydes = odsdes))
    wcratio_odsest <- wcratio_ods[[1]]
    wcratio_odsvar <- wcratio_ods[[2]]
    wcratio_odsse <- wcratio_ods[[3]]

    
    tryCatch(wsm <- rrvglm(-as.numeric(y_ods) ~ x1_ods + x2_ods, multinomial, weights = w_ods))
    sumwsm <- summary(wsm)
    wsm_odsest0 <- sumwsm@coef3[,1]
    inv_inf <- sumwsm@cov.unscaled
    scores <- weights(wsm, deriv = TRUE, type = "working")$deriv
    pwts <- weights(odsdes, "sampling")
    cons <- do.call(cbind, constraints(wsm))
    mmat <- model.matrix(wsm)
    case_index <- as.numeric(gsub("^(.+):.*", "\\1", rownames(mmat)))    
    mmatsum <- t(t(rowsum(mmat, case_index, reorder = FALSE )) / colSums(cons))
    scorethetabeta <- (((scores / pwts) %*% cons) * mmatsum)
    betahat <- coef(wsm)[K:(K+(p-1))]
    xmat <- as.matrix(wsm@x[,-1])
    xbeta <- xmat %*% betahat
    scorephi <- residuals(wsm, type = "response")[,c(-1,-K)] * as.vector(xbeta)
    scores <- cbind(scorephi, scorethetabeta)
    inffuns <- (scores) %*% inv_inf
    wsmvcov <- vcov(svytotal(inffuns, odsdes))
    wsm_odsvar0 <- diag(wsmvcov)
    wsm_odsse0 <- sqrt(wsm_odsvar0)
    
    vcovsm <- wsmvcov[c(1:3,8,9),c(1:3,8,9)]
    logodds1 <- c(1,wsm_odsest0[1:3]) * wsm_odsest0[8]
    var10 <- wsm_odsvar0[8]
    var11 <- wsm_odsest0[8]^2 * vcovsm[1,1] + wsm_odsest0[1]^2 * vcovsm[4,4] + 2 * wsm_odsest0[1] * wsm_odsest0[8] * vcovsm[1,4]
    var12 <- wsm_odsest0[8]^2 * vcovsm[2,2] + wsm_odsest0[2]^2 * vcovsm[4,4] + 2 * wsm_odsest0[2] * wsm_odsest0[8] * vcovsm[2,4]
    var13 <- wsm_odsest0[8]^2 * vcovsm[3,3] + wsm_odsest0[3]^2 * vcovsm[4,4] + 2 * wsm_odsest0[3] * wsm_odsest0[8] * vcovsm[3,4]
    logodds2 <- c(1,wsm_odsest0[1:3]) * wsm_odsest0[9]
    var20 <- wsm_odsvar0[9]
    var21 <- wsm_odsest0[9]^2 * vcovsm[1,1] + wsm_odsest0[1]^2 * vcovsm[5,5] + 2 * wsm_odsest0[1] * wsm_odsest0[9] * vcovsm[1,5]
    var22 <- wsm_odsest0[9]^2 * vcovsm[2,2] + wsm_odsest0[2]^2 * vcovsm[5,5] + 2 * wsm_odsest0[2] * wsm_odsest0[9] * vcovsm[2,5]
    var23 <- wsm_odsest0[9]^2 * vcovsm[3,3] + wsm_odsest0[3]^2 * vcovsm[5,5] + 2 * wsm_odsest0[3] * wsm_odsest0[9] * vcovsm[3,5]
    varlogodds1 <- c(var10,var11,var12,var13)
    varlogodds2 <- c(var20,var21,var22,var23)
    wsm_odsest2 <- c(logodds1, logodds2)
    wsm_odsvar2 <- c(varlogodds1, varlogodds2)
    wsm_odsse2 <- sqrt(wsm_odsvar2)
    wsm_odsor <- exp(c(logodds1, logodds2))
    names(wsm_odsest2) <- c(paste0("ods", "_wsm_est_phi", c(0:(k-2)), "beta1"), paste0("ods", "_wsm_est_phi", c(0:(k-2)), "beta2"))
    names(wsm_odsvar2) <- c(paste0("ods", "_wsm_var_phi", c(0:(k-2)), "beta1"), paste0("ods", "_wsm_var_phi", c(0:(k-2)), "beta2"))
    names(wsm_odsse2) <- c(paste0("ods", "_wsm_se_phi", c(0:(k-2)), "beta1"), paste0("ods", "_wsm_se_phi", c(0:(k-2)), "beta2"))
    wsm_odsest <- wsm_odsest0
    wsm_odsvar <- wsm_odsvar0
    wsm_odsse <- wsm_odsse0
    names(wsm_odsest) <- paste0("ods", "_wsm_est_", smparams)
    names(wsm_odsvar) <- paste0("ods", "_wsm_var_", smparams)
    names(wsm_odsse) <- paste0("ods", "_wsm_se_", smparams)
    
    mean_sampwt <- mean(ods_sample$invsampwt)
    sd_sampwt <- sd(ods_sample$invsampwt)
    cv_sampwt <- sd_sampwt/mean_sampwt 
    
    names(mean_sampwt) <- "mean_sampwt"
    names(sd_sampwt) <- "sd_sampwt"
    names(cv_sampwt) <- "cv_sampwt"
    
    names(sim) <- "sim"
    names(simseed) <- "simseed"
    
    resultvec <- c(sim, simseed,
                   trueprob, truebeta,
                   clm_popest, clm_popvar, clm_popse,
                   acat_popest, acat_popvar, acat_popse,
                   cratio_popest, cratio_popvar, cratio_popse, 
                   sm_popest, sm_popvar, sm_popse, sm_popest2, sm_popvar2, sm_popse2,   
                   clm_srsest, clm_srsvar, clm_srsse, 
                   acat_srsest, acat_srsvar, acat_srsse, 
                   cratio_srsest, cratio_srsvar, cratio_srsse, 
                   sm_srsest, sm_srsvar, sm_srsse, sm_srsest2, sm_srsvar2, sm_srsse2,   
                   clm_odsest, clm_odsvar, clm_odsse, 
                   acat_odsest, acat_odsvar, acat_odsse, 
                   cratio_odsest, cratio_odsvar, cratio_odsse, 
                   sm_odsest, sm_odsvar, sm_odsse, sm_odsest2, sm_odsvar2, sm_odsse2, 
                   wclm_odsest, wclm_odsvar, wclm_odsse, 
                   wacat_odsest, wacat_odsvar, wacat_odsse, 
                   wcratio_odsest, wcratio_odsvar, wcratio_odsse, 
                   wsm_odsest, wsm_odsvar, wsm_odsse, wsm_odsest2, wsm_odsvar2, wsm_odsse2,
                   mean_sampwt, sd_sampwt, cv_sampwt
    )
    resultmat[s, ] <- resultvec
    colnames(resultmat) <- names(resultvec)
    print(s)
    
  }
  
  return(resultmat)

}
  

resultmatout <- mclapply(1:ncores, ods_sim, mc.cores = ncores, mc.silent = F)

simout <- as.data.frame(do.call(rbind, resultmatout))

outfile <- paste0("simout_probit_N", N,
                  "_neach", neach,
                  "_K", K,
                  paste0("_trueprob", 1:K, trueprob, collapse = ""),
                  paste0("_truebeta", 1:p, truebeta, collapse = ""),
                  ".txt")

write.table(simout, outfile, col.names = T, row.names = F)




