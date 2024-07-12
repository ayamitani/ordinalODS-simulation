#----------------------------------------------------------------------------------------------------------------------------
#   Simulation program for ordinal ODS - stereotype model (SM)
#      
#   This program runs the simulation study presented in
#   "Applying survey-weighted proportional odds models for improved inference in outcome-dependent samples with ordinal outcomes"
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
  baseseed <- 4445
  numsim <- 2 # number of simulations to run
  N <- 10000 # sample size in population
  n <- 80 # sample n from each outcome k = 1,...,K
  K <- 5 # number of outcome categories
  phi0 <- 0
  phi1 <- 0.25
  phi2 <- 0.5
  phi3 <- 0.75
  scenarioL <- 0 # indicator for scenario L
  scenarioM <- 1 # indicator for scenario M
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
covs <- covs <- c("x1", "x2")
p <- length(covs)


### parameters


clmparams <- acatparams <- cratioparams <- c(rep(paste0("a", 1:K1)), "beta1", "beta2")
mlmparams <- c(rep(paste0("a", 1:K1)), rep(paste0("beta1", 1:K1)), rep(paste0("beta2", 1:K1)))
smparams <- c(rep(paste0("phi", 1:K2)), rep(paste0("a", 1:K1)), "beta1", "beta2")

ncores <- detectCores()
simtotal <- ncores * numsim

#ncores <- 1

load("estparam.RData")
#load(here("Simulation", "estparam.RData"))

runsim <- function(scenarionum){
  
  scenario <- if(scenarioL == 1){
    "L"
  }else if(scenarioM == 1){
    "M"
  }else{
    "H"
  }
  scenarionumroman <- c("i", "ii", "iii", "iv", "v")
  scen <- paste0(scenario, scenarionumroman[scenarionum])
  estparam %>%
    dplyr::filter(scenario == scen & model == "sm") -> trueparam
  beta1 <- as.numeric(trueparam["beta1"])
  beta2 <- as.numeric(trueparam["beta2"])
  p <- length(c(beta1, beta2))
  alpha1 <- as.numeric(trueparam["alpha1"])
  alpha2 <- as.numeric(trueparam["alpha2"])
  alpha3 <- as.numeric(trueparam["alpha3"])
  alpha4 <- as.numeric(trueparam["alpha4"])
  alpha <- c(alpha1, alpha2, alpha3, alpha4, 0)
  
  truealpha <- alpha
  names(truealpha) <- paste0("truealpha", 1:K)
  truebeta <- c(beta1, beta2)
  names(truebeta) <- paste0("truebeta", 1:p)
  phi1 <- as.numeric(trueparam["phi1"])
  phi2 <- as.numeric(trueparam["phi2"])
  phi3 <- as.numeric(trueparam["phi3"])
  truephi <- phi <- c(0, phi1, phi2, phi3, 1)
  names(truephi) <- paste0("truephi", 0:(K-1))
  truesm <- c((truealpha)[-K], 
              rev(truephi[-1])*truebeta[1], rev(truephi[-1])*truebeta[2])
  names(truesm) <- c(paste0("truesmalpha", 1:K1), paste0("truesm_phi", c(0:(K-2)), "beta1"), paste0("truesm_phi", c(0:(K-2)), "beta2"))
  #
  
  
  #######################################     Begin simulation     #########################################
  
  
  
  ods_sim <- function(x){
    
    from <- seq(1, simtotal, by = numsim)[x]
    to <- seq(numsim, simtotal, by = numsim)[x]
    
    resultmat <- matrix(NA, numsim, 259)
    
    for (s in 1:numsim){
      
      sim <- (from + s - 1)
      
      simseed <- baseseed * 2 * sim * scenarionum ^ 2
      set.seed(simseed)
      
      print(simseed)
      
      
      ### population data
      x1 <- runif(N)
      x2 <- rbinom(N, size = 1, prob = 0.5)
      
      xb <- x1 * beta1 + x2 * beta2 
      truebeta <- c(beta1, beta2)
      names(truebeta) <- paste0("truebeta", 1:p)
      
      # generate log odds 
      e2 <- exp(alpha[4] + phi[2] * (beta1 * x1 + beta2 * x2))
      e3 <- exp(alpha[3] + phi[3] * (beta1 * x1 + beta2 * x2))
      e4 <- exp(alpha[2] + phi[4] * (beta1 * x1 + beta2 * x2))
      e5 <- exp(alpha[1] + phi[5] * (beta1 * x1 + beta2 * x2))

      # generate probabilities
      p2 <- e2 / (1 + e2 + e3 + e4 + e5)
      p3 <- e3 / (1 + e2 + e3 + e4 + e5)
      p4 <- e4 / (1 + e2 + e3 + e4 + e5)
      p5 <- e5 / (1 + e2 + e3 + e4 + e5)
      p1 <- 1 - p2 - p3 - p4 - p5
      
      # create empty list of outcome
      ylist <- vector("list", N)
      ydummylist <- vector("list", N)
      
      
      for(i in 1:N){
        ydummy <- rmultinom(n = 1, size = 1, prob = c(p1[i], p2[i], p3[i], p4[i], p5[i])) # generate outcome 0/1
        for (q in 1:5)  # convert 0/1 to y=k
        { if (ydummy[q] == 1) { y = q } }
        ydummylist[[i]] <- ydummy
        ylist[[i]] <- y
      }
      
      y <- as.factor(unlist(ylist))
      
      
      propy <- prop.table(table(y))
      names(propy) <- c("simprop1", "simprop2", "simprop3", "simprop4", "simprop5")
      propy
      
      id <- c(1:N)
      popdata <- as.data.frame(cbind(id,y,x1,x2)) %>% 
        mutate(case = as.numeric(y > 1)) 
      
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
      sm_popse2 <- sqrt(sm_popvar2)
      sm_popor <- exp(c(logodds1, logodds2))
      names(sm_popest2) <- c(paste0("pop", "_sm_est_phi", c(0:(K-2)), "beta1"), paste0("pop", "_sm_est_phi", c(0:(K-2)), "beta2"))
      names(sm_popvar2) <- c(paste0("pop", "_sm_var_phi", c(0:(K-2)), "beta1"), paste0("pop", "_sm_var_phi", c(0:(K-2)), "beta2"))
      names(sm_popse2) <- c(paste0("pop", "_sm_se_phi", c(0:(K-2)), "beta1"), paste0("pop", "_sm_se_phi", c(0:(K-2)), "beta2"))
      sm_popest <- sm_popest0
      sm_popvar <- sm_popvar0
      sm_popse <- sm_popse0
      names(sm_popest) <- paste0("pop", "_sm_est_", smparams)
      names(sm_popvar) <- paste0("pop", "_sm_var_", smparams)
      names(sm_popse) <- paste0("pop", "_sm_se_", smparams)
      
      
      
      ### Outcome-dependent sampling
      
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
      names(sm_odsest2) <- c(paste0("ods", "_sm_est_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_est_phi", c(0:(K-2)), "beta2"))
      names(sm_odsvar2) <- c(paste0("ods", "_sm_var_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_var_phi", c(0:(K-2)), "beta2"))
      names(sm_odsse2) <- c(paste0("ods", "_sm_se_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_se_phi", c(0:(K-2)), "beta2"))
      sm_odsest <- sm_odsest0
      sm_odsvar <- sm_odsvar0
      sm_odsse <- sm_odsse0
      names(sm_odsest) <- paste0("ods", "_sm_est_", smparams)
      names(sm_odsvar) <- paste0("ods", "_sm_var_", smparams)
      names(sm_odsse) <- paste0("ods", "_sm_se_", smparams)
      
      
      
      ### with ipw
      odsdes <- svydesign(id = ~1,
                          strata = ~ y, 
                          weights = ~ invsampwt, 
                          data = ods_sample)
      
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
      names(wsm_odsest2) <- c(paste0("ods", "_wsm_est_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_est_phi", c(0:(K-2)), "beta2"))
      names(wsm_odsvar2) <- c(paste0("ods", "_wsm_var_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_var_phi", c(0:(K-2)), "beta2"))
      names(wsm_odsse2) <- c(paste0("ods", "_wsm_se_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_se_phi", c(0:(K-2)), "beta2"))
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
      
      
      # sm bias
      sm_odsest3 <- c(sm_odsest0[4:7], sm_odsest2)
      wsm_odsest3 <- c(wsm_odsest0[4:7], wsm_odsest2)
      sm_bias <- (sm_odsest3 - truesm)
      names(sm_bias) <- c(paste0("ods", "_sm_bias_", smparams)[4:7], paste0("ods", "_sm_bias_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_bias_phi", c(0:(K-2)), "beta2"))
      wsm_bias <- (wsm_odsest3 - truesm)
      names(wsm_bias) <- c(paste0("ods", "_wsm_bias_", smparams)[4:7], paste0("ods", "_wsm_bias_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_bias_phi", c(0:(K-2)), "beta2"))
      sm_absbias <- abs(sm_bias)
      names(sm_absbias) <- c(paste0("ods_sm", "_absbias_", smparams)[4:7], paste0("ods", "_sm_absbias_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_absbias_phi", c(0:(K-2)), "beta2"))
      wsm_absbias <- abs(wsm_bias)
      names(wsm_absbias) <- c(paste0("ods_wsm", "_absbias_", smparams)[4:7], paste0("ods", "_wsm_absbias_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_absbias_phi", c(0:(K-2)), "beta2"))
      sm_sqerr <- sm_bias ^ 2
      names(sm_sqerr) <- c(paste0("ods", "_sm_sqerr_", smparams)[4:7], paste0("ods", "_sm_sqerr_phi", c(0:(K-2)), "beta1"), paste0("ods", "_sm_sqerr_phi", c(0:(K-2)), "beta2"))
      wsm_sqerr <- wsm_bias ^ 2
      names(wsm_sqerr) <- c(paste0("ods", "_wsm_sqerr_", smparams)[4:7], paste0("ods", "_wsm_sqerr_phi", c(0:(K-2)), "beta1"), paste0("ods", "_wsm_sqerr_phi", c(0:(K-2)), "beta2"))
      
      
      resultvec <- c(sim, simseed,
                     truealpha, truebeta, truephi, truesm, 
                     propy,
                     sm_popest, sm_popvar, sm_popse, sm_popest2, sm_popvar2, sm_popse2,   
                     sm_odsest, sm_odsvar, sm_odsse, sm_odsest2, sm_odsvar2, sm_odsse2, 
                     wsm_odsest, wsm_odsvar, wsm_odsse, wsm_odsest2, wsm_odsvar2, wsm_odsse2,
                     sm_bias, sm_absbias, sm_sqerr,
                     wsm_bias, wsm_absbias, wsm_sqerr,
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
  
  outfile <- paste0("simout_sm_N", N,
                    "_neach", neach,
                    "_K", K,
                    "_scenario", scen,
                    ".txt")
  
  write.table(simout, outfile, col.names = T, row.names = F)
  
}



runsim(1)
runsim(2)
runsim(3)
runsim(4)
runsim(5)
