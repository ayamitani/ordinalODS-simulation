#----------------------------------------------------------------------------------------------------------------------------
#   Simulation program for ordinal ODS - cumulative logit model (CM)
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
  baseseed <- 1111
  numsim <- 6 # number of simulations to run
  N <- 10000 # sample size in population
  n <- 80 # sample n from each outcome k = 1,...,K
  K <- 5 # number of outcome categories
  scenarioL <- 1 # indicator for scenario L
  scenarioM <- 0 # indicator for scenario M
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
  cratiomodel <- vglm(ordered(y) ~ x1 + x2, data = datain, family = cratio(parallel = TRUE, reverse = T))
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
  cratiomodel <- svy_vglm(ordered(y) ~ x1 + x2, design = svydes, family = cratio(parallel = TRUE, reverse = T))
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

load("estparam.RData")

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
    dplyr::filter(scenario == scen & model == "cm") -> trueparam
  beta1 <- as.numeric(trueparam["beta1"])
  beta2 <- as.numeric(trueparam["beta2"])
  p <- length(c(beta1, beta2))
  alpha1 <- as.numeric(trueparam["alpha1"])
  alpha2 <- as.numeric(trueparam["alpha2"])
  alpha3 <- as.numeric(trueparam["alpha3"])
  alpha4 <- as.numeric(trueparam["alpha4"])
  alpha <- c(alpha1, alpha2, alpha3, alpha4)
  truealpha <- alpha
  truebeta <- c(beta1, beta2)
  names(truebeta) <- paste0("truebeta", 1:p)
  names(truealpha) <- paste0("truealpha", 1:(K-1))
  
  truecm <- c(truealpha, truebeta)
  
  ncores <- detectCores()
  simtotal <- ncores * numsim
  
  
  
  #######################################     Begin simulation     #########################################
  
  
  
  ods_sim <- function(x){
    
    from <- seq(1, simtotal, by = numsim)[x]
    to <- seq(numsim, simtotal, by = numsim)[x]
    
    resultmat <- matrix(NA, numsim, 106)
    
    for (s in 1:numsim){
      
      sim <- (from + s - 1)
      
      print(c(from, to, sim))
      
      simseed <- baseseed * 2 * sim * scenarionum ^ 2
      
      print(simseed)
      set.seed(simseed)
      
      
      
      
      ### population data
      x1 <- runif(N)
      x2 <- rbinom(N, size = 1, prob = 0.5)
      
      xb <- x1 * beta1 + x2 * beta2 
      truebeta <- c(beta1, beta2)
      names(truebeta) <- paste0("truebeta", 1:p)
      
      
      # generate log odds 
      #e1 <- exp(alpha[1] + phi[1] * (beta1 * x1 + beta2 * x2))  # e1 = 1
      e2 <- exp(truealpha[1] + (beta1 * x1 + beta2 * x2))
      e3 <- exp(truealpha[2] + (beta1 * x1 + beta2 * x2))
      e4 <- exp(truealpha[3] + (beta1 * x1 + beta2 * x2))
      e5 <- exp(truealpha[4] + (beta1 * x1 + beta2 * x2))
      
      # generate probabilities
      p1 <- e2 / (1 + e2)
      p2 <- e3 / (1 + e3) - e2 / (1 + e2) 
      p3 <- e4 / (1 + e4) - e3 / (1 + e3) 
      p4 <- e5 / (1 + e5) - e4 / (1 + e4) 
      p5 <- 1 - p4 - p3 - p2 - p1
      
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
      
      
      tryCatch(clm_pop <- outclm(datain = popdata, name = "pop"))
      clm_popest <- clm_pop[[1]]
      clm_popvar <- clm_pop[[2]]
      clm_popse <- clm_pop[[3]]
      
      
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
      
      
      tryCatch(clm_ods <- outclm(datain = ods_sample, name = "ods"))
      clm_odsest <- clm_ods[[1]]
      clm_odsvar <- clm_ods[[2]]
      clm_odsse <- clm_ods[[3]]
      
      
      ### with ipw
      odsdes <- svydesign(id = ~1,
                          strata = ~ y, 
                          weights = ~ invsampwt, 
                          data = ods_sample)
      
      
      tryCatch(wclm_ods <- outwclm(datain = ods_sample, name = "ods", svydes = odsdes))
      wclm_odsest <- wclm_ods[[1]]
      wclm_odsvar <- wclm_ods[[2]]
      wclm_odsse <- wclm_ods[[3]]
      
      
      mean_sampwt <- mean(ods_sample$invsampwt)
      sd_sampwt <- sd(ods_sample$invsampwt)
      cv_sampwt <- sd_sampwt/mean_sampwt 
      
      names(mean_sampwt) <- "cratio_mean_sampwt"
      names(sd_sampwt) <- "cratio_sd_sampwt"
      names(cv_sampwt) <- "clm_cv_sampwt"
      
      names(sim) <- "sim"
      names(simseed) <- "simseed"
      
      # bias
      clm_bias <- (clm_odsest - truecm)
      names(clm_bias) <- paste0("ods_clm", "_bias_", clmparams)
      wclm_bias <- (wclm_odsest - truecm)
      names(wclm_bias) <- paste0("ods_wclm", "_bias_", clmparams)
      clm_absbias <- abs(clm_bias)
      names(clm_absbias) <- paste0("ods_clm", "_absbias_", clmparams)
      wclm_absbias <- abs(wclm_bias)
      names(wclm_absbias) <- paste0("ods_wclm", "_absbias_", clmparams)
      clm_sqerr <- clm_bias ^ 2
      names(clm_sqerr) <- paste0("ods_clm", "_sqerr_", clmparams)
      wclm_sqrtt <- wclm_bias ^ 2
      names(wclm_sqrtt) <- paste0("ods_wclm", "_sqerr_", clmparams)
      
      resultvec <- c(sim, simseed,
                     propy, truealpha, truebeta,
                     clm_popest, clm_popvar, clm_popse, 
                     clm_odsest, clm_odsvar, clm_odsse, 
                     wclm_odsest, wclm_odsvar, wclm_odsse, 
                     clm_bias, clm_absbias, clm_sqerr, 
                     wclm_bias, wclm_absbias, wclm_sqrtt,
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
  
  outfile <- paste0("simout_cm_N", N,
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




