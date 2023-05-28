# require(dplyr)
# require(Iso)
# fun_BOINET(c(0.1, 0.2, 0.4, 0.6), c(0.1, 0.7, 0.2, 0.1), trueOBD = 2,
#           pT = 0.35, qE = 0.2, csize = 3, cN = 9)

fun_BOINET <- function(plist, qlist, trueOBD, pT, qE, csize, cN,  
                       phi = 0.3, delta = 0.6, 
                       onlyDM = TRUE, alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,
                       current = 1, doselimit = Inf, BOINEThyper = NULL, eta_et = NULL, xi_et = NULL)
{
  # do not consider dose elimination 
  BOINET_one <- function(doseDT, current, csize, decisionM)
  {
    n = doseDT$n[current] + csize
    x = doseDT$x[current] + rbinom(1, csize, prob = doseDT$tox[current])  # row
    y = doseDT$y[current] + rbinom(1, csize, prob = doseDT$eff[current])  # col
    
    doseDT$n[current] <- n 
    doseDT$x[current] <- x  # tox 
    doseDT$y[current] <- y  # eff
    
    decision <- decisionM[[round(n/csize)]][x + 1, y + 1]
    
    dN <- nrow(doseDT)
    above <- ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # NA or number
    below <- ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1])   # NA or number 
    # above <- ifelse(current == dN, dN, current + 1)
    # below <- ifelse(current == 1, 1, current - 1)
    
    switch(decision,
           EU = { 
             doseDT$keep[current] = 0                    # remove current dose 
             newdose <- above
           }, 
           DU_E = {
             doseDT$keep[current] = 0                    # remove current dose 
             newdose <- below
           }, 
           DU_T = {
             doseDT$keep[current:nrow(doseDT)] = 0                    # remove current and above doses 
             newdose <- below
           }, 
           E = {
             newdose <- ifelse(is.na(above), current, above)
           }, 
           S = {
             newdose <- current     # stay 
           },  
           D = {
             newdose <- ifelse(is.na(below), current, below)
           }, 
           ESD = {
             newdose <- NULL
             temp <- c(below, current, above)[!is.na(c(below, current, above))]   # find which are not NA (at least current)
             if(length(temp) == 1)
             {
               newdose <- temp
             } else
             {
               if(!is.na(above))   # above dose exists 
               { 
                 if(doseDT$n[above] == 0)   # above does has not been used 
                 {
                   newdose <- above   # use above 
                 } else    # above dose has been used 
                 {
                   qlist <- doseDT$y[temp]/doseDT$n[temp]        # select the ones with largest qhat 
                   left <- temp[qlist > (max(qlist) - 1e-4)]
                   newdose <- sample(left, 1)         # could be tie 
                 }
               } else   # above dose do not exist, i.e. below exist 
               {
                 qlist <- doseDT$y[temp]/doseDT$n[temp]    # still need to check 
                 left <- temp[qlist > (max(qlist) - 1e-4)]
                 newdose <- sample(left, 1)
               }
             }
           }, 
           ESDU = {                           # current should be removed due to futility rule 
             doseDT$keep[current] = 0
             newdose <- NULL
             temp <- c(below, above)[!is.na(c(below, above))]    # check is.na
             if(length(temp) == 0)
             {
               newdose <- NA
             } else if(length(temp) == 1)
             {
               newdose <- temp
             } else
             {
               qlist <- doseDT$y[temp]/doseDT$n[temp]    # still need to check 
               left <- temp[qlist > (max(qlist) - 1e-4)]
               newdose <- sample(left, 1)
             }
           }
    )
    return(list(doseDT = doseDT, newdose = newdose))
  }
  
  fun_BOINETdec <- function(n, pT, qE, lambda1, lambda2, eta1, eta_et, xi_et,
                            onlyDM = TRUE, alphaT = 1, betaT = 1, alphaE = 1, betaE = 1)
  {
    decisionM <- matrix("D", nrow = n+1, ncol = n+1)
    # decision matrix: row = toxicity, col = efficacy 
    notD <- sum((0:n)/n < lambda2)
    decisionM[1:notD, ] <- "S"
    notS <- sum((0:n)/n <= eta1)
    decisionM[1:notD, 1:notS] <- "ESD"
    E <- sum((0:n)/n <= lambda1)
    decisionM[1:E, 1:notS] <- "E"
    
    safetyPP <- NULL 
    futilPP <- NULL
    
    safetyPP <- sapply(0:n, FUN = function(x) 1-pbeta(0.4, x+alphaT, n-x+betaT)) %>% 
      matrix(., nrow = n+1, ncol = n+1) 
    # futility posterior prob
    futilPP <- sapply(0:n, FUN = function(x) 1-pbeta(0.2, x+alphaE, n-x+betaE)) %>% 
      matrix(., nrow = n+1, ncol = n+1) %>% t()
    
    if(is.null(xi_et))
    {
      decisionM[futilPP < 0.1 & decisionM %in% c("E")] <- "EU"      
      decisionM[futilPP < 0.1 & decisionM %in% c("ESD")] <- "ESDU"   
      decisionM[futilPP < 0.1 & decisionM %in% c("S", "D")] <- "DU_E"  # "DU_E"
    } else
    {
      decisionM[futilPP < xi_et & decisionM %in% c("E")] <- "EU"      
      decisionM[futilPP < xi_et & decisionM %in% c("ESD")] <- "ESDU"   
      decisionM[futilPP < xi_et & decisionM %in% c("S", "D")] <- "DU_E"  # "DU_E"
    }
    
    if(is.null(eta_et))
    {
      decisionM[safetyPP > 0.9] <- "DU_T" # "DU_T"  
    } else
    {
      decisionM[safetyPP > eta_et] <- "DU_T" # "DU_T"  
    }
    
    if(onlyDM == TRUE)
    {
      result <- decisionM
    } else
    {
      result <- list("DM" = decisionM, "SPP" = safetyPP, "FPP" = futilPP)
    }
    return(result)
  }
  
  fun_BOINET_cutpt <- function(nj, phi, delta, eta1 = NULL)
  {
    phi1 <- 0.1 * phi
    phi2 <- 1.4 * phi
    lambda1 <- seq(from = round(phi1, 2) + 0.01, to = round(phi, 2) - 0.01, by = 0.01)
    lambda2 <- seq(from = round(phi, 2) + 0.01, to = round(phi2, 2) - 0.01, by = 0.01)
    if(is.null(eta1))
    {
      delta1 <- 0.6 * delta
      eta1 <- seq(from = round(delta1, 2) + 0.01, to = round(delta, 2) - 0.01, by = 0.01)
      choice <- expand.grid(lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1) %>% as.matrix()
    } else
    {
      delta1 <- 0.6 * delta
      choice <- expand.grid(lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1) %>% as.matrix()
    }
    
    incorrect <- function(nj, delta, delta1, phi, phi1, phi2, param)
    {
      lambda1 = param[1]
      lambda2 = param[2]
      eta1 = param[3]
      1/6 * (pbinom(nj*lambda2-1, nj, phi1) * (1-pbinom(nj*eta1, nj, delta1)) + 
               2/3 * (pbinom(nj*lambda2-1, nj, phi1) - pbinom(nj*lambda1, nj, phi1)) * pbinom(nj*eta1, nj, delta1) + 
               1 - pbinom(nj*lambda2-1, nj, delta1) 
      ) + 
        1/6 * (pbinom(nj * lambda1, nj, phi1) * pbinom(nj*eta1, nj, delta) + 
                 2/3 * (pbinom(nj*lambda2 -1, nj, phi1) - pbinom(nj*lambda1, nj, phi1)) * pbinom(nj*eta1, nj, delta) + 
                 1 - pbinom(nj*lambda2-1, nj, delta)
        ) + 
        1/6 * (pbinom(nj*lambda1, nj, phi) * pbinom(nj*eta1, nj, delta) + 
                 2/3 * (pbinom(nj*lambda2 - 1, nj, phi) - pbinom(nj*lambda1, nj, phi)) * pbinom(nj*eta1, nj, delta) + 
                 1 - pbinom(nj*lambda2-1, nj, delta)
        ) + 
        1/6 * (pbinom(nj*lambda1, nj, phi2) * pbinom(nj*eta1, nj, delta1) + 
                 2/3 * (pbinom(nj*lambda2 - 1, nj, phi2) - pbinom(nj*lambda1, nj, phi2)) * pbinom(nj*eta1, nj, delta1) +
                 pbinom(nj*lambda2 - 1, nj, phi2) * (1- pbinom(nj*eta1, nj, delta1))
        ) + 
        1/6 * (pbinom(nj*lambda1, nj, phi2) * pbinom(nj*eta1, nj, delta) + 
                 2/3 * (pbinom(nj*lambda2 - 1, nj, phi2) - pbinom(nj*lambda1, nj, phi2)) * pbinom(nj*eta1, nj, delta) +
                 pbinom(nj*lambda2 - 1, nj, phi2) * (1- pbinom(nj*eta1, nj, delta))
        )
    }
    
    result <- sapply(1:nrow(choice), FUN = function(x) incorrect(nj, delta, delta1, phi, phi1, phi2, choice[x,])) %>% 
      which.min()
    return(choice[result,])
  }
  
  fun_BOINETOBD <- function(doseDT, pT)
  {
    DTremain <- doseDT[doseDT$n > 0, ]   # should call as used 
    rN <- nrow(DTremain)
    
    used <- which(doseDT$n > 0)     
    
    phat <- rep(pT, nrow(doseDT))
    phat[used] <- doseDT$x[used]/doseDT$n[used]
    ptilde <- pava(phat, w = 1/((doseDT$x + 0.05) * (doseDT$n - doseDT$x + 0.05)/((doseDT$n + 0.1)^2 * (doseDT$n + 0.1 + 1)))) + 0.001*(1:nrow(doseDT))
    
    q <- NULL
    z1 <- NULL
    z2 <- NULL
    for(i in used)
    {
      q <- c(q, rep(1, doseDT$y[i]), rep(0, doseDT$n[i] - doseDT$y[i]))
      z1 <- c(z1, rep(i, doseDT$n[i]))
    }
    
    if(min(ptilde) > pT)
      return(99)             # no does toxicity is acceptable 
    kset <- expand.grid(c(-2, -1, -0.5, 0, 0.5, 1, 2, 3), 
                        c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)) %>% as.matrix() 
    
    dev <- 1e8
    for(x in 1:nrow(kset))
    {
      if(kset[x,1] == 0)
      {
        x1 <- log(z1)
      } else
      {
        x1 <- (z1)^kset[x,1]
      }
      if(kset[x,2] == 0)
      {
        x2 <- log(z1)
      } else
      {
        x2 <- (z1)^kset[x,2]
      }
      lm1 <- glm(q ~ x1 + x2, family = "binomial")
      if(deviance(lm1) < dev)
      {
        dev <- deviance(lm1)
        est <- coef(lm1)
        bestK <- kset[x, ]
      }
    }
    if(bestK[1] == 0)
    {
      x1 <- log(1:nrow(doseDT))
    } else
    {
      x1 <- (1:nrow(doseDT))^bestK[1]
    }
    if(bestK[2] == 0)
    {
      x2 <- log(1:nrow(doseDT))
    } else
    {
      x2 <- (1:nrow(doseDT))^bestK[2]
    }
    logit <- est[1] + ifelse(is.na(est[2]), 0, est[2]) * x1 + ifelse(is.na(est[3]), 0, est[3]) * x2
    qtilde <- exp(logit)/(1 + exp(logit))
    
    J <- which.min(abs(ptilde - pT))
    OBD <- used[which.max(qtilde[used[used<= J]])]
    return(OBD)
    
  }
  
  dN <- length(plist)
  idlist <- 1:dN
  doseDT <- data.frame(id = 1:dN, tox = plist, eff = qlist, n = rep(0, dN), 
                       x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN)) %>% 
    arrange(tox) # sort by toxicity rate 
  
  if(is.null(BOINEThyper)) {
    hyperparam <- fun_BOINET_cutpt(csize * cN, phi, delta, eta1 = NULL)
  } else if(length(BOINEThyper)==1)
  {
    hyperparam <- fun_BOINET_cutpt(csize * cN, phi, delta, eta1 = BOINEThyper)
  } else
  {
    hyperparam <- BOINEThyper
  }
  
  decisionM <- list()
  for(i in 1:cN)
  {
    decisionM[[i]] <- fun_BOINETdec(n = i * csize, pT, qE, hyperparam[1], hyperparam[2], hyperparam[3], 
                                    eta_et, xi_et)
  }
  
  earlystop <- FALSE
  record <- rep(-1, cN)          # record dose selection
  for(i in 1:cN)
  {
    BOINET_onestep <- BOINET_one(doseDT, current, csize, decisionM)
    doseDT <- BOINET_onestep$doseDT
    record[i] <- current
    current <- BOINET_onestep$newdose   # next cohort dose 
    # check the next dose exists or not 
    if(is.na(current))
    {
      earlystop <- TRUE
      break
    } else if(doseDT$n[current] + csize > doselimit){  # check the next dose exceeds the dose limit or not 
      break   # not early stop 
    } 
  }
  
  if(earlystop == TRUE)
  {
    OBD <- 99   # 99: early stop
    rN <- 0
  } else
  {
    OBD <- fun_BOINETOBD(doseDT, pT)   # can be 99
    DTremain <- doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN <- nrow(DTremain)   # how many number of doses remain for selection
  }
  
  trueOBDone <- trueOBD[1]   # the best one 
  
  # with one OBD
  select_OBD <- ifelse(((trueOBDone %in% idlist) & (OBD %in% idlist) & (OBD == trueOBDone))  | ((!trueOBDone %in% idlist) & (!OBD %in% idlist)), 1, 0)
  num_at_OBD <- ifelse(trueOBDone %in% idlist, doseDT[trueOBDone, 4], NA)
  risk_allocate <- ifelse(trueOBDone %in% idlist, 
                          ifelse(doseDT[trueOBDone, 4] < csize * cN/5, 1, 0), NA)
  
  # with all possible case
  select_OBD_gen <- ifelse(((trueOBDone %in% idlist) & (OBD %in% idlist) & (OBD %in% trueOBD))  | ((!trueOBDone %in% idlist) & (!OBD %in% idlist)), 1, 0)
  num_at_OBD_gen <- ifelse(trueOBDone %in% idlist, sum(doseDT[trueOBD, 4]), NA)
  risk_allocate_gen <- ifelse(trueOBDone %in% idlist, 
                              ifelse(sum(doseDT[trueOBD, 4]) < csize * cN/5 * length(trueOBD), 1, 0), NA)
  
  num_overdose_OBD <- ifelse(trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  num_overdose_nOBD <- ifelse(!trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  
  
  
  # not include doseDT and record here
  return(data.frame(earlystop = earlystop, OBD = OBD, rN = rN, trueOBD = trueOBDone, 
                    select_OBD = select_OBD, num_at_OBD = num_at_OBD, 
                    select_OBD_gen = select_OBD_gen, num_at_OBD_gen = num_at_OBD_gen, 
                    num_overdose_OBD = num_overdose_OBD, 
                    num_overdose_nOBD = num_overdose_nOBD,
                    risk_allocate = risk_allocate,
                    risk_allocate_gen = risk_allocate_gen) %>% cbind(., t(data.frame(doseDT$n)))
  )
}




