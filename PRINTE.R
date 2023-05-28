# require(dplyr)
# require(Iso)
# fun_PRINTE(c(0.1, 0.2, 0.4, 0.6), c(0.1, 0.7, 0.2, 0.1), trueOBD = 2,
#            pT = 0.4, pE = 0.4, qE = 0.2, csize = 3, cN = 9)

fun_PRINTE <- function(plist, qlist, trueOBD, pT, pE, qE, csize, cN, mindelta = 0.1,
                       size = 1000, p1 = 0.2, p2 = 0.4, q1 = 0.2, q2 = 0.6, pgrad = 0.2, 
                       eta = 0.95, xi = 0.3, eps1 = 0.05, eps2 = 0.05, 
                       alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,
                       alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5, 
                       current = 1, doselimit = Inf)
{
  
  # PRINTE for one cohort 
  PRINTE_one <- function(doseDT, current, pT, pE, qE, csize, decisionM)
  {
    n = doseDT$n[current] + csize
    x = doseDT$x[current] + rbinom(1, csize, prob = doseDT$tox[current])  # row
    y = doseDT$y[current] + rbinom(1, csize, prob = doseDT$eff[current])  # col
    
    decision = decisionM[[round(n/csize)]][x + 1, y + 1]
    doseDT$n[current] <- n 
    doseDT$x[current] <- x
    doseDT$y[current] <- y
    
    # find the above and below closet available dose, can be NA if no dose found
    dN <- nrow(doseDT)
    above <- ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # NA or number
    below <- ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1])   # NA or number 
    # obtain dose for the next cohort, NA means early terminate 
    switch (decision,
            EU = { 
              doseDT$keep[current] = 0                    # remove current dose 
              newdose <- ifelse(is.na(above), below, above)
            }, 
            DU_E = {
              doseDT$keep[current] = 0                    # remove current dose 
              newdose <- below                            # terminate if no available dose below
            }, 
            DU_T = {
              doseDT$keep[current:dN] = 0                 # remove current and above doses 
              newdose <- below                            # terminate if no available dose below
            }, 
            E = {
              newdose <- ifelse(is.na(above), current, above)  # stay if no available dose above 
            }, 
            S = {
              newdose <- current     # stay 
            },  
            D = {
              newdose <- ifelse(is.na(below), current, below)  # stay if no available dose below 
            }
    )
    return(list(doseDT = doseDT, newdose = newdose))
  }
  
  fun_decision <- function(n, pT, pE, qE, eta = 0.95, xi = 0.3, eps1 = 0.05, eps2 = 0.05, 
                           alphaT = 1, betaT = 1, alphaE = 1, betaE = 1, onlyDM = TRUE)
  {
    if(pT %in% c(0.3, 0.4))
    {
      bL <- fun_bound(pT, pE, eps1, eps2)
    } else if(pT == 0.35)
    {
      bL <- list(bT_up = (1:10)/10,
                 bT_low = (0:9)/10,
                 bE_up = c(qE, pE, 0.6, 0.8, 1.0),
                 bE_low = c(0, qE, pE, 0.6, 0.8),
                 dT = c(rep("ND", 4), rep("D", 6)),
                 dE = c("E", "E", "S", "S", "S")
      )
    }
    dT <- sapply(0:n, FUN = function(x) 
      bL$dT[which.max((pbeta(bL$bT_up, x + alphaT, n-x+betaT) - pbeta(bL$bT_low, x+alphaT, n-x+betaT))/
                        (bL$bT_up - bL$bT_low))])
    dE <- sapply(0:n, FUN = function(x)
      bL$dE[which.max((pbeta(bL$bE_up, x + alphaE, n-x+betaE) - pbeta(bL$bE_low, x+alphaE, n-x+betaE))/
                        (bL$bE_up - bL$bE_low))])
    
    # decision matrix: row = toxicity, col = efficacy 
    decisionM <- rbind(rep(dE, each = sum(dT=="ND")) %>% matrix(., ncol = n+1),
                       matrix("D", nrow = sum(dT == "D"), ncol = n+1))
    # safety posterior prob
    safetyPP <- sapply(0:n, FUN = function(x) 1-pbeta(pT, x+alphaT, n-x+betaT)) %>% 
      matrix(., nrow = n+1, ncol = n+1) 
    # futility posterior prob
    futilPP <- sapply(0:n, FUN = function(x) 1-pbeta(qE, x+alphaE, n-x+betaE)) %>% 
      matrix(., nrow = n+1, ncol = n+1) %>% t()
    
    decisionM[futilPP < xi & decisionM == "E"] <- "EU"      
    decisionM[futilPP < xi & decisionM %in% c("S", "D")] <- "DU_E"
    decisionM[safetyPP > eta] <- "DU_T"  
    
    if(onlyDM == TRUE)
    {
      result <- decisionM
    } else
    {
      result <- list("DM" = decisionM, "SPP" = safetyPP, "FPP" = futilPP)
    }
    return(result)
  }
  
  fun_bound <-function(pT,pE,eps1,eps2){
    # toxicity
    bound_TE_up <- seq(from=pT-eps1, to=0.0001, by=-(eps1+eps2))
    bound_TE_low <- pmax(bound_TE_up-(eps1+eps2),0)
    
    bound_TD_low<- seq(from=pT+eps2, to=0.9999, by=(eps1+eps2))
    bound_TD_up <- pmin(bound_TD_low+(eps1+eps2),1)
    
    bound_T_low <- sort(c(bound_TE_low,pT-eps1,bound_TD_low))
    bound_T_up <- sort(c(bound_TE_up,pT+eps2,bound_TD_up))
    # efficacy 
    n.decimal<-nchar(pE,type="bytes")-2   #calculate number of decimal
    gcd.pE<- round(FRACTION::gcd(pE*10^n.decimal,(1-pE)*10^n.decimal)/(10^n.decimal), n.decimal)    #greatest common divisor
    while(gcd.pE >= 0.25)
    {
      gcd.pE <- gcd.pE-0.1
    }
    
    bound_E_up<-seq(from=pE,to=1,by=gcd.pE)
    if(pE > gcd.pE + 1e-4)
    {
      for(i in 1:floor(pE/gcd.pE))
      {
        if(bound_E_up[1] - gcd.pE > 1e-4)
        {
          tmp <- bound_E_up[1] - gcd.pE
          bound_E_up <- c(tmp, bound_E_up)
        }
      }
    }
    bound_E_low<- c(0, bound_E_up)[-(1+length(bound_E_up))]
    
    # toxicity decision 
    dT <- c(rep("ND", length(bound_TE_up) + 1),     # ND: not de-escalate
            rep("D", length(bound_TD_up)))          # D: de-escalate
    
    # efficacy decision 
    
    dE <- c(rep("E", which.min(abs(bound_E_up - pE))),       # E: escalate
            rep("S", length(bound_E_up) - which.min(abs(bound_E_up - pE))))  # S: stay
    
    return(list("bT_up"=bound_T_up,"bT_low"=bound_T_low,
                "bE_up"=bound_E_up,"bE_low"=bound_E_low, 
                "dT" = dT, "dE" = dE))
  }
  
  fun_OBD <- function(doseDT, pT, qE, mindelta = 0.1, size = 1000, p1, p2, q1, q2, pgrad = 0.2, 
                      alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5)
  {
    DTremain <- doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN <- nrow(DTremain)
    # 
    if(rN == 0)                  # no available dose has patients
    {
      return(0)           
    } else                        # at least one dose left, need to check the utility of each dose (PRINTE)
    {
      alphaTpost <- alphaTO + DTremain$x
      betaTpost <- betaTO + DTremain$n - DTremain$x
      alphaEpost <- alphaEO + DTremain$y
      betaEpost <- betaEO + DTremain$n - DTremain$y
      
      # posterior mean and var 
      meanTpost <- alphaTpost/(alphaTpost + betaTpost)
      varTpost <- meanTpost * (1-meanTpost)/(1+alphaTpost+betaTpost)
      
      if(rN > 1)
      { 
        D <- length(alphaTpost)
        # sample efficacy (all)
        sampleE <- sapply(1:D, FUN = function(x) rbeta(size, alphaEpost[x], betaEpost[x])) # size * D
        
        utilE <- sampleE 
        utilE[sampleE >= q2] <- 1
        utilE[sampleE <= q1] <- 0
        utilE[sampleE < q2 & sampleE > q1] <- (sampleE[sampleE < q2 & sampleE > q1] - q1)/(q2-q1)
        # sample toxicity (all)
        sampleT <- sapply(1:D, FUN = function(x) rbeta(size, alphaTpost[x], betaTpost[x])) %>% 
          apply(., MARGIN = 1, FUN = pava, w = 1/varTpost) %>% 
          t() # size * D 
        utilT <- sampleT
        utilT[sampleT >= p2] <- 0
        utilT[sampleT <= p1] <- 1
        utilT[sampleT < p2 & sampleT > p1] <- 1- (sampleT[sampleT<p2 & sampleT>p1] - p1)/(p2-p1)
        
        EU <- colMeans(utilE, na.rm = T) * colMeans(utilT, na.rm = T)
        
        # use all dose 
        opt <- which.max(EU)
        pin <- mean((sampleT[,opt] <= pT) * (sampleE[,opt] >= qE + mindelta), na.rm = T)
      } else
      {
        opt <- 1
        
        # use all 
        sampleE <- rbeta(size, alphaEpost, betaEpost)
        sampleT <- rbeta(size, alphaTpost, betaTpost)
        pin <- mean((sampleT <= pT) * (sampleE >= (qE + mindelta)), na.rm = T)
      }
      
      if(pin < pgrad)   # only compare with the threshold if we use PRINTE
      {
        return(0)     # no dose has utility larger than p_graduate 
      } else
      {
        return(DTremain$id[opt])    # select OBD 
      }
    }
  }
  
  dN <- length(plist)
  idlist <- 1:dN
  doseDT <- data.frame(id = 1:dN, tox = plist, eff = qlist, n = rep(0, dN), 
                       x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN)) %>% 
    arrange(tox) # sort by toxicity rate 
  
  decisionM <- list()
  for(i in 1:cN)
  {
    decisionM[[i]] <- fun_decision(n = i * csize, pT, pE, qE, eta, xi, eps1, eps2, 
                                   alphaT, betaT, alphaE, betaE)
  }
  
  earlystop <- FALSE
  record <- rep(-1, cN)          # record dose selection
  for(i in 1:cN)
  {
    PRINTE_onestep <- PRINTE_one(doseDT, current, pT, pE, qE, csize, decisionM)
    doseDT <- PRINTE_onestep$doseDT
    record[i] <- current
    current <- PRINTE_onestep$newdose   # next cohort dose 
    # check the next dose exists or not 
    if(is.na(current))
    {
      earlystop <- TRUE
      break
    } else if(doseDT$n[current] + csize > doselimit){  # check the next dose exceeds the dose limit or not 
      break   # not early stop 
    }
  }
  
  # select optimal dose if more than two doses left after the trial 
  if(earlystop == FALSE)
  {
    OBD <- fun_OBD(doseDT, pT, qE, mindelta, size, p1, p2, q1, q2, pgrad, 
                   alphaTO, betaTO, alphaEO, betaEO)
    DTremain <- doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN <- nrow(DTremain)   # how many number of doses remain for selection
  } else
  {
    OBD <- 99   # 99: early stop
    rN <- 0
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



