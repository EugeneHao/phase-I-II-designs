# require(dplyr)
# require(Iso)
fun_UBI(c(0.1, 0.2, 0.4, 0.6), c(0.1, 0.7, 0.2, 0.1), trueOBD = 2,
        pT = 0.4, qE = 0.2, csize = 3, cN = 9)

fun_UBI <- function(plist, qlist, trueOBD, pT, qE, csize, cN, mindelta = 0.1,
                    size = 1000, p1 = 0.15, p2 = 0.4, q1 = 0.2, q2 = 0.6,
                    eta = 0.95, xi = 0.3, alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,   
                    current = 1, doselimit = Inf)
{
  
  # UBI for one cohort 
  UBI_one <- function(doseDT, current, pT, pE, qE, csize, decisionM)
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
  # 
  
  fun_decision <- function(n, pT, qE, eta = 0.95, xi = 0.3, 
                           alphaT = 1, betaT = 1, alphaE = 1, betaE = 1, onlyDM = TRUE)
  {
    phi1 = 0.6 * pT
    phi2 = 1.4 * pT
    
    lambda_e = log((1-phi1)/(1-pT))/log(pT * (1-phi1)/phi1/(1-pT))  # lower boundary
    lambda_d = log((1-pT)/(1-phi2))/log(phi2 * (1-pT)/pT/(1-phi2))  # upper boundary 
    
    qhat <- (0:n)/n
    qhat[qhat > 0.66] <- 0
    
    phat <- (0:n)/n
    U_E <-  matrix(rep(qhat, each = n+1), nrow = n+1)
    phat_1 <- c(rep(0, sum(phat <= lambda_e)), 
                phat[phat > lambda_e & phat < lambda_d], 
                rep(1, sum(phat >= lambda_d)))   # when qhat <= eff
    phat_2 <- c(rep(0, sum(phat <= 0.15)), 
                phat[phat > 0.15 & phat < lambda_d]/3, 
                rep(1, sum(phat >= lambda_d)))   # when qhat > eff
    U_T <- matrix(c(rep(phat_1, which.max(qhat)), 
                    rep(phat_2, n + 1 - which.max(qhat))), nrow = n+1)
    U <- U_E - 2 * U_T
    
    decisionM <- matrix("S", nrow = n + 1, ncol = n+1)
    decisionM[U >= 0 ] <- 'E'
    decisionM[U < -1/3] <- "D"
    
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
  
  fun_UBIOBD <- function(doseDT, pT, qE, p1, p2, q1, q2)
  {
    DTremain <- doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN <- nrow(DTremain)
    
    if(rN == 0)                  # no available dose has patients
    {
      return(0)           
    } else if(rN == 1)           # just return the only one dose  
    {
      opt <- 1
      return(DTremain$id[opt])    # select OBD 
    } else if(rN > 1)                       # at least two  doses left, need to check the utility 
    {
      used <- DTremain$id
      alphaTpost <- 0.5 + DTremain$x
      betaTpost <- 0.5 + DTremain$n - DTremain$x
      
      meanTpost <- alphaTpost/(alphaTpost + betaTpost)
      varTpost <- meanTpost * (1-meanTpost)/(1+alphaTpost+betaTpost)
      
      phat <- doseDT$x[used]/doseDT$n[used]
      ptilde <- pava(phat, w = 1/varTpost)  
      
      phat <- ptilde 
      qhat <- doseDT$y[DTremain$id]/doseDT$n[DTremain$id]
      
      phat[phat > p2] <- p2
      phat[phat < p1] <- p1
      qhat[qhat > q2] <- q2
      qhat[qhat < q1] <- q1
      u <- qhat - 2 * phat
      opt <- DTremain$id[which.max(u)]
      return(opt)
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
    decisionM[[i]] <- fun_decision(n = i * csize, pT, qE, eta, xi, alphaT, betaT, alphaE, betaE)
  }
  
  earlystop <- FALSE
  record <- rep(-1, cN)          # record dose selection
  for(i in 1:cN)
  {
    UBI_onestep <- UBI_one(doseDT, current, pT, pE, qE, csize, decisionM)
    doseDT <- UBI_onestep$doseDT
    record[i] <- current
    current <- UBI_onestep$newdose   # next cohort dose 
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
    OBD <- fun_UBIOBD(doseDT, pT, qE, p1, p2, q1, q2)
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




