# require(dplyr)
# require(Iso)
# fun_uTPI(c(0.1, 0.2, 0.4, 0.6), c(0.1, 0.7, 0.2, 0.1), trueOBD = 2, pT = 0.4, qE = 0.2, 
#          eps = 0.1, del = 0.1, w1 = 0.7, w4 = 0.3, nstar = 9, 
#          cutoff_tox = 0.95, cutoff_eff = 0.9, csize = 3, cN = 9)

fun_uTPI = function(plist, qlist, trueOBD, pT, qE, eps, del, w1, w4, nstar,
                    cutoff_tox = 0.95, cutoff_eff = 0.9, csize, cN, current = 1, doselimit = Inf,
                    alphaE = 1, betaE = 1, alphaT = 1, betaT = 1) 
{

  uTPI_one = function(doseDT, current, pT, qE, del, csize, decisionM, kU, w1, w4, nstar,
                      alphaE, betaE, cutoff_eff) {
    
    n = doseDT$n[current] + csize # total number of subjects that have been treated at the current dose level
    x = doseDT$x[current] + rbinom(1, csize, prob = doseDT$tox[current])  # col: total number of DLT events at the current dose level
    y = doseDT$y[current] + rbinom(1, csize, prob = doseDT$eff[current])  # row: total number of responses at the current dose level
    
    phat = x/n # cumulative phat at the current dose level
    qhat = y/n # cumulative qhat at the current dose level
    
    doseDT$n[current] = n  # update total number of patients have been treated at the current dose level
    doseDT$x[current] = x  # update total number of DLT events 
    doseDT$y[current] = y  # update total number of eff events
    
    dN = nrow(doseDT) # number of dose levels
    above = ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # the next available higher dose, can be NA if no available dose
    below = ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1]) # the next available lower dose, can be NA if no available dose
    
    decision = decisionM[[round(n/csize)]][y+1, x+1] # col: tox; row: res
    
    u.l = seq(0, (1/del-1)/(1/del), by=del) # lower bounds of desirability intervals
    u.u = seq(del, 1, by=del) # upper bounds of desirability intervals
    
    #-----------------compute kU's for admissible dose set-------------------#
    
    # calculate futility posterior prob at current dose
    futilPP = pbeta(qE, y+alphaE, n-y+betaE)
    
    # posterior desirability for current dose
    if(futilPP > cutoff_eff) {
      kU_current = -Inf # current dose will not be considered due to futility
      doseDT$keep[current] = 0 # eliminate current dose level due to futility
    } else {
      # compute sum(uhat)
      if(n >= nstar) {
        uhat = (w1*qhat+(1-phat)*w4)*n 
      } else {
        uhat = (w1*qhat+w4)*n
      }
      
      # update kU for the current dose
      kU[current] = which.max(pbeta(u.u,1+uhat,1+n-uhat)-pbeta(u.l,1+uhat,1+n-uhat))
      kU[current] = kU[current] + (1-pbeta(u.u[kU[current]],1+uhat,1+n-uhat)) # Avoid ties
      kU_current = kU[current]
    }
    
    # get kU for the dose below
    if(is.na(below)) {
      kU_below = -Inf
    } else {
      kU_below = kU[below]
    }
    
    # get kU for the dose above
    if(is.na(above)) {
      kU_above = -Inf
    } else {
      kU_above = kU[above]
    }
    #------------------------------------------------------------------------#
    
    # obtain dose for the next cohort, NA means early terminate 
    switch (decision,
            DU_E = {
              doseDT$keep[current] = 0    # remove current dose 
              newdose = below             # terminate if no available dose below
            }, 
            DU_T = {
              doseDT$keep[current:dN] = 0 # remove current and above doses 
              newdose = below             # terminate if no available dose below
            }, 
            D = {
              # stay if no available dose below 
              newdose = ifelse(is.na(below), current, below)  
            },
            AL = {
              # evaluate kU's for each dose in the local admissible dose set
              post = c(kU_below, kU_current, kU_above) + runif(3)*(10)^(-15)
              select = c(below, current, above)
              newdose = select[which.max(post)]
            },
            AS = {
              # evaluate kU's for each dose in the local admissible dose set
              post = c(kU_below, kU_current) + runif(2)*(10)^(-15)
              select = c(below, current)
              newdose = select[which.max(post)]
            }
    )
    
    return(list(doseDT = doseDT, newdose = newdose, kU = kU))
  }
  
  fun_uTPI_dec = function(n, pT, qE, eps = 0.1, nstar, cutoff_tox = 0.95, cutoff_eff = 0.90,
                          alphaT = 1, betaT = 1, alphaE = 1, betaE = 1, onlyDM = T) {
    
    # Identify k*
    kstar = min(floor(pT*10) + 1, 1/eps)
    
    x = 0:n # number of DLTs at the given dose j
    y = 0:n # number of responses at the given dose j
    
    tox.l = seq(0, (1/eps-1)/(1/eps), by=eps) # lower bound of toxicity rate intervals
    tox.u = seq(eps, 1, by=eps) # upper bound of toxicity rate intervals
    
    
    # Calculate strongest toxicity index for all possible number of DLTs
    kT = sapply(0:n, FUN = function(x) which.max((pbeta(tox.u, x + alphaT, n-x+betaT) - pbeta(tox.l, x+alphaT, n-x+betaT))))
    kTM = matrix(rep(kT, each = n+1), nrow = n+1) # kT matrix
    
    
    #------------Decision-making steps--------------#
    
    # Decision table initialization
    decisionM = matrix("TBD", nrow = n+1, ncol = n+1) 
    
    # Step (a): De-escalate if kT > k*
    decisionM[kTM > kstar] = "D" 
    
    # Step (b): Consider large/small admissible dose set
    if(n < nstar) {
      decisionM[kTM <= kstar] = "AL" # Large admissible dose set {j-1, j, j+1}
    } else {
      decisionM[kTM < kstar] = "AL" 
      decisionM[kTM == kstar] = "AS" # Small admissible dose set {j-1, j}
    }
    
    
    # Dose elimination due to futility control
    futilPP = sapply(0:n, FUN = function(x) pbeta(qE, x+alphaE, n-x+betaE)) %>% 
      matrix(., nrow = n+1, ncol = n+1)
    decisionM[futilPP > cutoff_eff & decisionM == "D"] = "DU_E"
    
    # Dose elimination due to safety control
    safetyPP = sapply(0:n, FUN = function(x) 1-pbeta(pT, x+alphaT, n-x+betaT)) %>% 
      matrix(., nrow = n+1, ncol = n+1) %>% t()
    decisionM[safetyPP > cutoff_tox] = "DU_T"
    
    
    
    if(onlyDM == TRUE){
      result = decisionM
    } else {
      result = list("DM" = decisionM, "SPP" = safetyPP, "FPP" = futilPP)
    }
    
    return(result)
  }
  
  fun_uTPIOBD = function(doseDT, pT, w1, w4){
    
    # Obtain the isotonically transformed values for toxicity probs
    PT = (doseDT$x + 0.05)/(doseDT$n + 0.1)
    ptilde = pava(PT, doseDT$n+0.1)+0.001*seq(1,nrow(doseDT))
    
    # Obtain the unimodal-isotonically transformed values for efficacy probs
    qtilde = peestimate(doseDT$y, doseDT$n)
    
    # Compute utility
    U = 100*(w1*qtilde + w4*(1-ptilde))
    U[which(doseDT$n == 0)] = -100
    U[which(doseDT$keep == 0)] = -100
    
    # Estimate MTD and OBD
    MTD = which.min(abs(ptilde-pT))
    OBD = which.max(U[1:MTD])
    
    # Evaluation of the selected OBD
    if(U[OBD] < 0) { # Consider the situation where all Utilities are < 0
      OBD = 99 # if all utilities are below 0, we conclude that no OBD is selected
    }
    
    return(OBD)
  }
  
  peestimate = function(yE,n) {
    ndose = length(yE)
    lik = rep(0,ndose)
    pe = (yE+0.05)/(n+0.1)
    p.e = matrix(NA,ndose,ndose)
    
    # unimodal isotonic regression
    for (i in 1:ndose) {
      
      if (i==1) {x=seq(ndose,1,by=-1)} else {x = c(1:(i-1),seq(ndose,i))}
      
      p.e[i,] = ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
      lik[i] = prod(dbinom(yE,n,p.e[i,]))		
    }
    lik = lik/sum(lik)
    pe = t(p.e)%*%lik+0.01*seq(1,ndose)
    
    return(pe)
  }
  
  decisionM = list()
  for (i in 1:cN) {
    decisionM[[i]] = fun_uTPI_dec(n=i*csize, pT, qE, eps, nstar, cutoff_tox, cutoff_eff)
  }
  
  dN = length(plist) # the total number of dose levels
  idlist = 1:dN # id index for each dose level
  
  # create a data frame to record observed data at each dose level
  doseDT = data.frame(id = 1:dN, tox = plist, eff = qlist, n = rep(0, dN), 
                      x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN)) %>% 
    arrange(tox) # sort by toxicity rate 
  
  earlystop = F
  record = rep(-1, cN)  # record dose selection
  
  # assign initial posterior desirability to each dose level
  kU = rep((w1*0.5+w4)*(1/del), dN)
  
  for (i in 1:cN) {
    uTPI_onestep = uTPI_one(doseDT, current, pT, qE, del, csize, decisionM, 
                            kU, w1, w4, nstar, alphaE, betaE, cutoff_eff)
    doseDT = uTPI_onestep$doseDT
    record[i] = current
    current = uTPI_onestep$newdose   # next cohort dose
    kU = uTPI_onestep$kU # updated kU for all doses
    
    # check whether the next dose exists or not 
    if(is.na(current)) {
      earlystop = T
      break
    } else if(doseDT$n[current] + csize > doselimit) {  # check whether the next dose exceeds the dose limit or not 
      break   # not early stop 
    } 
  }
  
  # OBD selection
  if(earlystop == F) {
    OBD = fun_uTPIOBD(doseDT, pT, w1, w4)
    DTremain = doseDT[doseDT$keep == 1 & doseDT$n > 0, ]  # need to select doses applied to patients
    rN = nrow(DTremain)   # how many number of doses remain for selection
  } else {
    OBD = 99   # 99: early stop
    rN = 0
  }
  
  trueOBDone = trueOBD[1]   # the best one 
  
  # with one OBD
  select_OBD <- ifelse(((trueOBDone %in% idlist) & (OBD %in% idlist) & (OBD == trueOBDone))  | ((!trueOBDone %in% idlist) & (!OBD %in% idlist)), 1, 0)
  num_at_OBD <- ifelse(trueOBDone %in% idlist, doseDT[trueOBDone, 4], NA)
  risk_allocate <- ifelse(trueOBDone %in% idlist, 
                          ifelse(doseDT[trueOBDone, 4] < csize * cN/5, 1, 0), NA)
  
  num_overdose_OBD <- ifelse(trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  num_overdose_nOBD <- ifelse(!trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT[which(plist > (pT+0.1)), 4]), 0), NA)
  
  
  # not include doseDT and record here
  return(data.frame(earlystop = earlystop, OBD = OBD, rN = rN, trueOBD = trueOBDone, 
                    select_OBD = select_OBD, num_at_OBD = num_at_OBD, 
                    num_overdose_OBD = num_overdose_OBD, 
                    num_overdose_nOBD = num_overdose_nOBD,
                    risk_allocate = risk_allocate) %>% cbind(., t(data.frame(doseDT$n)))
  )
} 
