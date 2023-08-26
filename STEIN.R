# require(dplyr)
# require(Iso)
# fun_STEIN(c(0.1, 0.2, 0.4, 0.6), c(0.1, 0.7, 0.2, 0.1), trueOBD = 2, pT = 0.4, psi1 = 0.3, 
#           psi2 = 0.8, w1 = 0.33, w2 = 1.09, cutoff_tox = 0.95, cutoff_eff = 0.9, csize = 3, cN = 9)

fun_STEIN = function(plist, qlist, trueOBD, pT, psi1, psi2, w1, w2, cutoff_tox, cutoff_eff, 
                     csize, cN, current = 1, doselimit = Inf,
                     futil_threshold = 0.2, alphaT = 1, betaT = 1, alphaE = 1, betaE = 1) 
{

  STEIN_one = function(doseDT, current, pT, psi, phi_L, phi_U, 
                       csize, decisionM, futil_threshold,
                       alphaE=1, betaE=1, cutoff_eff) {
    
    n = doseDT$n[current] + csize # total number of subjects that have been treated at the current dose level
    x = doseDT$x[current] + rbinom(1, csize, prob = doseDT$tox[current])  # col: total number of DLT events at the current dose level
    y = doseDT$y[current] + rbinom(1, csize, prob = doseDT$eff[current])  # row: total number of responses at the current dose level
    
    doseDT$n[current] = n  # update total number of patients have been treated at the current dose level
    doseDT$x[current] = x  # update total number of DLT events 
    doseDT$y[current] = y  # update total number of eff events
    
    phat = x/n # cumulative phat at the current dose level
    
    dN = nrow(doseDT) # number of dose levels
    above = ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # the next available higher dose, can be NA if no available dose
    below = ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1]) # the next available lower dose, can be NA if no available dose
    
    decision = decisionM[[round(n/csize)]][y+1, x+1] # col: tox; row: res
    
    # obtain dose for the next cohort, NA means early terminate 
    switch (decision,
            DU_E = {
              doseDT$keep[current] = 0                    # remove current dose 
              newdose = below                             # terminate if no available dose below
            }, 
            DU_T = {
              doseDT$keep[current:dN] = 0                 # remove current and above doses 
              newdose = below                             # terminate if no available dose below
            }, 
            S = {
              newdose = current     # stay 
            },  
            D = {
              newdose = ifelse(is.na(below), current, below)  # stay if no available dose below 
            },
            
            TBD = {
              # calculate futility posterior prob at current dose
              futilPP = pbeta(futil_threshold, y+alphaE, n-y+betaE)
              
              # posterior prob for current dose
              if(futilPP > cutoff_eff) {
                post_current = -Inf # current dose will not be considered due to futility
                doseDT$keep[current] = 0  # remove current dose due to futility
              } else {
                post_current = 1 - pbeta(psi, y+1, n-y+1) 
              }
              
              # posterior prob for the dose below
              if(is.na(below)) {
                post_below = -Inf
              } else {
                post_below = 1 - pbeta(psi, doseDT$y[below]+1, doseDT$n[below]-doseDT$y[below]+1)
              }
              
              # posterior prob for the dose above
              if(is.na(above)) {
                post_above = -Inf
              } else {
                post_above = 1 - pbeta(psi, doseDT$y[above]+1, doseDT$n[above]-doseDT$y[above]+1)
              }
              
              # evaluate posterior probs for each dose in the local admissible dose set
              if(phat <= phi_L) { # yellow region, we consider three possible doses {j-1, j, j+1}
                post = c(post_below, post_current, post_above) + c(0, 0.001, 0.002)
                select = c(below, current, above)
                newdose = select[which.max(post)]
              } else if (phat < phi_U & phat > phi_L) { # orange region, we only consider two possible doses {j-1, j}
                post = c(post_below, post_current) + c(0, 0.001)
                select = c(below, current)
                newdose = select[which.max(post)]
              }
            }
            
    )
    return(list(doseDT = doseDT, newdose = newdose))
  }
  
  fun_STEIN_dec = function(n, pT, phi1, phi2, psi1, psi2, cutoff_tox = 0.95, cutoff_eff = 0.95, 
                           alphaT = 1, betaT = 1, alphaE = 1, betaE = 1, 
                           futil_threshold = 0.2, onlyDM = T) {
    
    # Calculate phi_L and phi_U based on pre-specified phi1, phi2 and pT (phi0)
    phi_L = log((1-phi1)/(1-pT))/log((pT*(1-phi1))/(phi1*(1-pT))) # lower bound
    phi_U = log((1-pT)/(1-phi2))/log((phi2*(1-pT))/(pT*(1-phi2))) # upper bound
    
    # Calculate psi based on pre-specified psi1 and psi2
    psi = log((1-psi1)/(1-psi2))/log((psi2*(1-psi1))/(psi1*(1-psi2)))
    
    qhat = (0:n)/n # all possible observed efficacy probability
    phat = (0:n)/n # all possible observed toxicity probability
    
    qhatM = t(matrix(rep(qhat, each = n+1), nrow = n+1)) # qhat matrix
    phatM = matrix(rep(phat, each = n+1), nrow = n+1) # phat matrix
    
    #------------Decision-making steps--------------#
    
    # Decision table initialization
    decisionM = matrix("TBD", nrow = n + 1, ncol = n+1) 
    
    # Step (a): De-escalate if phat >= phi_U
    decisionM[phatM >= phi_U] = "D" 
    
    # Step (b): Retain current dose level if phat < phi_U and qhat > psi
    decisionM[phatM < phi_U & qhatM >= psi] = "S" 
    
    # Safety posterior prob
    safetyPP = sapply(0:n, FUN = function(x) 1-pbeta(pT, x+alphaT, n-x+betaT)) %>% 
      matrix(., nrow = n+1, ncol = n+1) %>% t()
    # Futility posterior prob
    futilPP = sapply(0:n, FUN = function(x) pbeta(futil_threshold, x+alphaE, n-x+betaE)) %>% 
      matrix(., nrow = n+1, ncol = n+1)
    
    # Apply elimination rules due to safety and futility concerns
    decisionM[futilPP > cutoff_eff & decisionM %in% c("S", "D")] = "DU_E"      
    decisionM[safetyPP > cutoff_tox] = "DU_T"
    #-----------------------------------------------#
    
    if(onlyDM == TRUE){
      result = decisionM
    } else {
      result = list("DM" = decisionM, "SPP" = safetyPP, "FPP" = futilPP)
    }
    return(result)
  }
  
  fun_STEINOBD = function(doseDT, pT, w1 = 0.33, w2 = 1.09) {
    
    # Obtain the isotonically transformed values for toxicity probs
    PT = (doseDT$x + 0.05)/(doseDT$n + 0.1)
    ptilde = pava(PT, doseDT$n+0.1)+0.001*seq(1,nrow(doseDT))
    ptilde[which(doseDT$n == 0)] = 20
    ptilde[which(doseDT$keep == 0)] = 20
    
    # Obtain the unimodal-isotonically transformed values for efficacy probs
    qtilde = peestimate(doseDT$y, doseDT$n)
    qtilde[which(doseDT$keep == 0)] = 0
    qtilde[which(doseDT$n == 0)] = 0
    
    # Utility computation
    u = qtilde - w1*ptilde - w2*ptilde*(ptilde > pT)
    opt = doseDT$id[which.max(u)]
    return(opt)
    
  }
  
  peestimate = function(yE,n) {
    ndose = length(yE)
    lik = rep(0,ndose)
    pe = (yE+0.05)/(n+0.1)
    p.e = matrix(NA,ndose,ndose)
    
    for (i in 1:ndose) {
      
      if (i==1) {x=seq(ndose,1,by=-1)} else {x = c(1:(i-1),seq(ndose,i))}
      
      p.e[i,] = ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
      lik[i] = prod(dbinom(yE,n,p.e[i,]))		
    }
    lik = lik/sum(lik)
    pe = t(p.e)%*%lik+0.01*seq(1,ndose)
    
    return(pe)
  }
  
  
  phi1 = 0.75*pT
  phi2 = 1.25*pT
  psi = log((1-psi1)/(1-psi2))/log((psi2*(1-psi1))/(psi1*(1-psi2)))
  phi_L = log((1-phi1)/(1-pT))/log((pT*(1-phi1))/(phi1*(1-pT))) # lower bound
  phi_U = log((1-pT)/(1-phi2))/log((phi2*(1-pT))/(pT*(1-phi2))) # upper bound
  
  # Get decision matrix for each possible sample size n
  decisionM = list()
  for (i in 1:cN) {
    decisionM[[i]] <- fun_STEIN_dec(n = i*csize, pT, phi1, phi2, psi1, psi2, cutoff_tox, cutoff_eff)
  }
  
  dN = length(plist) # the total number of dose levels
  idlist = 1:dN # id index for each dose level
  
  doseDT = data.frame(id = 1:dN, tox = plist, eff = qlist, n = rep(0, dN), 
                      x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN)) %>% 
    arrange(tox) # sort by toxicity rate 
  
  # u = doseDT$eff - w1*doseDT$tox - w2*doseDT$tox*(doseDT$tox > pT)
  # doseDT$u = u
  earlystop = FALSE
  record = rep(-1, cN)  # record dose selection
  
  for (i in 1:cN) {
    STEIN_onestep = STEIN_one(doseDT, current, pT, psi, phi_L, phi_U, csize, decisionM, futil_threshold,
                              alphaE, betaE, cutoff_eff = cutoff_eff)
    doseDT = STEIN_onestep$doseDT
    record[i] = current
    current = STEIN_onestep$newdose   # next cohort dose 
    # check whether the next dose exists or not 
    if(is.na(current)) {
      earlystop = TRUE
      OBD = 99
      break
    } else if(doseDT$n[current] + csize > doselimit) {  # check whether the next dose exceeds the dose limit or not 
      break   # not early stop 
    } 
  }
  
  # OBD selection
  if(earlystop == FALSE) {
    OBD = fun_STEINOBD(doseDT, pT, w1, w2)
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
