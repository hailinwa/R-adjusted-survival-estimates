diff.adj.surv.v3 <- function(indata, coxf, seednum=1986, n.sim=2000, cb.alpha=.05, starttime=NULL, endtime=NULL){

  #include libraries
  libs <- c("survival", "plyr", "tictoc")
  lapply(libs, require, character.only = TRUE)
  
tic("Total time")
####data manipulation####
  #check abnormal arguments
  # check "Surv" object must be there
  if (cb.alpha <= 0 | cb.alpha >= 1) stop("alpha must be within (0,1)")
  
tic("check data")
  #extract strata variable from formula
  scoxf <- deparse(substitute(coxf[[3]]), width.cutoff = 500)
  strata <- gsub("strata\\((.*?)\\)", "\\1", regmatches(scoxf, regexpr("strata\\(.*?\\)", scoxf)))
  #extract time / event variable from formula
  terms.y <- function (x){
    if (class(x) == "formula")
      c(terms.y(x[[2]]))
    else if (class(x) == "call") {
      if (x[[1]] == as.name("Surv") )
        unlist(lapply(x[-1], terms.y))
    }
    else (deparse(x))
  }
  tm.evt <- terms.y(coxf)

  
  if (length(tm.evt) == 2){#regular cox
    time <- tm.evt[1]
    event <- tm.evt[2]
    indata <- subset(indata, indata[[time]] > 0 & indata[[event]] >= 0)
    }
  else if (length(tm.evt) == 3){#left truncated cox
    ltime <- tm.evt[1]
    time <- tm.evt[2]
    event <- tm.evt[3]
    indata <- subset(indata, indata[[time]] > 0 & indata[[event]] >= 0 & indata[[ltime]] >= 0)
    }
  else {stop("Incorrect number of time and event parameter")}  

toc()
tic("coxph")
  # fit stratified cox model
  cox <- coxph(coxf, data = indata, ties = "breslow", model = TRUE)
  coxd <- coxph.detail(cox, riskmat = TRUE)
toc()

####general matrix / vector set up####
tic("matrix setup")
      nt <- NROW(coxd$time)             #number of distinct time points by strata (coxph output)
      nz <- NCOL(coxd$x)                #number of risk factor
      nl <- NROW(coxd$x)                #number of cases
      ns <- NROW(coxd$strata)           #number of strata
      nst <- sum(1 : ns)
      sc <- t(combn(ns, 2))
      nsc <- NROW(sc)                   #number of strata combination
      
  #create pooled time list
      bytime <- sort(unique(indata[[time]][which(indata[[event]] == 1)]))
      nt.unique <- length(bytime)
  #obtain sorted value of strata variable
      bystrata <- sort(unique(indata[[strata]]))
  #create mapping matrix pointing to the location of adj.surv and other computed matrices
      byloc <- matrix(, nrow = nt.unique, ncol = ns)
      for (s in 1 : ns){
        for (t in 1 : nt.unique){
          cur.time <- suppressWarnings(max(indata[[time]][which(indata[[strata]] == bystrata[s] & indata[[time]] <= bytime[t] & indata[[event]] == 1)]))
          byloc[t, s] <- suppressWarnings(min(which(indata[[strata]] == bystrata[s] & indata[[time]] == cur.time & indata[[event]] == 1)))
        }
      }
  #order coxd$x by case (to match with coxd$riskmat and case order from original data)
      if (nz == 1){coxd$x <- matrix(coxd$x[order(as.numeric(names(coxd$x)))], nrow = nl, ncol = 1)}
      else {coxd$x <- coxd$x[order(as.numeric(row.names(coxd$x))), ]}
  #histmat & riskmat (nl x nl, original case order)
      histmat <- matrix(0, nrow = nl, ncol = nl)
      riskmat <- matrix(0, nrow = nl, ncol = nl)
      for (l in 1 : nl){
        histmat[l, which(indata[[strata]] == indata[[strata]][l] & indata[[time]] <= indata[[time]][l])] = 1
        if (exists("ltime") == TRUE)
          {riskmat[l, which(indata[[strata]] == indata[[strata]][l] & indata[[time]] >= indata[[time]][l] & indata[[time]][l] > indata[[ltime]])] = 1}
        else if (exists("ltime") == FALSE)
          {riskmat[l, which(indata[[strata]] == indata[[strata]][l] & indata[[time]] >= indata[[time]][l])] = 1}
      }
toc()

tic("DAS")
####Adjusted survival w/ variance####
  #zbeta (by case): nlx1
      zbeta <- coxd$x %*% cox$coefficients
  #weighted zbeta, combination of z and calculated zbeta: nl.w x (nz+zbeta+count)
      zbeta.comb <- plyr::count(cbind(coxd$x, coxd$x %*% cox$coefficients))
  #number of combination of risk factor
      nl.w <- NROW(zbeta.comb)
  #z (weighted): nl.w x nz
      if (nz == 1){z.w <- zbeta.comb}
      else {z.w <- zbeta.comb[, 1 : nz]}
  #zbeta (weighted): nl.w x 1
      zbeta.w <- zbeta.comb[, nz + 1]
  #weight (number in each z combination): nl.w x 1
      weight <- zbeta.comb[, nz + 2]

  #s0beta (by case): nl x 1
      s0beta <- riskmat %*% exp(zbeta)
  #s1beta (by case): nl x nz
      s1beta <- riskmat %*% (coxd$x * as.vector(exp(zbeta)))
  #ebeta (by case): nl x nz
      ebeta <- s1beta / as.vector(s0beta)

  #baseline hazard: nl x 1
      h0 <- indata[[event]] / s0beta
  #baseline cumulative hazard (Breslow's estimator): nl x 1
      ch0 <- histmat %*% h0
  #Predicted survival probability by Z: nlxnl.w (weighted)
      s.w <- exp(-ch0 %*% t(exp(zbeta.w)))
toc()
tic("V1+V2")
  #variance of adj. survival
    #V1
      v1core <- histmat %*% (indata[[event]] / (s0beta ^ 2))
      v1 <- (s.w %*% (weight * exp(zbeta.w))) ^ 2 * v1core / (nl ^ 2)
    #V2
      ssh <- matrix(0, nrow = nl, ncol = nz)
      for (i in 1 : nz){
        z.w.rep <- matrix(rep(t(z.w[ ,i]), each = nl), nrow = nl)
        ebeta.rep <- matrix(rep(ebeta[ ,i], each = nl.w), ncol = nl.w, byrow = TRUE)
        ssh[, i] = (s.w * ( histmat %*% ((z.w.rep - ebeta.rep) * (h0 %*% t(exp(zbeta.w)))))) %*% weight
      }
      v2 <- diag(ssh %*% cox$var %*% t(ssh)) / (nl ^ 2)
toc()
    #sum
      adjvar <- v1 + v2
      adjse <- sqrt(adjvar)
    #Directly adjusted survival prob
      adjsurv <- (s.w %*% weight) / nl
      
tic("CB sim setup")
####Confidence band for adjusted survival####
      adj.s <- as.matrix(cbind(indata[[strata]], indata[[time]], indata[[event]]))
      adj.s <- subset(adj.s, adj.s[, 3] == 1)
      adj.s <- plyr::count(adj.s[order(adj.s[, 1], adj.s[, 2]), ])
      tmin.bystrata <- NULL
      tmax.bystrata <- NULL
      for (s in 1 : ns){
        tmin.bystrata[s] <- min(coxd$time[which(adj.s[1] == bystrata[s] & coxd$nevent > 0)])
        tmax.bystrata[s] <- max(coxd$time[which(adj.s[1] == bystrata[s] & coxd$nrisk >= 15)])
      }
      t1 <- ifelse(is.null(starttime), max(tmin.bystrata), starttime)
      t2 <- ifelse(is.null(endtime), min(tmax.bystrata), endtime)
      if (t1 >= t2) stop("Start time >= stop time for CB construction")
      
      set.seed(seednum)
      gcase <- matrix(rnorm(as.integer(nl * n.sim), mean = 0, sd = 1), nrow = nl, ncol = n.sim)
      gevent <- gcase * as.vector(indata[[event]])
      h0.sim <- gevent / as.vector(s0beta)
      ch0.sim <- histmat %*% h0.sim
      w1.sim <- as.vector(-(s.w %*% (weight * exp(zbeta.w)))) * ch0.sim / nl
      betadiff.sim <- cox$var %*% (t(coxd$x - ebeta) %*% gevent)
      w2.sim <- -(ssh %*% betadiff.sim) / nl
toc()
tic("DAS CB")
    #compute DAS CB
      w.sim <- w1.sim + w2.sim
      cband <- matrix(, nrow = nl , ncol = 2)
      qb <- matrix(, nrow = n.sim, ncol = nst)
      ca <- matrix(, nrow = 1, ncol = nst)
      for (s in 1 : ns){
        loc.t12 <- which(indata[[time]] >= tmin.bystrata[s] & indata[[time]] <= tmax.bystrata[s] & indata[[strata]] == bystrata[s])
        qb[, s] <- sort(apply(abs(w.sim[loc.t12, ] / adjse[loc.t12]), 2, max, na.rm = TRUE))
        ca[1, s] <- quantile(qb[, s], probs = c(1 - cb.alpha))
        cband[loc.t12 , 1] = (adjsurv - ca[1, s] * adjse)[loc.t12]
        cband[loc.t12 , 2] = (adjsurv + ca[1, s] * adjse)[loc.t12]
      }
    #output DAS
      adj.s <- as.matrix(cbind(indata[[strata]], indata[[time]], indata[[event]], adjsurv, adjse, cband))
      adj.s <- subset(adj.s, adj.s[, 3] == 1)
      adj.s <- plyr::count(adj.s[order(adj.s[, 1], adj.s[, 2]), ])
      adj.s.out <- cbind(adj.s[, 1 : 2], coxd$nrisk, coxd$nevent, adj.s[, 4 : 7])
      colnames(adj.s.out) <- c('strata', 'time', 'nrisk', 'nevent', 'adjsurv', 'adjSE', 'LCB', 'UCB')
      das <- list()
      das[["DAS"]] <- adj.s.out
toc()

tic("DAS diff+var")
####Difference of adjusted survival w/ variance####
    #Diff between adjusted survival
        adjvdiff <- matrix(, nrow = nt.unique, ncol = nst)
        wdiff <- matrix(, nrow = nt.unique, ncol = nst)
    #DAS, variance for individual strata
        for (s in 1 : ns){
          wdiff[, s] <- adjsurv[byloc[, s]]
          adjvdiff[, s] <- adjvar[byloc[, s]]
        }
    #variance of survival difference between 2 selected groups
        for (s in 1 : nsc){
          s1 <- sc[s, 1]
          s2 <- sc[s, 2]
          #survival difference between 2 selected groups: nt.unique x combination of strata
          wdiff[, ns + s] <- adjsurv[byloc[, s1]]-adjsurv[byloc[, s2]]
          #variance
          adjvdiff[, ns+s] <- suppressWarnings(v1[byloc[, s1]] + v1[byloc[, s2]] +
          diag((ssh[byloc[, s1], ] - ssh[byloc[, s2], ]) %*% cox$var %*% t(ssh[byloc[, s1], ] - ssh[byloc[, s2], ])) / (nl^2))
          }
        adjsediff <- sqrt(adjvdiff)
toc()

####Simulated confidence band####
tic("CB")
    #map to position of pooled time list
    byloc.t12 <- which(bytime >= t1 & bytime <= t2)
    bytime.t12 <- bytime[byloc.t12]
    nt.sim <- NROW(bytime.t12)
    cband.pool <- matrix(, nrow = nt.unique, ncol = (2 * nst))
    pvalue <- matrix(, nrow = 1, ncol = nsc)
    #DAS CB map to pooled time
    for (s in 1 : ns){
      cband.pool[byloc.t12, s * 2 - 1] <- wdiff[byloc.t12, s] - ca[1, s] * adjsediff[byloc.t12, s]
      cband.pool[byloc.t12, s * 2]     <- wdiff[byloc.t12, s] + ca[1, s] * adjsediff[byloc.t12, s]
    }
    #Diff CB
    for (s in 1 : nsc){
      s1 <- sc[s, 1]
      s2 <- sc[s, 2]
      wdiff.sim <- suppressWarnings(w1.sim[byloc[, s1], ] + w1.sim[byloc[, s2], ] - (w2.sim[byloc[, s1], ] - w2.sim[byloc[, s2], ]))
      qb[, ns + s] <- sort(apply(abs(wdiff.sim[byloc.t12, ] / adjsediff[byloc.t12, ns+s]), 2, max, na.rm = TRUE))
      ca[1, ns + s] <- quantile(qb[, ns + s], probs = c(1 - cb.alpha))
      cband.pool[byloc.t12, ((ns + s) * 2 - 1)] <- wdiff[byloc.t12, ns + s] - ca[1, ns + s] * adjsediff[byloc.t12, ns+s]
      cband.pool[byloc.t12, ((ns + s) * 2)]     <- wdiff[byloc.t12, ns + s] + ca[1, ns + s] * adjsediff[byloc.t12, ns+s]
      qmax <- max((abs(wdiff[byloc.t12, ns + s] / adjsediff[byloc.t12, ns + s])), na.rm = TRUE)
      pvalue[s] <- length(which(qb[, ns + s] > qmax)) / n.sim
    }
toc()

####compute N at risk at pooled event time for each strata####
tic("N at risk")
    indata$oneind <- 1
    nrcox <- suppressWarnings(coxph(Surv(indata[[time]],oneind)~strata(indata[[strata]])+factor(indata[[event]]), data = indata, ties = "breslow"))
    nrcoxd <- coxph.detail(nrcox)
    nrdat <- as.matrix(cbind(indata[[strata]], indata[[time]]))
    nrdat <- plyr::count(nrdat[order(nrdat[, 1], nrdat[, 2]), ])
    n.risk <- matrix(, nrow = nt.unique , ncol = ns)
    for (s in 1 : ns){
      for (t in 1 : nt.unique){
        n.risk[t, s] <- suppressWarnings(max(nrcoxd$nrisk[which(nrdat[1] == bystrata[s] & nrdat[2] >= bytime[t])]))
      }
    }
toc()
tic("output dataset")
####output overall pooled dataset####
  das.pool <- NULL
  das.pool.colname <- NULL
  p.names <- NULL
  for (s in 1 : ns){
    das.pool <- data.frame(cbind(das.pool, n.risk[, s], wdiff[, s], adjsediff[, s], cband.pool[, (2 * s - 1) : (2 * s)]))
    das.pool.colname <- c(das.pool.colname, 
                          paste("Nrisk", strata, bystrata[s], sep = "_"),
                          paste("DAS",   strata, bystrata[s], sep = "_"),
                          paste("SE",    strata, bystrata[s], sep = "_"),
                          paste("LCB",   strata, bystrata[s], sep = "_"),
                          paste("UCB",   strata, bystrata[s], sep = "_"))
  }
  for (s in 1 : nsc){
    s1 <- sc[s, 1]
    s2 <- sc[s, 2]
    das.pool <- data.frame(cbind(das.pool, wdiff[, ns + s], adjsediff[, ns + s], cband.pool[, (2 * ns + 2 * s - 1): (2 * ns + 2 * s)]))
    das.pool.colname <- c(das.pool.colname,
                          paste("Diff", strata, bystrata[s1], bystrata[s2], sep = "_"),
                          paste("SE",   strata, bystrata[s1], bystrata[s2], sep = "_"),
                          paste("LCB",  strata, bystrata[s1], bystrata[s2], sep = "_"),
                          paste("UCB",  strata, bystrata[s1], bystrata[s2], sep = "_"))
    p.names <- c(p.names, paste("p",    strata, bystrata[s1], bystrata[s2], sep = "_"))  
  }
  colnames(das.pool) <- das.pool.colname
  das.pool <- data.frame(cbind(bytime, das.pool))
  colnames(das.pool)[1] <- c('Time')
  das[["DASpooled"]] <- das.pool
  #output p-value of confidence band (within selected time range)
  colnames(pvalue) <- p.names
  rownames(pvalue) <- "p-value"
  adj.s.d.pvalue <- data.frame(pvalue)
  das[["diff.p"]] <- adj.s.d.pvalue
  #output simulation parameter
  adj.s.d.cb.sim <- data.frame(cbind(n.sim, t1, t2))
  colnames(adj.s.d.cb.sim) <- c('nsim', 't1', 't2')
  das[["sim.param"]] <- adj.s.d.cb.sim
toc()
toc()
View(cbind(indata[[strata]],indata[[time]],indata[[event]],zbeta,s0beta,h0))
return(das)
}


# das.a.v2<-diff.adj.surv.v2(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+factor(donorgp),seednum=2018,n.sim=2000,starttime=,endtime=)
# das.a.v2<-diff.adj.surv.v2(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(donorgp)+factor(kps)+factor(leuk2),seednum=2018,n.sim=2000,starttime=,endtime=)
# View(das.a.v2$DASpooled)
# View(das.a.v2$DAS)
# rm(das.b.v2)

# das.b.v3<-diff.adj.surv.v3(indata=checked,coxf=Surv(intdxtx,indxsurv,dead)~strata(condted)+factor(yeartx)+factor(karnofcat),
#                            seednum=2018,n.sim=2000,starttime=,endtime=)
das.b.v3<-diff.adj.surv.v3(indata=checkrel,coxf=Surv(intdxtx,indxsurv,dead)~strata(condted)+factor(rel),
                           seednum=2018,n.sim=2000,starttime=,endtime=)
View(das.b.v3$DAS)
View(das.b.v3$DASpooled)


das.a.v2<-diff.adj.surv.v2(indata=m11n10000,coxf=Surv(time,event)~strata(strata)+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,seednum=2018,n.sim=1000,starttime=,endtime=)


cox <- coxph(Surv(indxsurv,dead)~strata(condted)+factor(yeartx)+factor(karnofcat), data = checked, ties = "breslow", model = TRUE)
coxd <- coxph.detail(cox, riskmat = TRUE)

coxlt <- coxph(Surv(intdxtx,indxsurv,dead)~strata(condted)+factor(yeartx)+factor(karnofcat), data = checked, ties = "breslow", model = TRUE)
coxdlt <- coxph.detail(coxlt, riskmat = TRUE)








