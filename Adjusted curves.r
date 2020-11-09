diff.adj.surv.v2 <- function(indata, coxf, seednum=1986, n.sim=2000, cb.alpha=.05, starttime=NULL, endtime=NULL){

tic("Total time")
####data manipulation####
  #include libraries
  libs <- c("survival", "plyr", "tictoc")
  lapply(libs, require, character.only = TRUE)
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
  if (length(tm.evt) != 2) stop("Must specify one time and one event in Surv function")
  time <- tm.evt[1]
  event <- tm.evt[2]

  #exclude case with missing event or time
  indata <- subset(indata, !is.na(indata[[time]]) & !is.na(indata[[event]]))
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
        riskmat[l, which(indata[[strata]] == indata[[strata]][l] & indata[[time]] >= indata[[time]][l])] = 1
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
return(das)
}


# das.a.v2<-diff.adj.surv.v2(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+factor(donorgp),seednum=2018,n.sim=2000,starttime=,endtime=)
das.a.v2<-diff.adj.surv.v2(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(donorgp)+factor(kps)+factor(leuk2),seednum=2018,n.sim=2000,starttime=,endtime=)
View(das.a.v2$DASpooled)
View(das.a.v2$DAS)

rm(das.b.v2)
das.b.v2<-diff.adj.surv.v2(indata=checked,coxf=Surv(intxrel,dfs)~strata(condted)+factor(yeartx)+factor(karnofcat),
                           seednum=2018,n.sim=2000,starttime=,endtime=5)
View(das.b.v2$DASpooled)
View(das.b.v2$DAS)

das.a.v2<-diff.adj.surv.v2(indata=m11n10000,coxf=Surv(time,event)~strata(strata)+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,seednum=2018,n.sim=1000,starttime=,endtime=)


seq(1, 9, by = 2)

das.diff.plot <- function(indata, strata, gp1, gp2, ci = TRUE, ci.alpha = 0.05, cb = TRUE, xmax = NULL, ymax = .5) {
  libs <- c("ggplot2")
  lapply(libs, require, character.only = TRUE)
    if ((ci.alpha <= 0 | ci.alpha >= 1) & ci == TRUE) stop("alpha must be within (0,1)")
    x.lim <- ifelse(is.null(xmax), NA, xmax)
    y.lim <- ifelse(is.null(ymax), .5, ymax)
    
    se.data <- indata$diff
    diff.name <- paste("Diff", strata, gp1, gp2, sep="_")
    se.data$diff.name <- se.data[, diff.name]
    if (ci == TRUE){
      z <- qnorm(1-ci.alpha/2)
      se.name <- paste("SE", strata, gp1, gp2, sep="_")
      se.data$se.name <- se.data[, se.name]
      se.data$lci <- se.data$diff.name - z*se.data$se.name
      se.data$uci <- se.data$diff.name + z*se.data$se.name
    }
    if (cb == TRUE){
      cb.data <- indata$diff.cb
      p.data <- indata$diff.p
      sim.data <- indata$sim.param
      lcb.name <- paste("LCB", strata, gp1, gp2, sep = "_")
      ucb.name <- paste("UCB", strata, gp1, gp2, sep = "_")
      p.name <- paste("pvalue", strata, gp1, gp2, sep = "_")
      cb.data$lcb.name <- cb.data[, lcb.name]
      cb.data$ucb.name <- cb.data[, ucb.name]
      p.data$p.name <- p.data[, p.name]
      label.p <- paste("p-value = ", p.data$p.name, " for t in [", sim.data$t1, ", ", sim.data$t2, "]", sep = "")
    }
    fig <- ggplot()+
      {if(cb)geom_ribbon(data = cb.data, aes(x = time, ymin = lcb.name, ymax = ucb.name), fill = 'cornflowerblue', alpha = 0.5)}+
      geom_step(data = se.data, mapping = aes(x = time, y = diff.name))+
      {if(ci)geom_step(data = se.data, mapping = aes(x = time, y = lci), linetype = 3)}+
      {if(ci)geom_step(data = se.data, mapping = aes(x = time, y = uci), linetype = 3)}+
      {if(cb)annotate(geom = "text", x = x.lim/6, y = y.lim, label = label.p)}+
      scale_y_continuous(name = "Probability", limits = c(-y.lim, y.lim))+
      scale_x_continuous(name = "Time", limits = c(0, x.lim))+
      geom_hline(yintercept = 0, linetype = 5)
  return(fig)
}


View(cbind(das.b.v2$time,das.b.v2$nevent,das.b.v2$nrisk))




fig.a<-das.diff.plot(indata=das.a.v2,strata="trtgp",gp1=0,gp2=1,ci=TRUE,ci.alpha=0.05,cb=TRUE,xmax=60)
plot(fig.a)

View(cbind(das.a.v2$diff.cb,das.a$diff.cb))

View(as.data.frame(das.a.v2))

write.csv(das.c,"E:/Users/hwang/SAScode/04275 Applied survival analysis/Project/output.csv")

cox<-coxph(Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+factor(donorgp),data=lk1601,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)

cox<-coxph(Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+age,data=lk1601,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)

cox<-coxph(Surv(intxrel,dfs)~strata(condted)+factor(sex),data=checked,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)
View(coxd$x)



coxd$x <- coxd$x[order(as.numeric(row.names(coxd$x))), ]

View(as.data.frame(coxd$x))

cb=as.data.frame(a$adj.s.cb)
cb2=as.data.frame(b$adj.s.cb)


