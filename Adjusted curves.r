diff.adj.surv <- function(indata,coxf,diff=FALSE,cb=FALSE,
                          seednum=1986,n.sim=2000,cb.alpha=.05,starttime=NULL,endtime=NULL){
####data manipulation####
  #check abnormal arguments
      # check "Surv" object must be there
  if (diff == FALSE & cb == TRUE) stop("DIFF argument must be TRUE for computing confidence band")
  if (cb.alpha <= 0 | cb.alpha >= 1) stop("alpha must be within (0,1)")
  
  #start runtime counting
  ptm<-proc.time()
  #include libraries
  libs<-c("survival","plyr","Matrix")
  lapply(libs, require, character.only = TRUE)
  
  #extract strata variable from formula
  scoxf<-deparse(substitute(coxf[[3]]),width.cutoff = 500)
  strata<-gsub("strata\\((.*?)\\)","\\1",regmatches(scoxf,regexpr("strata\\(.*?\\)",scoxf)))
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
  tm.evt<-terms.y(coxf)
  if (length(tm.evt)!=2) stop("Must specify one time and one event in Surv function")
  time<-tm.evt[1]
  event<-tm.evt[2]

  #exclude case with missing event or time
  if (is.na(indata[[time]]) | is.na(indata[[event]])) warning("Deleted cases w/ missing event or time")
  indata<-subset(indata,!is.na(indata[[time]]) & !is.na(indata[[event]]))
  
  # fit stratified cox model
  cox<-coxph(coxf,data = indata,ties = "breslow",model = TRUE)
  # details of Cox model
  coxd<-coxph.detail(cox,riskmat = TRUE)

####Adjusted survival w/ variance####
      nt<-NROW(coxd$time)             #number of time points (coxph output)
      nz<-NCOL(coxd$x)                #number of risk factor
      nl<-NROW(coxd$x)                #number of cases
      ns<-NROW(coxd$strata)           #number of strata

  #order coxd$x by case (to match with coxd$riskmat and case order from original data)
      coxd$x<-coxd$x[order(as.numeric(row.names(coxd$x))),]
  #obtain trtgp column: nt x 1
      trtlst<-subset(plyr::count(subset(indata,indata[[event]] == 1,select = c(strata,time,event))),select = c(strata))
  #integral matrix: nt x nt (lower triangle of 1 by strata)
      intmat<-as.matrix(Matrix::bdiag(lapply(as.vector(coxd$strata),function(x) lower.tri(matrix(NA,x,x),diag = TRUE)))*1L)
  #event matrix: nt x nl, 1 for a case at event time, 0 otherwise
      evtmat<-matrix(0,nrow = nt,ncol = nl)
      for (t in 1:nt){
        evtmat[t,which(indata[[event]] == 1 & indata[[time]] == coxd$time[t] & indata[[strata]] == trtlst[t,1])] = 1}
  #histmat & riskmat (nl x nl, original case order)
      histmat<-matrix(0,nrow = nl,ncol = nl)
      riskmat<-matrix(0,nrow = nl,ncol = nl)
      for (l in 1:nl){
        histmat[l,which(indata[[strata]] == indata[[strata]][l] & indata[[time]] <= indata[[time]][l])] = 1
        riskmat[l,which(indata[[strata]] == indata[[strata]][l] & indata[[time]] >= indata[[time]][l])] = 1
      }

  #zbeta (by case): nlx1
      zbeta<-coxd$x%*%cox$coefficients
  #weighted zbeta, combination of z and calculated zbeta: nl.w x (nz+zbeta+count)
      zbeta.comb<-plyr::count(cbind(coxd$x,coxd$x%*%cox$coefficients))
  #number of combination of risk factor
      nl.w<-NROW(zbeta.comb)
  #z (weighted): nl.w x nz
      z.w<-zbeta.comb[,1:nz]
  #zbeta (weighted): nl.w x 1
      zbeta.w<-zbeta.comb[,nz+1]
  #weight (number in each z combination): nl.w x 1
      weight<-zbeta.comb[,nz+2]

  #s0beta (by time): ntx1
      s0beta<-t(coxd$riskmat)%*%exp(zbeta)
  #s0beta (by case): nlx1
      s0beta.l<-riskmat%*%exp(zbeta)
  #s1beta: ntxnz
      s1beta<-t(coxd$riskmat)%*%(coxd$x*as.vector(exp(zbeta)))
  #s1beta (by case): nlxnz
      s1beta.l<-riskmat%*%(coxd$x*as.vector(exp(zbeta)))
  #ebeta: ntxnz
      ebeta<-s1beta/as.vector(s0beta)
  #ebeta (by case): nlxnz
      ebeta.l<-s1beta.l / as.vector(s0beta.l)

  #baseline hazard: ntx1
      h0<-coxd$nevent/s0beta
  #baseline cumulative hazard (Breslow's estimator): ntx1
      ch0<-intmat%*%h0
  #Predicted survival probability by Z: ntxnl (unweighted) / ntxnl.w (weighted)
      s.w<-exp(-ch0%*%t(exp(zbeta.w)))
  #Directly adjusted survival prob
      adjsurv<-(s.w%*%weight)/nl

  #variance of adj. survival
    #V1
      v1core<-intmat%*%(coxd$nevent/(s0beta^2))
      v1<-(s.w%*%(weight*exp(zbeta.w)))^2*v1core/(nl^2)
    #V2
      ssh<-matrix(0,nrow = nt,ncol = nz)
      for (i in 1:nl.w){
        z.w.rep<-matrix(unlist(rep(z.w[i,],each = nt)),nrow = nt)
        hmat<-intmat%*%(exp(zbeta.w[i])*(z.w.rep-ebeta)*as.vector(h0))
        ssh<-ssh+weight[i]*s.w[,i]*hmat
        }
      v2<-diag(ssh%*%cox$var%*%t(ssh))/(nl^2)
    #sum
      adjvar<-v1+v2
      adjse<-sqrt(adjvar)

  #output adjsurv w/ var
      adj.s<-data.frame(cbind(trtlst,coxd$time,coxd$nrisk,coxd$nevent,adjsurv,adjse))
      colnames(adj.s)<-c('strata','time','nrisk','nevt','adjsurv','adjse')
      das<-list()
      das[["surv"]]<-adj.s

#####Difference of adjusted survival w/ variance####
  #create pooled time list
      p1<-subset(indata,indata[[event]] == 1,select = c(time,strata,event))
      #create pooled time list with nevent for each strata as a column
      p2<-reshape2::dcast(p1,p1[[time]] ~ p1[[strata]],value.var = event,fun.aggregate = length)
      bytime<-p2[,1]
      nt.unique<-NROW(p2)
      #obtain sorted value of strata variable
      bystrata<-sort(unique(p1[,strata]))
      bystrata.ref<-matrix(unlist(rep(bystrata,each = nt.unique)),nrow = nt.unique)
      #create mapping matrix pointing to the location of adj.surv and other computed matrices
      byloc<-matrix(,nrow = nt.unique,ncol = ns)
      for (s in 1:ns){
        cur.strata<-bystrata[s]
        for (t in 1:nt.unique){
          cur.time<-suppressWarnings(max(adj.s[which(adj.s[,1] == cur.strata & adj.s[,2] <= p2[t,1]),2]))
          byloc[t,s]<-suppressWarnings(max(which(adj.s[,1] == cur.strata & adj.s[,2] == cur.time)))
        }
      }

  if (diff == TRUE){
    #Diff between adjusted survival
        sc<-t(combn(ns,2))
        nsc<-NROW(sc) #sum(1:(ns-1))
        adjvdiff<-matrix(,nrow = nt.unique,ncol = sum(1:ns))
        w<-matrix(,nrow = nt.unique,ncol = nsc)
        w.names<-NULL
        se.names<-NULL
    #variance for individual strata
        for (s in 1:ns){
          adjvdiff[,s]<-adjvar[byloc[,s]]
          se.names<-c(se.names,paste("SE",strata,bystrata[s],sep = "_"))
        }

    #variance of survival difference between 2 selected groups
        for (s in 1:nsc){
          s1<-sc[s,1]
          s2<-sc[s,2]
          #survival difference between 2 selected groups: nt.unique x combination of strata
          w[,s]<-adjsurv[byloc[,s1]]-adjsurv[byloc[,s2]]
          w.names<-c(w.names,paste("Diff",strata,bystrata[s1],bystrata[s2],sep = "_"))
          se.names<-c(se.names,paste("SE",strata,bystrata[s1],bystrata[s2],sep = "_"))
          #variance
          adjvdiff[,ns+s]<-suppressWarnings(v1[byloc[,s1]] + v1[byloc[,s2]] +
          diag((ssh[byloc[,s1],] - ssh[byloc[,s2],]) %*% cox$var %*% t(ssh[byloc[,s1],] - ssh[byloc[,s2],])) / (nl^2)
          )}
        adjsediff<-sqrt(adjvdiff)

  #output difference of adjusted survival
      colnames(w)<-w.names
      colnames(adjsediff)<-se.names
      adj.s.d<-data.frame(cbind(bytime,w,adjsediff))
      colnames(adj.s.d)[1]<-c('time')
      das[["diff"]]<-adj.s.d
  }
#####Simulated confidence band##########
  if (cb == TRUE){
    set.seed(seednum)
    t1<-ifelse(is.null(starttime),min(bytime),starttime)
    t2<-ifelse(is.null(endtime),max(bytime),endtime)
    #map to position of pooled time list
    loc.t12<-which(bytime >= t1 & bytime <= t2)
    bytime.t12<-bytime[loc.t12]
    nt.sim<-NROW(bytime.t12)
    #generate G for each case
    gcase<-matrix(rnorm(as.integer(nl*n.sim), mean = 0, sd = 1), nrow = nl, ncol = n.sim)
    #aggregate G to unique time x strata, match with coxph output
    gevent<-evtmat%*%gcase
    gevent.l<-gcase*as.vector(indata[[event]])
    

    h0.sim<-gevent/as.vector(s0beta)
    ch0.sim<-intmat%*%h0.sim
    w1.sim<-as.vector(-(s.w%*%(weight*exp(zbeta.w))))*ch0.sim/nl
    ubeta.sim<-t(coxd$x-ebeta.l)%*%gevent.l
    
    w.t12<-matrix(,nrow = bytime.t12,ncol = nsc)
    cband<-matrix(,nrow = nt.sim,ncol = (2*nsc))
    pvalue<-matrix(,nrow = 1,ncol = nsc)
    cb.names<-NULL
    p.names<-NULL
    for (s in 1:nsc){
      s1<-sc[s,1]
      s2<-sc[s,2]
      w.sim<-suppressWarnings(w1.sim[byloc[,s1],]+w1.sim[byloc[,s2],]+
            ((ssh[byloc[,s1],]-ssh[byloc[,s2],])%*%cox$var%*%ubeta.sim)/nl)
      qb<-sort(apply(abs(w.sim[loc.t12,]/adjsediff[loc.t12,ns+s]),2,max,na.rm = TRUE))
      c.alpha<-qb[round(n.sim*(1-cb.alpha))]

      cband[,(s*2-1)]<-w[loc.t12,s]-c.alpha*adjsediff[loc.t12,ns+s]
      cband[,(s*2)]  <-w[loc.t12,s]+c.alpha*adjsediff[loc.t12,ns+s]
      qmax<-max((abs(w[loc.t12,s]/adjsediff[loc.t12,ns+s])),na.rm = TRUE)
      pvalue[s]<-length(which(qb>qmax))/n.sim
      cb.names<-c(cb.names,paste("LCB",strata,bystrata[s1],bystrata[s2],sep = "_"),paste("UCB",strata,bystrata[s1],bystrata[s2],sep = "_"))
      p.names<-c(p.names,paste("pvalue",strata,bystrata[s1],bystrata[s2],sep = "_"))
    }

  #output confidence band
    w.t12<-as.matrix(w[loc.t12,])
    colnames(w.t12)<-w.names
    colnames(cband)<-cb.names
    adj.s.d.cb<-data.frame(cbind(bytime.t12,w.t12,cband))
    colnames(adj.s.d.cb)[1]<-c('time')
    das[["diff.cb"]]<-adj.s.d.cb
  #output p-value of confidence band (within selected time range)
    colnames(pvalue)<-p.names
    rownames(pvalue)<-"p-value"
    adj.s.d.pvalue<-data.frame(pvalue)
    das[["diff.p"]]<-adj.s.d.pvalue
  }

####post misc.####
  #end runtime counting
  ptm<-proc.time()-ptm
  print(ptm)


  # return(cbind(coxd$time,h0.sim,ch0.sim))
  # return(as.data.frame(valimat))
  return(das)
  }
#################################################################################################
das.diff.plot <- function(indata,strata,gp1,gp2,ci = TRUE,ci.alpha = 0.05,cb = TRUE,xmax = NULL) {
    library("ggplot2")
    if ((ci.alpha <= 0 | ci.alpha >= 1) & ci == TRUE) stop("alpha must be within (0,1)")
    se.data<-indata$diff
    cb.data<-indata$diff.cb
    p.data<-indata$diff.p
    x.lim<-ifelse(is.null(xmax),NA,xmax)
    
    diff.name<-paste("Diff",strata,gp1,gp2,sep = "_")
    se.data$diff.name<-se.data[,diff.name]
    if (ci == TRUE){
      z<-qnorm(1-ci.alpha/2)
      se.name<-paste("SE",strata,gp1,gp2,sep = "_")
      se.data$se.name<-se.data[,se.name]
      se.data$lci<-se.data$diff.name - z*se.data$se.name
      se.data$uci<-se.data$diff.name + z*se.data$se.name
    }
    if (cb == TRUE){
      lcb.name<-paste("LCB",strata,gp1,gp2,sep = "_")
      ucb.name<-paste("UCB",strata,gp1,gp2,sep = "_")
      p.name<-paste("pvalue",strata,gp1,gp2,sep = "_")
      cb.data$lcb.name<-cb.data[,lcb.name]
      cb.data$ucb.name<-cb.data[,ucb.name]
      p.data$p.name<-p.data[,p.name]
    }
    fig<-ggplot()+
      {if(cb)geom_ribbon(data = cb.data,aes(x = time,ymin = lcb.name,ymax = ucb.name),fill = 'cornflowerblue',alpha = 0.5)}+
      geom_step(data = se.data, mapping = aes(x = time, y = diff.name))+
      {if(ci)geom_step(data = se.data, mapping = aes(x = time, y = lci),linetype = 3)}+
      {if(ci)geom_step(data = se.data, mapping = aes(x = time, y = uci),linetype = 3)}+
      scale_y_continuous(name = "Probability",limits = c(-.5,.5))+
      scale_x_continuous(name = "Time",limits = c(0,x.lim))+
      geom_hline(yintercept = 0,linetype = 5)
  print(p.data$p.name)
  return(fig)
}
#################################################################################################
das.a<-diff.adj.surv(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+factor(donorgp),
                     diff=TRUE,cb=TRUE,seednum=2018628,n.sim=2000,starttime=,endtime=)

das.b<-diff.adj.surv(indata=lk1601,coxf=Surv(intxsurv,dead)~strata(donorgp)+factor(kps)+factor(leuk2),
                     diff=TRUE,cb=TRUE,seednum=2018628,n.sim=2000,starttime=,endtime=)

das.c<-diff.adj.surv(indata=checked,coxf=Surv(intxrel,dfs)~strata(condted)+factor(karnofcat)+factor(sex)+factor(donorgp),
                     diff=TRUE,cb=TRUE,seednum=2018,n.sim=5000,starttime=,endtime=60)

das.d<-diff.adj.surv(indata=checked,coxf=Surv(intxrel,dfs)~strata(condted)+factor(karnofcat)+factor(sex)+age,
                     diff=TRUE,cb=TRUE,seednum=2018628,n.sim=2000,starttime=,endtime=60)

fig.a<-das.diff.plot(indata=das.c,strata="condted",gp1=1,gp2=2,ci=TRUE,ci.alpha=0.05,cb=TRUE,xmax=60)
plot(fig.a)

keny.das=as.data.frame(das.d$surv)
keny.dasdiff=as.data.frame(das.d$diff)
keny.dasdiffcb=as.data.frame(das.c$diff.cb)

write.csv(das.c,"E:/Users/hwang/SAScode/04275 Applied survival analysis/Project/output.csv")



cox<-coxph(Surv(intxsurv,dead)~strata(trtgp)+factor(kps)+factor(leuk2)+factor(donorgp),data=lk1601,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)

cox<-coxph(Surv(intxsurv,dead)~strata(donorgp)+factor(kps)+factor(leuk2),data=lk1601,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)

cb=as.data.frame(a$adj.s.cb)
cb2=as.data.frame(b$adj.s.cb)

cox<-coxph(Surv(intxrel,dfs)~strata(condted)+age,data=checked,ties="breslow",model=TRUE)
coxd<-coxph.detail(cox,riskmat=TRUE)
