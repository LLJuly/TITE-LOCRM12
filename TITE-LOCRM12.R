################################################################################
# TITE-LOCRM12: A Local Continual Reassessment Method for Drug Combination
# Optimization Based on Late-Onset Toxicity and Efficacy Outcomes
#
# Authors:
#   Li Liu
#   Ruitao Lin
#   Nolan A. Wages
################################################################################

# -------------------------------------------------------------------------------------------------------------------
# Inputs of function tite.locrm12
# -------------------------------------------------------------------------------------------------------------------

# p.true.tox        ---> True toxicity probabilities of the dose combinations
# p.true.eff        ---> True efficacy probabilities of the dose combinations
# target.tox        ---> Target DLT probability
# pdc.eff           ---> A dose combination with an efficacy probability greater than or equal to pdc.eff
#                        is considered promising
# ntrial            ---> Number of simulated trials
# nmax              ---> Maximum sample size
# cohort.size       ---> Number of participants in each cohort
# rate              ---> Participant arrival rate, defined as the expected number of cohorts per DLT assessment window
#                        Example: if obswin = 60 and rate = 3, then 3 cohorts are expected to arrive
#                        within 60 time units
# surv              ---> Distribution of time to toxicity ("uniform" or "Weibull")
# surv.eff          ---> Distribution of time to efficacy ("uniform" or "Weibull")
# weight.scheme     ---> Weighting scheme for toxicity ("linear" or "piecewise.linear")
# weight.scheme.eff ---> Weighting scheme for efficacy ("linear", "adaptive", or "interval")
# first.prob        ---> Probability that a DLT occurs during the first half of the toxicity assessment window
#                        when surv = "Weibull" conditional on an event occurring
# first.prob.eff    ---> Probability that an efficacy outcome occurs during the first half of the efficacy
#                        assessment window when surv.eff = "Weibull" conditional on an event occurring
# obswin            ---> Length of the DLT assessment window
# obswin.eff        ---> Length of the efficacy assessment window
# cutoff.tox        ---> Cutoff probability used for the safety stopping rule
# t1                ---> First cut-off time point for the toxicity weight when weight.scheme = "piecewise.linear"
# t2                ---> Second cut-off time point for the toxicity weight when weight.scheme = "piecewise.linear"
# q1                ---> Weight assigned at follow-up time t1 when weight.scheme = "piecewise.linear"
# q2                ---> Weight assigned at follow-up time t2 when weight.scheme = "piecewise.linear"
# random.seed       ---> Random seed for the simulation
# var.a             ---> Variance of the prior distribution for parameter a

# -------------------------------------------------------------------------------------------------------------------

tite.locrm12 <- function(p.true.tox, 
                         p.true.eff, 
                         target.tox=0.35, 
                         pdc.eff=0.45, 
                         ntrial=5000, 
                         nmax=51, 
                         cohort.size=3,
                         rate=3, 
                         surv= "uniform", 
                         surv.eff="uniform",
                         weight.scheme="linear", 
                         weight.scheme.eff="adaptive", 
                         first.prob=0.3, 
                         first.prob.eff=0.3, 
                         obswin=60, 
                         obswin.eff=180,
                         cutoff.tox=0.85, 
                         t1=20, 
                         t2=40,
                         q1=0.1, 
                         q2=0.3, 
                         random.seed=23219, 
                         var.a=1.34){
  set.seed(random.seed)
  
  library(pocrm)
  library(rjags)
  
  if (cohort.size==1){
    stop("Cohort size should be greater than 1.")
  }
  
  if(nmax %% cohort.size != 0){
    stop("Sample size is not divisible by cohort size! Execution stopped.")
  }
  
  skeleton <- getprior(0.05,target.tox,3,5)
  skeleton2 <- getprior(0.05,target.tox,2,2)
  skeleton3 <- getprior(0.05,target.tox,2,3)
  skeleton4 <- getprior(0.05,target.tox,3,4)
  
  crmh = function(a,x,y,w,s) {
    lik = exp(-a^2/2/s^2) 
    for (i in 1:length(x)) {lik = lik * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))}
    return(lik)
  }
  crmht = function(a,x,y,w,s) {
    lik = a * exp(-a^2/2/s^2) 
    for (i in 1:length(x)) {lik = lik * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))}
    return(lik)
  }
  
  ndrugA <- dim(p.true.tox)[2]
  ndrugB <- dim(p.true.tox)[1]
  comb_select_matrix <- matrix(rep(0, length(p.true.tox)), nrow=ndrugA) 
  tox_array <- eff_array <- pts_array <- array(0, c(ndrugA, ndrugB, ntrial))
  
  locationC1 <- NULL
  for(id in 1:ndrugB){locationC1 <- c(locationC1,rep(id,ndrugA))}
  locationC <- cbind(locationC1, rep(1:ndrugA,ndrugB))            
  locationM <- matrix(1:length(p.true.tox), nrow=ndrugB, byrow=T) 
  
  trial.duration <- numeric()
  for(trial in 1:ntrial){ 
    if (trial %% 100 == 0) {cat(paste("trial = ", trial, "\n"))}
    
    comb.curr = 1  
    comb.curr.second <- 1 
    y.true <- numeric() 
    u.true <- numeric() 
    y.true.eff <- numeric() 
    u.true.eff <- numeric() 
    arrival <- c(rep(0, cohort.size)) 
    pat.dose <- numeric() 
    next.arrival <- 0 
    
    for(icohort in 1:(nmax/cohort.size)){
      rand.block <- sample(seq(1,cohort.size), 1) 
      for (i.size in 1:cohort.size) {
        if (i.size == rand.block) {
          dose <- comb.curr.second
        } else {
          dose <- comb.curr
        }
        
        pat.dose <- c(pat.dose, dose)
        p.true <- t(p.true.tox)[dose]
        p.true.e <- t(p.true.eff)[dose]
        
        # toxicity
        if (surv == "uniform") {
          ynew <- rbinom(1,1, p.true)
          unew <- ifelse(ynew, runif(1, 0, obswin), Inf)
        } else if (surv == "Weibull") {
          shape.Weibull <- log(log(1 - p.true)/log(1 - first.prob*p.true))/log(2)
          scale.Weibull <- obswin / ((-log(1 - p.true))^(1/shape.Weibull))
          unew <- rweibull(1, shape.Weibull, scale.Weibull)
          ynew <- as.integer(unew <= obswin)
          if (!ynew) unew <- Inf
        }
        y.true <- c(y.true, ynew)
        u.true <- c(u.true, unew)
        
        # efficacy
        if (surv.eff == "uniform") {
          ynew.eff <- rbinom(1,1, p.true.e)
          unew.eff <- ifelse(ynew.eff, runif(1, 0, obswin.eff), Inf)
        } else if (surv.eff == "Weibull") {
          shape.Weibull.eff <- log(log(1 - p.true.e)/log(1 - first.prob.eff*p.true.e))/log(2)
          scale.Weibull.eff <- obswin.eff / ((-log(1 - p.true.e))^(1/shape.Weibull.eff))
          unew.eff <- rweibull(1, shape.Weibull.eff, scale.Weibull.eff)
          ynew.eff <- as.integer(unew.eff <= obswin.eff)
          if (!ynew.eff) unew.eff <- Inf
        }
        y.true.eff <- c(y.true.eff, ynew.eff)
        u.true.eff <- c(u.true.eff, unew.eff)
      } # Data of a single cohort obtained.
      
      utox <- u.true + arrival 
      ueff <- u.true.eff + arrival 
      
      if (icohort == nmax/cohort.size){
        next.arrival <- next.arrival + obswin.eff
      }
      else if (icohort < nmax/cohort.size) {
        next.arrival <- next.arrival + obswin/rate
      }
      
      B <- rep(0, length(y.true))  
      B[utox<=next.arrival] <- 1 
      B.eff <- rep(0, length(y.true.eff)) 
      B.eff[ueff<=next.arrival] <- 1 
      censor <- pmin(next.arrival, utox) - arrival 
      censor.eff <- pmin(next.arrival, ueff) - arrival 
      followup <- pmin(censor, obswin) 
      followup.eff <- pmin(censor.eff, obswin.eff)
      arrival <- c(arrival, rep(next.arrival, cohort.size))
      
      orders <- matrix(nrow=4, ncol=5)
      d.curr <- locationC[comb.curr,]
      adj <- c(0,0,comb.curr,0,0)
      if(d.curr[1]-1>0){adj[1] <- locationM[d.curr[1]-1,d.curr[2]]} 
      if(d.curr[2]-1>0){adj[2] <- locationM[d.curr[1],d.curr[2]-1]}
      if(d.curr[1]+1<=ndrugB){adj[5] <- locationM[d.curr[1]+1,d.curr[2]]}
      if(d.curr[2]+1<=ndrugA){adj[4] <- locationM[d.curr[1],d.curr[2]+1]}
      orders[1,]<-c(adj[1],adj[2],adj[3],adj[4],adj[5]) 
      orders[2,]<-c(adj[2],adj[1],adj[3],adj[4],adj[5])
      orders[3,]<-c(adj[1],adj[2],adj[3],adj[5],adj[4])
      orders[4,]<-c(adj[2],adj[1],adj[3],adj[5],adj[4])
      # eliminate duplicate orders
      if((d.curr[1]==1 & d.curr[2]==ndrugA) | (d.curr[1]==ndrugB & d.curr[2]==1)){
        orders <- orders[1,]
      }else{
        if(d.curr[1]==1 | d.curr[2]==1){orders <- orders[c(1,3),] }
        if(d.curr[1]==ndrugB | d.curr[2]==ndrugA ){orders <- orders[c(1,2),] }
      }
      # if a single ordering is input as a vector, convert it to a matrix
      if(is.vector(orders)){orders <- t(as.matrix(orders))}
      
      na.index <- which(is.na(t(p.true.tox)))
      orders[orders %in% na.index] <- NA
      orders.re <- matrix(orders[which(orders>0)], nrow=nrow(orders))
      orders.re <- unique(orders.re)
      
      ske.use <- NULL
      if(sum(orders.re[1,]>0)==2){ske.use <- skeleton2} 
      if(sum(orders.re[1,]>0)==3){ske.use <- skeleton3}
      if(sum(orders.re[1,]>0)==4){ske.use <- skeleton4}
      if(sum(orders.re[1,]>0)==5){ske.use <- skeleton}
      ske <- getwm(orders.re, ske.use)
      
      comb.use <- orders.re[1,]
      pat.number.use <- which(pat.dose %in% comb.use)
      pat.dose.use <- pat.dose[pat.number.use] 
      local.rank <- match(pat.dose.use, sort(comb.use)) 
      
      # add suspension rule (wait until all participants' toxicity outside the local space are fully observed)
      pat.number.out <- which(pat.dose %in% setdiff(pat.dose, comb.use))
      tox.out <- B[pat.number.out]
      followup.out <- followup[pat.number.out]
      if (any(tox.out == 0 & followup.out < obswin)){ 
        next.arrival <- next.arrival + max(pmin(u.true[pat.number.out], obswin) - followup.out) 
        B <- rep(0,length(y.true))  
        B[utox<=next.arrival] <- 1
        B.eff <- rep(0, length(y.true.eff)) 
        B.eff[ueff<=next.arrival] <- 1
        censor <- pmin(next.arrival,utox) - arrival[1:(length(arrival)-cohort.size)] 
        censor.eff <- pmin(next.arrival, ueff) - arrival[1:(length(arrival)-cohort.size)]
        followup <- pmin(censor, obswin) 
        followup.eff <- pmin(censor.eff, obswin.eff)
        arrival[(length(arrival)-(cohort.size-1)):length(arrival)] <- rep(next.arrival, cohort.size) 
      }
      
      time.intrial.eff <- pmin(next.arrival - arrival, obswin.eff) 
      
      # check if the lowest does is too toxic
      B.lowest.dose <- B[which(pat.dose == 1)]
      if(pbeta(target.tox,1+sum(B.lowest.dose),1+length(B.lowest.dose)-sum(B.lowest.dose),lower.tail = FALSE)>cutoff.tox){comb.curr <- -1; break}
      
      # use the data in the local space  
      B.use <- B[pat.number.use]
      followup.use <- followup[pat.number.use]
      marginal.vec <- numeric()
      post.tox.prob <- matrix(0, nrow = nrow(ske), ncol = ncol(ske))
      
      if (weight.scheme == "linear"){
        weight = followup.use/obswin
      } else if (weight.scheme == "piecewise.linear"){
        weight = (q1/t1)*pmin(followup.use,t1)+((q2-q1)/(t2-t1))*pmax(0, pmin(followup.use,t2)-t1)+((1-q2)/(obswin-t2))*pmax(0, followup.use-t2)
      }
      
      for (iorder in 1:nrow(orders.re)) { 
        x.use <- NULL 
        for (i in 1:length(local.rank)){
          x.use <- c(x.use, ske[iorder, local.rank[i]])
        }
        marginal=integrate(crmh, lower=-Inf, upper=Inf, x=x.use, y=B.use, w=weight, s=sqrt(var.a), abs.tol = 0)$value;
        est=integrate(crmht, lower=-10, upper=10, x=x.use, y=B.use, w=weight, s=sqrt(var.a), abs.tol = 0)$value/marginal
        post.tox.prob[iorder,] <- ske[iorder, ]**exp(est)
        marginal.vec <- c(marginal.vec, marginal)
      }
      
      mprior.prob = rep(1/nrow(ske), nrow(ske));  
      mpost.prob = (marginal.vec*mprior.prob)/sum((marginal.vec)*mprior.prob)
      mpost.prob.matrix <- t(matrix(rep(mpost.prob, each=ncol(ske)), nrow=ncol(ske)))
      avg.tox <-colSums(mpost.prob.matrix*post.tox.prob)
      
      modelstring.E <- "
      model{
        for (i in 1:length(pat.dose.eff)){
          pe[i] <- pt(alpha + beta1 * doseA[i] + beta2 * doseB[i] + gamma1 * doseA[i]^2 + gamma2 * doseB[i]^2, mean.t, var.t, dof)
          yE[i] ~ dbern(pe[i]*wE[i]) 
        }
          alpha ~ dnorm(mean.alpha, 1/var.alpha)
          beta1 ~ dnorm(mean.beta1, 1/var.beta1)
          beta2 ~ dnorm(mean.beta2, 1/var.beta2)
          gamma1 ~ dnorm(mean.gamma1, 1/var.gamma1)
          gamma2 ~ dnorm(mean.gamma2, 1/var.gamma2)
          dof ~ dunif(dof1,dof2)
      }
      "
      
      adj.nz <- adj[which(adj>0)]
      cand.set9 <- NULL
      for (idose in 1:length(adj.nz)) {
        cand.set9 <- rbind(cand.set9, locationC[adj.nz[idose], ])
      }
      local.eff.dose <- locationM[min(cand.set9[, 1]):max(cand.set9[,1]), min(cand.set9[, 2]):max(cand.set9[,2])] 
      local.eff.dose <- as.vector(local.eff.dose) 
      pat.number.eff <- which(pat.dose %in% local.eff.dose)
      pat.dose.eff <- pat.dose[pat.number.eff]
      
      row.doseA <- c(0.08,0.16,0.24,0.32,0.40)
      doseAunique <- ((row.doseA[1:ndrugA]-mean(row.doseA[1:ndrugA]))/sqrt(var(row.doseA[1:ndrugA]))) 
      doseBuniqie <- ((row.doseA[1:ndrugB]-mean(row.doseA[1:ndrugB]))/sqrt(var(row.doseA[1:ndrugB])))
      
      doseA <- doseAunique[locationC[pat.dose.eff, ][,2]] 
      doseB <- doseBuniqie[locationC[pat.dose.eff, ][,1]] 
      yE.effic <- B.eff[pat.number.eff]
      followup.effic <- followup.eff[pat.number.eff]
      u.true.effic <- u.true.eff[pat.number.eff]
      time.intrial.effic <- time.intrial.eff[pat.number.eff]
      
      # calculate weights for efficacy
      wE.effic <- numeric(length(yE.effic))
      if (weight.scheme.eff == "linear"){
        wE.effic <- followup.effic/obswin.eff
      } else if (weight.scheme.eff == "adaptive"){
        for (i in 1:length(pat.dose.eff)){
          i.dose <- pat.dose.eff[i]
          i.same.pat <- which(pat.dose.eff == i.dose)
          nume <- sum(u.true.effic[i.same.pat] <= followup.effic[i] & time.intrial.effic[i.same.pat] == obswin.eff) + followup.effic[i]/obswin.eff
          denom <- sum(yE.effic[i.same.pat]==1 & time.intrial.effic[i.same.pat]==obswin.eff) + 1
          wE.effic[i] <- nume/denom
        }
      } else if (weight.scheme.eff == "interval"){
        for (i in 1:length(pat.dose.eff)){
          if (followup.effic[i] < obswin.eff / 3) {
            wE.effic[i] <- 0
          } else if (followup.effic[i] < 2 * obswin.eff / 3) {
            wE.effic[i] <- 1/3
          } else if (followup.effic[i] < obswin.eff) {
            wE.effic[i] <- 2/3
          }
        }
      }
      wE.effic[yE.effic == 1 | followup.effic >= obswin.eff] <- 1
      
      jags.data <- list(
        pat.dose.eff = pat.dose.eff,
        doseA = doseA,
        doseB = doseB,
        yE = yE.effic,
        wE = wE.effic,
        mean.alpha=0,   var.alpha =1.3,
        mean.beta1=0.8, var.beta1 =1.3,
        mean.beta2=0.8, var.beta2 =1.3,
        mean.gamma1=0,  var.gamma1 =1.3,
        mean.gamma2=0,  var.gamma2 =1.3,
        mean.t=0,       var.t =1,
        dof1=2,         dof2=10
      )
      
      jags <- jags.model(textConnection(modelstring.E), data =jags.data, n.chains=1, n.adapt=5000, quiet=TRUE)
      e.sample <- coda.samples(jags,c('alpha','beta1','beta2','dof','gamma1','gamma2'), n.iter=2000, progress.bar="none")
      pos <- colMeans(as.matrix(e.sample))
      
      eff <- matrix(nrow = ndrugB, ncol = ndrugA)
      pos <- as.vector(pos)
      for (i in 1:ndrugB){ 
        for (j in 1:ndrugA){
          eff[i,j]= pt(pos[1] + pos[2] * doseAunique[j] + pos[3] * doseBuniqie[i] + pos[5] * doseAunique[j]^2 + pos[6] * doseBuniqie[i]^2, df=pos[4])
        }
      }
      
      # dose assignment       
      loss=abs(avg.tox-target.tox)
      min.loss.loc <- which(loss == min(loss))
      if (length(min.loss.loc) >1) {
        min.loss.loc <- sample(min.loss.loc, 1)
      } 
      MTDC.dlt <- avg.tox[min.loss.loc]
      
      accep.dose <- orders.re[1, which(avg.tox <= MTDC.dlt)]
      
      eff.values <- t(eff)[accep.dose]
      max.eff <- max(eff.values)
      most.eff.dose <- accep.dose[eff.values == max.eff]
      if (length(most.eff.dose) > 1) {
        most.eff.dose <- sample(most.eff.dose, 1)
      }
      
      most.eff.dose.dlt <- avg.tox[which(orders.re[1, ]== most.eff.dose)]
      lower.dose <- orders.re[1, which(avg.tox <= most.eff.dose.dlt)]
      
      if (length(lower.dose)==1) {
        comb.curr.second <- most.eff.dose
      } else if (length(lower.dose) > 1){
        cand.dose <- setdiff(lower.dose, most.eff.dose) 
        cand.dose.dlt <- avg.tox[which(orders.re[1,] %in% cand.dose)]
        max.cand.dose.dlt <- max(cand.dose.dlt)
        comb.curr.second <- orders.re[1, which(avg.tox == max.cand.dose.dlt)]
        comb.curr.second <- setdiff(comb.curr.second, most.eff.dose)
        if (length(comb.curr.second) > 1){
          comb.curr.second <- sample(comb.curr.second,1)
        }
      }
      
      comb.curr <- most.eff.dose
    } # All cohorts are enrolled.
    
    y_matrix <- eff_matrix <- n_matrix <- matrix(rep(0, length(p.true.tox)), nrow = ndrugA) 
    for (i in 1:length(pat.dose)){
      y_matrix[pat.dose[i]] <- y_matrix[pat.dose[i]] + y.true[i]
      eff_matrix[pat.dose[i]] <- eff_matrix[pat.dose[i]] + y.true.eff[i]
      n_matrix[pat.dose[i]] <- n_matrix[pat.dose[i]] + 1
    }
    tox_array[,,trial] <- y_matrix
    eff_array[,,trial] <- eff_matrix
    pts_array[,,trial] <- n_matrix
    trial.duration <- c(trial.duration, next.arrival) 
    
    # Isotonic regression
    phat = (y_matrix + 0.05)/(n_matrix + 0.1)
    phat = Iso::biviso(phat, n_matrix + 0.1, warn = TRUE)[,]
    phat = phat * (n_matrix != 0) + (1e-05) * (matrix(rep(1:dim(n_matrix)[1], each = dim(n_matrix)[2], len = length(n_matrix)), dim(n_matrix)[1],
                                                      byrow = T) + matrix(rep(1:dim(n_matrix)[2], each = dim(n_matrix)[1], len = length(n_matrix)), dim(n_matrix)[1]))
    phat[n_matrix == 0] = 10
    comb.curr2 = which(abs(phat - target.tox) == min(abs(phat - target.tox)))
    if(length(comb.curr2)>1){
      comb.curr2 <- comb.curr2[sample(1:length(comb.curr2),1)]}
    mtd <- locationC[comb.curr2,]
    

    # Efficacy analysis
    # If the trial stops after the first cohort, then the code will not reach modelstring.E, row.doseA, doseAunique, doseBuniqie
    # so we need to define them again here.
    modelstring.E <- "
      model{
        for (i in 1:length(pat.dose.eff)){
          pe[i] <- pt(alpha + beta1 * doseA[i] + beta2 * doseB[i] + gamma1 * doseA[i]^2 + gamma2 * doseB[i]^2, mean.t, var.t, dof)
          yE[i] ~ dbern(pe[i]*wE[i]) 
        }
          alpha ~ dnorm(mean.alpha, 1/var.alpha)
          beta1 ~ dnorm(mean.beta1, 1/var.beta1)
          beta2 ~ dnorm(mean.beta2, 1/var.beta2)
          gamma1 ~ dnorm(mean.gamma1, 1/var.gamma1)
          gamma2 ~ dnorm(mean.gamma2, 1/var.gamma2)
          dof ~ dunif(dof1,dof2)
      }
      "
    row.doseA <- c(0.08,0.16,0.24,0.32,0.40)
    doseAunique <- ((row.doseA[1:ndrugA]-mean(row.doseA[1:ndrugA]))/sqrt(var(row.doseA[1:ndrugA]))) 
    doseBuniqie <- ((row.doseA[1:ndrugB]-mean(row.doseA[1:ndrugB]))/sqrt(var(row.doseA[1:ndrugB])))

    doseA <- doseAunique[locationC[pat.dose, ][,2]] 
    doseB <- doseBuniqie[locationC[pat.dose, ][,1]] 
    yE.effic <- B.eff
    followup.effic <- followup.eff
    wE.effic <- rep(1, length(pat.dose))
    
    jags.data <- list(
      pat.dose.eff = pat.dose,
      doseA = doseA,
      doseB = doseB,
      yE = yE.effic,
      wE = wE.effic,
      mean.alpha=0,   var.alpha =1.3,
      mean.beta1=0.8, var.beta1 =1.3,
      mean.beta2=0.8, var.beta2 =1.3,
      mean.gamma1=0,  var.gamma1 =1.3,
      mean.gamma2=0,  var.gamma2 =1.3,
      mean.t=0,       var.t =1,
      dof1=2,         dof2=10
    ) 
    jags <- jags.model(textConnection(modelstring.E), data =jags.data, n.chains=1, n.adapt=5000, quiet=TRUE)
    e.sample <- coda.samples(jags,c('alpha','beta1','beta2','dof','gamma1','gamma2'), n.iter=2000, progress.bar="none")
    pos <- colMeans(as.matrix(e.sample))
    
    # get the posterior efficacy estimate of dose pairs
    eff <- matrix(nrow = ndrugB, ncol = ndrugA)
    pos <- as.vector(pos)
    for (i in 1:ndrugB){ 
      for (j in 1:ndrugA){
        eff[i,j]= pt(pos[1] + pos[2] * doseAunique[j] + pos[3] * doseBuniqie[i] + pos[5] * doseAunique[j]^2 + pos[6] * doseBuniqie[i]^2, df=pos[4])
      }
    }
    
    for(i in 1:ndrugB){  
      for(j in 1:ndrugA){  
        if (t(phat)[i,j] > t(phat)[mtd[1], mtd[2]]) {
          eff[i, j] <- -1000
        }
      }
    }
    
    comb.curr3 <- which(t(eff)==max(eff))
    
    if (comb.curr > 0){comb_select_matrix[comb.curr3] <- comb_select_matrix[comb.curr3] + 1} # OBDC
    else if (comb.curr == -1) {comb_select_matrix[comb.curr3] <- comb_select_matrix[comb.curr3]}
  } # All trials finished.
  
  # obtain all the metrics based on all simulated trials 
  comb_select_matrix <- t(comb_select_matrix)
  comb_select_matrix[which(is.na(p.true.tox))] <- NA
  comb_select_matrix_pct <- round(100*comb_select_matrix/ntrial, 2)
  
  tox_matrix <- t(round(apply(tox_array, c(1,2), mean), 2))
  tox_matrix[which(is.na(p.true.tox))] <- NA
  ntox <- sum(tox_matrix, na.rm = TRUE) 
  
  eff_matrix <- t(round(apply(eff_array, c(1,2), mean), 2))
  eff_matrix[which(is.na(p.true.tox))] <- NA 
  neff <- sum(eff_matrix, na.rm = TRUE) 
  
  pts_matrix <- t(round(apply(pts_array,c(1,2), mean), 2))
  pts_matrix[which(is.na(p.true.tox))] <- NA
  npts <-  sum(pts_matrix, na.rm = TRUE) 
  
  nstop <- ntrial-sum(comb_select_matrix, na.rm = TRUE)
  stop_percentage <- round(100*nstop/ntrial, 2) 
  
  trial_duration <- round(mean(trial.duration), 2)
  
  #############################################################
  #############    Summary  Evaluation Metrics    #############
  # 1. OBDC metrics
  safe_mask <- p.true.tox <= target.tox
  if (any(safe_mask)) {
    max.eff <- max(p.true.eff[safe_mask])
    obdc_mask <- (p.true.eff == max.eff) & (p.true.tox <= target.tox)
    PCS <- round(sum(comb_select_matrix_pct[obdc_mask]), 2)
    pts.OBDC <- sum(pts_matrix[obdc_mask])
    pct.pts.OBDC <- round(100 * pts.OBDC / nmax, 2)
  } else {
    max.eff <- NA
    PCS <- NA
    pts.OBDC <- NA
    pct.pts.OBDC <- NA
    warning("No safe combinations")
  }
  
  # 2. PDC metrics
  pdc_mask <- (p.true.tox <= target.tox) & (p.true.eff >= pdc.eff)
  if (any(pdc_mask)) {
    PDC.pct <- round(sum(comb_select_matrix_pct[pdc_mask]), 2)
    pts.PDC <- sum(pts_matrix[pdc_mask])
    pct.pts.PDC <- round(100 * pts.PDC / nmax, 2)
  } else {
    PDC.pct <- NA
    pts.PDC <- NA
    pct.pts.PDC <- NA
    warning("No safe promising combinations")
  }
  
  # 3. Overdose metrics
  overdose_mask <- p.true.tox > target.tox
  if (any(overdose_mask)) {
    overdose.pct <- round(sum(comb_select_matrix_pct[overdose_mask]), 2)
    pts.overdose <- sum(pts_matrix[overdose_mask])
    pct.pts.overdose <- round(100 * pts.overdose / nmax, 2)
  } else {
    overdose.pct <- NA
    pts.overdose <- NA
    pct.pts.overdose <- NA
    warning("No overly toxic combinations")
  }
  #############################################################
  
  result.list <- list(pct.sel = comb_select_matrix_pct,
                      tox = tox_matrix, 
                      ntox = ntox, 
                      eff = eff_matrix,
                      neff = neff,
                      pts = pts_matrix,
                      npts = npts,
                      nstop = nstop,
                      stop.pct = stop_percentage,
                      duration = trial_duration,
                      Optimal_Biological_Dose_Combo = "------------------ OBDC ----------------------------",
                      PCS.OBDC = PCS,              
                      pts.OBDC= pts.OBDC,         
                      pct.pts.OBDC = pct.pts.OBDC,
                      Promising_Dose_Combo  = "---------------- PDC -------------------------------",
                      PCS.PDC = PDC.pct,           
                      pts.PDC = pts.PDC,           
                      pct.pts.PDC = pct.pts.PDC,
                      Overdose = "------------------- Overdose ------------------------",
                      PS.overdose = overdose.pct, 
                      pts.overdose = pts.overdose, 
                      pct.pts.overdose = pct.pts.overdose)
  return(result.list)
}




