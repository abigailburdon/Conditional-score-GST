cs_var_indep_sigma <- function(surv.data, long.data){
  #Some necessary parts of the conditional score function do not involve gamma, eta
  #and sigma - for root finding, these functions only need to be evaluated once, not
  #for varying values of gamma, eta, sigma.
  n <- nrow(surv.data) #number of patients
  
  
  
  #Find all event times which will contribute to the conditional score sum
  t_vals <- sapply(1:n, function(i){
    #Event times only contribute if the observation is exact and there are enough
    #longitudinal measurements before t to do regression
    no_censoring <- surv.data$status[i]==1
    all_long_i <- long.data[long.data$patient==i,]
    enough_measures <- nrow(all_long_i)>=2
    return(no_censoring&enough_measures)})
  
  t_vals <- data.frame(patient=which(t_vals),
                       time = surv.data[which(t_vals),"surv_times"])
  
  
  
  
  #A function returning all patients at risk at a particular time
  #A patient is at risk at time t if they are still alive and also if they have
  #enough longitudinal measurements before t to do regression
  risk_set <- function(t, surv.data, long.data){
    ds <- surv.data[surv.data$surv_times>=t,"patient"]
    ds_obs <- table(factor(long.data$patient, lev = ds))
    return(ds[ds_obs>=2])
  }
  rs <- sapply(1:nrow(t_vals), function(j){
    risk_set(t_vals[j, "time"], surv.data, long.data)
  })
  
  
  #variance of the predicted observation
  theta <- function(i, t, surv.data, long.data){
    #function takes as inputs: patient number and time,
    #returns variance of the predicted longitudinal measurement
    ss <- long.data[long.data$patient == i, ]
    ts <- ss[ss$obs_time <=t , "obs_time"]
    1/length(ts) + ((t-mean(ts))^2) / sum((ts-mean(ts))^2)
  }
  theta_Vec <- Vectorize(theta, vectorize.args = "i")
  theta_vals <- sapply(1:nrow(t_vals), function(j){
    theta_Vec(rs[[j]], t_vals[j, "time"], surv.data, long.data)
  })
  
  
  X_hat <- function(i, t, long.data){
    #function takes as inputs: patient number and time,
    #returns value of predicted longitudinal  measurement
    set1 <- long.data[long.data$patient==i,]
    set2 <- set1[set1$obs_time<t,]
    x_obs <- set2$obs_time
    y_obs <- set2$long_obs
    
    x_bar <- mean(x_obs)
    y_bar <- mean(y_obs)
    
    #ordinary least squares
    b1_hat <- sum((x_obs-x_bar)*(y_obs-y_bar))/sum(((x_obs-x_bar))^2)
    b0_hat <- y_bar-b1_hat*x_bar
    
    return(b0_hat+b1_hat*t)}
  
  
  X_vals <- list(NULL)
  for(j in 1:nrow(t_vals)){
    X_vals[[j]] <- sapply(rs[[j]], function(i) X_hat(i, t_vals[j, "time"], long.data))}
  
  
  
  #return all values that will be used again in the part dependent on gamma, eta
  return(list(t_vals, rs, theta_vals, X_vals))
}

sigma_guess <- function(surv.data, long.data, n,
                        t_vals, rs, theta_vals, X_vals){
  
  #Residual sum of squares needed for estimating sigma_sq
  RSS <- function(i, long.data){
    #Function resturning the resisdual sum of squares for patient i
    set <- subset(long.data, patient==i)
    x_obs <- set$obs_time  #some cases result in error when all observations are zero
    
    if(nrow(set)>2 & sum(x_obs)!=0){
      #Patient i only contributes if they have more than two longitudinal measurement
      #this is for regression
      y_obs <- set$long_obs
      
      x_bar <- mean(x_obs)
      y_bar <- mean(y_obs)
      
      b1_hat <- sum((x_obs-x_bar)*(y_obs-y_bar))/sum(((x_obs-x_bar))^2)
      b0_hat <- y_bar-b1_hat*x_bar
      
      num <- sum((y_obs-b0_hat-b1_hat*x_obs)^2)
      denom <- length(y_obs)-2
    }else{
      num <- 0
      denom <- 0
    }
    return(c(num, denom))
  }
  rss <- rowSums(sapply(1:n, function(i) RSS(i, long.data)))
  
  return(rss[1]/rss[2])
  
}

cs_var_dep_sigma <- function(gamma, eta, sigma_sq,
                             surv.data, long.data,
                             t_vals, rs, theta_vals, X_vals){
  colSums(cs_var_dep_rows_sigma(gamma, eta, sigma_sq,
                                surv.data, long.data,
                                t_vals, rs, theta_vals, X_vals))
}

cs_var_dep_rows_sigma <- function(gamma, eta, sigma_sq,
                                  surv.data, long.data,
                                  t_vals, rs, theta_vals, X_vals){
  #This part of the function includes gamma and eta so this part will be evaluated
  #multiple times for varying gamma, eta.
  
  
  #Define empty lists
  S <- list(NULL)
  S_i <- rep(0, nrow(t_vals))
  E0 <- list(NULL)
  Z <- list(NULL)
  RSS <- rep(0,nrow(t_vals))
  
  for(j in 1:nrow(t_vals)){
    
    #keep a marker for the patient who's event time we're interested in
    #where they lie in the risk set
    patient_marker <- which(rs[[j]]== t_vals[j,"patient"])
    
    
    
    #sufficient statistic
    S[[j]] <- X_vals[[j]]
    S_i[j] <- S[[j]][patient_marker]+gamma*sigma_sq*theta_vals[[j]][patient_marker]
    S[[j]][patient_marker] <- S_i[j]
    
    
    #covariates
    Z[[j]] <- surv.data[surv.data$patient%in%rs[[j]],"treatment"]
    
    
    
    #condtitional density
    E0[[j]] <- exp(gamma*S[[j]]-
                     gamma^2*sigma_sq*theta_vals[[j]]/2 +
                     eta*Z[[j]])
    
    #RSS for estimation of sigma_sq
    set <- subset(long.data, patient==t_vals[j,"patient"])
    RSS[j] <- ifelse(nrow(set)<=2, 0, {
      Wi <- matrix(set$long_obs, ncol = 1)
      Ai <- matrix(c(rep(1, nrow(set)), set$obs_time), ncol = 2)
      alphai <- solve(t(Ai)%*%Ai)%*%t(Ai)%*%Wi
      
      t(Wi-Ai%*%alphai)%*%(Wi-Ai%*%alphai)-sigma_sq*nrow(set)+2*sigma_sq
    })
  }
  
  #Conditional score function - create vector where each element is the contribution to
  #the conditional score function at an event time
  
  #Vector E1 at each event time
  E1_1 <- mapply(function(x,y) x*y, S, E0)
  E1_2 <- mapply(function(x,y) x*y, Z, E0)
  
  #contribution at event time T_i
  vals <- matrix(0, nrow = nrow(t_vals), ncol = 3)
  for(i in 1:nrow(t_vals)){
    vals[i,1] <- S_i[i]-sum(E1_1[[i]])/sum(E0[[i]])
    vals[i,2] <- surv.data[surv.data$patient==t_vals[i,"patient"],"treatment"]-sum(E1_2[[i]])/sum(E0[[i]])
    vals[i,3] <- RSS[i]
  }
  
  
  #Return the sum over all event times
  return(vals)
}

sand_var_sigma <- function(gamma, eta, sigma_sq,
                           surv.data, long.data, n,
                           t_vals, rs, theta_vals, X_vals){
  
  
  
  #Define empty lists
  #denominator
  E0 <- rep(0, nrow(t_vals))
  #numerator 1
  E11 <- rep(0, nrow(t_vals))
  #numerator 2
  E12 <- rep(0, nrow(t_vals))
  #dE0/d gamma
  dE01 <- rep(0, nrow(t_vals))
  #dE0/d eta
  dE02 <- rep(0, nrow(t_vals))
  #dE0/d sigma_sq
  dE03 <- rep(0, nrow(t_vals))
  #dE11/d gamma
  dE111 <- rep(0, nrow(t_vals))
  #dE11/d eta
  dE112 <- rep(0, nrow(t_vals))
  #dE11/d sigma_sq
  dE113 <- rep(0, nrow(t_vals))
  #dE12/d gamma
  dE121 <- rep(0, nrow(t_vals))
  #dE12/d eta
  dE122 <- rep(0, nrow(t_vals))
  #dE12/d sigma_sq
  dE123 <- rep(0, nrow(t_vals))
  #dS/d gamma _i
  dS_gamma_i <- rep(0, nrow(t_vals))
  #dS/d sigma_sq _i
  dS_sigma_sq_i <- rep(0, nrow(t_vals))
  
  #mi counts for sigma_sq
  mi <- rep(0, nrow(t_vals))
  
  
  for(j in 1:nrow(t_vals)){
    
    #keep track of which patient in the risk set is the one who's event time we're interested in.
    patient_marker <- which(rs[[j]]== t_vals[j,"patient"])
    
    #number of longitudinal observations
    set <- subset(long.data, patient == t_vals[j,"patient"])
    mi[j] <- ifelse(nrow(set)>2, 2-nrow(set), 0)
    
    #sufficient statistic
    S <- X_vals[[j]]
    S[patient_marker] <- S[patient_marker]+gamma*sigma_sq*theta_vals[[j]][patient_marker]
    
    #derivatives
    dS_gamma_i[j] <- sigma_sq*theta_vals[[j]][patient_marker]
    dS_gamma <- rep(0, length(S))
    dS_gamma[patient_marker] <- dS_gamma_i[j]
    
    dS_sigma_sq_i[j] <- gamma*theta_vals[[j]][patient_marker]
    dS_sigma_sq <- rep(0, length(S))
    dS_sigma_sq[patient_marker] <- dS_sigma_sq_i[j]
    
    
    
    #covariates
    Z <- surv.data[surv.data$patient%in%rs[[j]],"treatment"]
    
    
    #condtitional density
    E0i <- exp(gamma*S-
                 gamma^2*sigma_sq*theta_vals[[j]]/2 +
                 eta*Z)
    
    
    E0[j] <- sum(E0i)
    E11[j] <- sum(S*E0i)
    E12[j] <- sum(Z*E0i)
    
    
    dE0j_gamma <- (S+gamma*dS_gamma-gamma*sigma_sq*theta_vals[[j]])*E0i
    dE0j_sigma_sq <- (gamma*dS_sigma_sq-gamma^2*theta_vals[[j]]/2)*E0i
    
    
    
    #derivatives of conditional density
    #dE0/d gamma
    dE01[j] <- sum(dE0j_gamma)
    #dE0/d eta
    dE02[j] <- sum(Z*E0i)
    #dE0/d sigma_sq
    dE03[j] <- sum(dE0j_sigma_sq)
    #dE11/d gamma
    dE111[j] <- sum(S*dE0j_gamma+dS_gamma*E0i)
    #dE11/d eta
    dE112[j] <- sum(S*Z*E0i)
    #dE11/d sigma_sq
    dE113[j] <- sum(S*dE0j_sigma_sq+dS_sigma_sq*E0i)
    #dE12/d gamma
    dE121[j] <- sum(Z*dE0j_gamma)
    #dE12/d eta
    dE122[j] <- sum(Z^2*E0i)
    #dE12/d sigma_sq
    dE123[j] <- sum(Z*dE0j_sigma_sq)
  }
  
  
  
  
  #Derivatives of score statistic
  #dU1/d gamma
  dU11 <- sum(dS_gamma_i + (E11/E0)*(dE01/E0)-dE111/E0)
  #dU1/d eta
  dU12 <- sum((E11/E0)*(dE02/E0)-dE112/E0)
  #dU1/d sigma_sq
  dU13 <- sum(dS_sigma_sq_i + (E11/E0)*(dE03/E0)-dE113/E0)
  #dU2/d gamma
  dU21 <- sum((E12/E0)*(dE01/E0)-dE121/E0)
  #dU2/d eta
  dU22 <- sum((E12/E0)*(dE02/E0)-dE122/E0)
  #dU2/d sigma_sq
  dU23 <- sum((E12/E0)*(dE03/E0)-dE123/E0)
  #dU3/d gamma
  dU31 <- 0
  #dU3/d eta
  dU32 <- 0
  #dU3/d sigma_sq
  dU33 <- sum(mi)
  
  #A matrix
  A <- matrix(c(dU11, dU21, dU31,
                dU12, dU22, dU32,
                dU13, dU23, dU33), nrow=3)/n
  
  
  vals <- matrix(0, n, ncol=3)
  vals[1:nrow(t_vals),] <- cs_var_dep_rows_sigma(gamma, eta, sigma_sq,
                                                 surv.data, long.data,
                                                 t_vals, rs, theta_vals, X_vals)
  B <- cov(vals)
  
  
  
  return(list(A,B))}

sand_var <- function(gamma, eta, sigma_sq,
                           surv.data, long.data, n,
                           t_vals, rs, theta_vals, X_vals){
  
  
  
  #Define empty lists
  #denominator
  E0 <- rep(0, nrow(t_vals))
  #numerator 1
  E11 <- rep(0, nrow(t_vals))
  #numerator 2
  E12 <- rep(0, nrow(t_vals))
  #dE0/d gamma
  dE01 <- rep(0, nrow(t_vals))
  #dE0/d eta
  dE02 <- rep(0, nrow(t_vals))
  #dE0/d sigma_sq
  dE03 <- rep(0, nrow(t_vals))
  #dE11/d gamma
  dE111 <- rep(0, nrow(t_vals))
  #dE11/d eta
  dE112 <- rep(0, nrow(t_vals))
  #dE11/d sigma_sq
  dE113 <- rep(0, nrow(t_vals))
  #dE12/d gamma
  dE121 <- rep(0, nrow(t_vals))
  #dE12/d eta
  dE122 <- rep(0, nrow(t_vals))
  #dE12/d sigma_sq
  dE123 <- rep(0, nrow(t_vals))
  #dS/d gamma _i
  dS_gamma_i <- rep(0, nrow(t_vals))
  #dS/d sigma_sq _i
  dS_sigma_sq_i <- rep(0, nrow(t_vals))
  
  #mi counts for sigma_sq
  mi <- rep(0, nrow(t_vals))
  
  
  for(j in 1:nrow(t_vals)){
    
    #keep track of which patient in the risk set is the one who's event time we're interested in.
    patient_marker <- which(rs[[j]]== t_vals[j,"patient"])
    
    #number of longitudinal observations
    set <- subset(long.data, patient == t_vals[j,"patient"])
    mi[j] <- ifelse(nrow(set)>2, 2-nrow(set), 0)
    
    #sufficient statistic
    S <- X_vals[[j]]
    S[patient_marker] <- S[patient_marker]+gamma*sigma_sq*theta_vals[[j]][patient_marker]
    
    #derivatives
    dS_gamma_i[j] <- sigma_sq*theta_vals[[j]][patient_marker]
    dS_gamma <- rep(0, length(S))
    dS_gamma[patient_marker] <- dS_gamma_i[j]
    
    dS_sigma_sq_i[j] <- gamma*theta_vals[[j]][patient_marker]
    dS_sigma_sq <- rep(0, length(S))
    dS_sigma_sq[patient_marker] <- dS_sigma_sq_i[j]
    
    
    
    #covariates
    Z <- surv.data[surv.data$patient%in%rs[[j]],"treatment"]
    
    
    #condtitional density
    E0i <- exp(gamma*S-
                 gamma^2*sigma_sq*theta_vals[[j]]/2 +
                 eta*Z)
    
    
    E0[j] <- sum(E0i)
    E11[j] <- sum(S*E0i)
    E12[j] <- sum(Z*E0i)
    
    
    dE0j_gamma <- (S+gamma*dS_gamma-gamma*sigma_sq*theta_vals[[j]])*E0i
    dE0j_sigma_sq <- (gamma*dS_sigma_sq-gamma^2*theta_vals[[j]]/2)*E0i
    
    
    
    #derivatives of conditional density
    #dE0/d gamma
    dE01[j] <- sum(dE0j_gamma)
    #dE0/d eta
    dE02[j] <- sum(Z*E0i)
    #dE0/d sigma_sq
    dE03[j] <- sum(dE0j_sigma_sq)
    #dE11/d gamma
    dE111[j] <- sum(S*dE0j_gamma+dS_gamma*E0i)
    #dE11/d eta
    dE112[j] <- sum(S*Z*E0i)
    #dE11/d sigma_sq
    dE113[j] <- sum(S*dE0j_sigma_sq+dS_sigma_sq*E0i)
    #dE12/d gamma
    dE121[j] <- sum(Z*dE0j_gamma)
    #dE12/d eta
    dE122[j] <- sum(Z^2*E0i)
    #dE12/d sigma_sq
    dE123[j] <- sum(Z*dE0j_sigma_sq)
  }
  
  
  
  
  #Derivatives of score statistic
  #dU1/d gamma
  dU11 <- sum(dS_gamma_i + (E11/E0)*(dE01/E0)-dE111/E0)
  #dU1/d eta
  dU12 <- sum((E11/E0)*(dE02/E0)-dE112/E0)
  #dU1/d sigma_sq
  dU13 <- sum(dS_sigma_sq_i + (E11/E0)*(dE03/E0)-dE113/E0)
  #dU2/d gamma
  dU21 <- sum((E12/E0)*(dE01/E0)-dE121/E0)
  #dU2/d eta
  dU22 <- sum((E12/E0)*(dE02/E0)-dE122/E0)
  #dU2/d sigma_sq
  dU23 <- sum((E12/E0)*(dE03/E0)-dE123/E0)
  #dU3/d gamma
  dU31 <- 0
  #dU3/d eta
  dU32 <- 0
  #dU3/d sigma_sq
  dU33 <- sum(mi)
  
  #A matrix
  A <- matrix(c(dU11, dU21, dU31,
                dU12, dU22, dU32,
                dU13, dU23, dU33), nrow=3)/n
  
  
  vals <- matrix(0, n, ncol=3)
  vals[1:nrow(t_vals),] <- cs_var_dep_rows_sigma(gamma, eta, sigma_sq,
                                                 surv.data, long.data,
                                                 t_vals, rs, theta_vals, X_vals)
  B <- cov(vals)
  
  
  
  return(list(A,B))}