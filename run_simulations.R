library(parallel)
library(Rlab)
library(mvtnorm)
library(nlme)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path ))
source("generate_patients.R")
source("conditional_score.R")
source("numerical_integration.R")
source("error_spending.R")

gamma <- 0.03
phi2 <- 6.25
sigma_sq <- 1
eta <- -0.5
delta <- -eta
N <- 50
alpha <- 0.025
beta <- 0.1
b0 <- 6
b1 <- 3
K <- 5
phi1 <- 3.5

cl <- makeCluster(8)  
invisible(
  clusterEvalQ(cl,{
    library(parallel)
    library(Rlab)
    library(mvtnorm)
    library(nlme)
    library(survival)
    source("generate_patients.R")
    source("conditional_score.R")
    source("numerical_integration.R")
    source("error_spending.R")
  }))
clusterExport(cl, c("gamma", "phi2", "sigma_sq", "eta", "delta", "alpha",
                    "beta", "b0", "b1", "K", "phi1"))

###################################
### Calibrate design parameters ###
###################################
{
  ### Baseline hazard
  #Choose value such that median survival on control arm, at mean values of random
  #effects, is 3 years
  h0_val <- uniroot(function(h0){
      #Cumulative hazard function
      if(gamma!=0){
        Hi <- exp(gamma*b0+gamma*b1*3)/(h0*gamma*b1)
      }else{
        Hi <- 3/h0
      }
      #Survival function
      Si <- exp(-Hi)
      #Median survival - Si=0.5
      Si-0.5
    }, lower = 0.1, upper = 100)$root
  
  ### Analysis times
  #Large dataset
  all.data <-  generate_ps(1e3, gamma_v=gamma,
                           eta_f=eta,
                           sigma_sq=sigma_sq,
                           h0_val = h0_val,
                           b_sigma = matrix(c(phi1^2, 0, 0, phi2^2), 2,2))
  surv.data <- all.data[[1]]
  long.data <- all.data[[2]]
  #Analysis times chosen to be that quantiles of event times are equally spaced
  #with final analysis having 60% if events occurred
  ts <- surv.data[surv.data$status==1,"surv_times"]+surv.data[surv.data$status==1,"arrival_times"]
  analysis_times <- as.numeric(quantile(ts, (1:K)*0.6/K))
  
  print(paste("Basline hazard constant is ", h0_val))
  print(paste("Analysis times", analysis_times))
}
  
######################################
### Starting guess for sample size ###
######################################
{
  ### Covariance matrix
  #Truncate data at final analysis
  all.data_k <- data_cs(surv.data, long.data, analysis_times[K])
  surv.data_k <- all.data_k[[1]]
  long.data_k <- all.data_k[[2]]
  #Conditional score items
  cs_lists <- cs_var_indep_sigma(surv.data_k, long.data_k)
  t_vals <- cs_lists[[1]]
  rs <- cs_lists[[2]]
  theta_vals <- cs_lists[[3]]
  X_vals <- cs_lists[[4]]
  #Sandwich variance 
  AB <- sand_var(gamma, eta, sigma_sq,
                       surv.data_k, long.data_k, 1e3,
                       t_vals, rs, theta_vals, X_vals)
  V <- solve(AB[[1]])%*%AB[[2]]%*%solve(AB[[1]])
  
  #Required information and sample size in fixed test
  n_f <- ceiling(V[2,2]*(qnorm(0.975)+qnorm(0.9))^2/(delta^2))
  I_f <- n_f/V[2,2]
  
  #Inflation factor
  R_calc <- error_spending_constants_1s(
    K1 = K,         
    alpha1 = alpha,
    beta1 = beta, 
    delta1 =delta, 
    rho1 = 2)
  #Information and sample size scaled up by inflation factor
  I_max <- I_f*R_calc
  n0 <- ceiling(n_f*R_calc)
  
  print(paste("Maximum information required ", I_max))
  print(paste("Starting guess for sample size", n0))
  
  clusterExport(cl, c("analysis_times","I_max", "h0_val"))
}


########################
### Simulation study ###
########################

### Conditional score method
n_vals <- c()
p_vals <- c()
n <- n0
#We stop when we have 3 estimates for n with:
#(1) All power estimates between 0.85 and 0.95
#(2) monotonic increasing in n against power.
p1 <- p2 <- p3 <- n1 <- n2 <- 1
n3 <- 0
while(p1 >= p2 || p2 >= p3 || n1 >= n2 || n2 >= n3){
  clusterExport(cl, "n")
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
  
    ### Generate data
    all.data <- generate_ps(n, eta_f = eta, gamma_v = gamma, sigma_sq = sigma_sq, h0_val=h0_val,
                            b_sigma = matrix(c(3.5^2,0,0,phi2), 2,2))   
    surv.data <- all.data[[1]]
    long.data <- all.data[[2]]
    
    #Store results
    I_k <- c()
    Z <- c()
    n_k <- c()
    d_k <- c()
    
    #Find treatment effect estimate and information at each analysis of group sequential test
    for(k in 1:K){
      
      ### Data at analysis k
      all.data_k <- data_cs(surv.data, long.data, analysis_times[k])
      surv.data_k <- all.data_k[[1]]
      long.data_k <- all.data_k[[2]]
      n_k[k] <- nrow(surv.data_k)
      
      #Track number of events
      ds <- subset(surv.data_k, status == 1)$patient
      ds_obs <- table(factor(long.data_k$patient, lev = ds))
      d_k[k] <- sum(ds_obs>=2)
      
      ### Conditional score
      #Objects not needed in root finding
      cs_lists <- cs_var_indep_sigma(surv.data_k, long.data_k)
      t_vals <- cs_lists[[1]]
      rs <- cs_lists[[2]]
      theta_vals <- cs_lists[[3]]
      X_vals <- cs_lists[[4]]
      sigma_hat <- sigma_guess(surv.data_k, long.data_k, nrow(surv.data_k),
                               t_vals, rs, theta_vals, X_vals)
      #Conditional score starting point
      cox_model <- data.frame(coefficients=c(gamma,eta))
      #Root of conditional score
      root_find <- optim(
        par = as.numeric(cox_model$coefficients),
        function(x) {
          out <- cs_var_dep_sigma(x[1], x[2], sigma_hat,
                                  surv.data_k, long.data_k,
                                  t_vals, rs, theta_vals, X_vals)
          out[1]^2+out[2]^2})
      
      
      
      ### Sandwich variance
      AB <- sand_var_sigma(root_find$par[1], root_find$par[2], sigma_hat,
                           surv.data_k, long.data_k, n_k[k],
                           t_vals, rs, theta_vals, X_vals)
      A <- AB[[1]]
      B <- AB[[2]]
      A_inv <- try(solve(A), silent=T)
      if((class(A_inv)=="try-error")||sum(is.na(A_inv))>0){
        I_k[k] <- 0
        Z[k] <- 0
      }else{
        #Treatment effect estimate at analysis k
        theta_hat <- root_find$par[2]
        #Information at analysis k
        I_k[k] <- n_k[k]/(A_inv%*%B%*%A_inv)[2,2]
        #Z-statistic at analysis k
        Z[k] <- -theta_hat*sqrt(I_k[k])
      }
      
    }
    
    ### Analysis of trial
    #Boundary points for group sequential test
    bounds <- error_spending_func_1s(
      K=K,      #number of groups
      alpha = alpha, #required significance
      beta = beta,   #required power
      delta = -eta,        #minimum clinically significant difference
      I = I_k,            #Vector of observed information levels
      I_max,
      #Choice of form for the error spending functions - taking input t which is the
      #obtained fraction of the maximum information.
      error_spend_f = function(t, a){
        return(min(a*t^2,a))},
      error_spend_g = function(t, b){
        return(min(b*t^2,b))})
    cross_lower <- (Z < bounds[[1]])
    cross_lower_which <- ifelse(sum(cross_lower)==0, K+1, min(which(cross_lower)))
    cross_upper <- (Z > bounds[[2]])
    cross_upper_which <- ifelse(sum(cross_upper)==0, K+1, min(which(cross_upper)))
    result_gs <- cross_upper_which < cross_lower_which
    which_k <- min(cross_lower_which, cross_upper_which)
    
    #Other results to save
    #Follow-up time
    fut_tot <- apply(surv.data, 1, function(x){
      min(analysis_times[which_k],x[3]+x[5])
    })
    fut <- fut_tot-surv.data$arrival_times
    fut_arrived_gs <- fut[fut>0]
    #Number of hospital visits
    mean_visits_gs <- mean(sapply(surv.data[surv.data$arrival_times<analysis_times[which_k], "patient"], function(i){
      sub.i <- long.data[long.data$patient==i,]
      sum(sub.i$obs_time+surv.data[i, "arrival_times"] <= analysis_times[which_k])
    }))
    
    return(c(Z, I_k, result_gs, which_k,
             bounds[[3]], bounds[[4]],
             mean(fut_arrived_gs), mean_visits_gs))
  
  })
  end_time <- Sys.time()
  print(paste("n=",n,"p=",mean(vals[11,]),"time is ",end_time-start_time))
  
  #Next value for estimate of sample size n - find inflation factor from previous
  #estimates of n
  n_vals <- c(n_vals, n)
  p_vals <- c(p_vals, mean(vals[11,]))
  ratio <- (qnorm(1-alpha)+qnorm(1-beta))/(qnorm(1-alpha)+qnorm(p_vals))
  n <- round(mean(n_vals*ratio^2))
  
  #Update conditions for stoping search over n
  p1 <- p_vals[p_vals==min(p_vals[p_vals > 0.85])]
  n1 <- n_vals[p_vals==min(p_vals[p_vals > 0.85])]
  if(length(p1)==0){
    p1 <- 1
    n1 <- 1e5
  }
  p2 <- p_vals[abs(p_vals-0.9)==min(abs(p_vals-0.9))]
  n2 <- n_vals[abs(p_vals-0.9)==min(abs(p_vals-0.9))]
  p3 <- p_vals[p_vals==max(p_vals[p_vals < 0.95])]
  n3 <- n_vals[p_vals==max(p_vals[p_vals < 0.95])]
}

### Repeat under "naive" Cox model method
n_vals <- c()
p_vals <- c()
n <- n0
#We stop when we have 3 estimates for n with:
#(1) All power estimates between 0.85 and 0.95
#(2) monotonic increasing in n against power.
p1 <- p2 <- p3 <- n1 <- n2 <- 1
n3 <- 0
while(p1 >= p2 || p2 >= p3 || n1 >= n2 || n2 >= n3){
  clusterExport(cl, "n")
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
    
    ### Generate data
    all.data <- generate_ps(n, eta_f = eta, gamma_v = gamma, sigma_sq = sigma_sq, h0_val=h0_val,
                            b_sigma = matrix(c(3.5^2,0,0,phi2), 2,2))   
    surv.data <- all.data[[1]]
    long.data <- all.data[[2]]
    
    #Store results
    I_k <- c()
    Z <- c()
    n_k <- c()
    d_k <- c()
    
    #Find treatment effect estimate and information at each analysis of group sequential test
    for(k in 1:K){
      
      ### Data at analysis k
      all.data_k <- data_cs(surv.data, long.data, analysis_times[k])
      surv.data_k <- all.data_k[[1]]
      long.data_k <- all.data_k[[2]]
      n_k[k] <- nrow(surv.data_k)
      
      #Track number of events
      ds <- subset(surv.data_k, status == 1)$patient
      ds_obs <- table(factor(long.data_k$patient, lev = ds))
      d_k[k] <- sum(ds_obs>=2)
      
      ### Naive estimate
      ### Cox model
      m.cox <- coxph(Surv(surv_times, status)~treatment, data = surv.data_k)
      #Store information and standardised result
      I_k[k] <- 1/summary(m.cox)$coefficients[3]^2
      Z[k] <- -summary(m.cox)$coefficients[4]
      
    }
    
    ### Analysis of trial
    #Boundary points for group sequential test
    bounds <- error_spending_func_1s(
      K=K,      #number of groups
      alpha = alpha, #required significance
      beta = beta,   #required power
      delta = -eta,        #minimum clinically significant difference
      I = I_k,            #Vector of observed information levels
      I_max,
      #Choice of form for the error spending functions - taking input t which is the
      #obtained fraction of the maximum information.
      error_spend_f = function(t, a){
        return(min(a*t^2,a))},
      error_spend_g = function(t, b){
        return(min(b*t^2,b))})
    cross_lower <- (Z < bounds[[1]])
    cross_lower_which <- ifelse(sum(cross_lower)==0, K+1, min(which(cross_lower)))
    cross_upper <- (Z > bounds[[2]])
    cross_upper_which <- ifelse(sum(cross_upper)==0, K+1, min(which(cross_upper)))
    result_gs <- cross_upper_which < cross_lower_which
    which_k <- min(cross_lower_which, cross_upper_which)
    
    #Other results to save
    #Follow-up time
    fut_tot <- apply(surv.data, 1, function(x){
      min(analysis_times[which_k],x[3]+x[5])
    })
    fut <- fut_tot-surv.data$arrival_times
    fut_arrived_gs <- fut[fut>0]
    #Number of hospital visits
    mean_visits_gs <- mean(sapply(surv.data[surv.data$arrival_times<analysis_times[which_k], "patient"], function(i){
      sub.i <- long.data[long.data$patient==i,]
      sum(sub.i$obs_time+surv.data[i, "arrival_times"] <= analysis_times[which_k])
    }))
    
    return(c(Z, I_k, result_gs, which_k,
             bounds[[3]], bounds[[4]],
             mean(fut_arrived_gs), mean_visits_gs))
    
  })
  end_time <- Sys.time()
  print(paste("n=",n,"p=",mean(vals[11,]),"time is ",end_time-start_time))
  
  #Next value for estimate of sample size n - find inflation factor from previous
  #estimates of n
  n_vals <- c(n_vals, n)
  p_vals <- c(p_vals, mean(vals[11,]))
  ratio <- (qnorm(1-alpha)+qnorm(1-beta))/(qnorm(1-alpha)+qnorm(p_vals))
  n <- round(mean(n_vals*ratio^2))
  
  #Update conditions for stoping search over n
  p1 <- (p_vals[p_vals==min(p_vals[p_vals > 0.85])])[1]
  n1 <- (n_vals[p_vals==min(p_vals[p_vals > 0.85])])[1]
  if(is.na(p1)){
    p1 <- 1
    n1 <- 1e5
  }
  p2 <- (p_vals[abs(p_vals-0.9)==min(abs(p_vals-0.9))])[1]
  n2 <- (n_vals[abs(p_vals-0.9)==min(abs(p_vals-0.9))])[1]
  p3 <- (p_vals[p_vals==max(p_vals[p_vals < 0.95])])[1]
  n3 <- (n_vals[p_vals==max(p_vals[p_vals < 0.95])])[1]
}









