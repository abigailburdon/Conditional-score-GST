#Function generating entire dataset of patients - produces both longitudinal and
#survival outcomes
generate_ps <- function(n,
                        gamma_v=0.03,
                        eta_f=-0.5,
                        lambda=0.022,
                        sigma_sq=10,
                        h0_val = 5.5,
                        obs_times = c(seq(0, 0.25, length = 6), seq(0.5, 5, by = 0.25)),
                        recruit_end = 2,
                        b_mean = c(6, 3),
                        b_sigma = matrix(c(3.5^2,0,0,2.5^2), 2,2))
  
  
{
  
  
  ###Variables to fix
  b_mat <- rmvnorm(n, b_mean, b_sigma)
  
  #fixed covariates
  Z <- rbern(n, prob = 0.5)
  
  #Generate survival times for each patient
  T_val <- sapply(1:n, function(i){
    
    H <- function(t){
      if(gamma_v!=0){
        v1 <- exp(gamma_v*b_mat[i,1]+eta_f*Z[i])
        v2 <- exp(gamma_v*b_mat[i,2]*t)-1
        v3 <- h0_val*gamma_v*b_mat[i,2]
        H_val <- v1*v2/v3
      }
      else{
        H_val <- t*exp(gamma_v*b_mat[i,1]+eta_f*Z[i])/h0_val
      }
      return(H_val)
    }
    
    #### inversion method e.g T = H^-1(u) ##########
    
    #generate uniform random variable for use in the inversion method
    u <- runif(1,0,1)
    
    if(H(10*obs_times[length(obs_times)]) < -log(1-u)){
      Ts <- Inf
    }else{
      Ts <- uniroot(function(t) H(t) + log(1-u), lower=0,
                    upper = 2*obs_times[length(obs_times)],
                    extendInt = "upX")$root
    }
    
    return(Ts)})
  
  
  
  ####Censoring times ######
  C_val <- rexp(n, lambda)
  status <-  as.numeric(as.numeric(T_val) < as.numeric(C_val))
  T_val[status==0] <- C_val[status==0] #replace event times with censoring times for
  T_val <- as.numeric(T_val)           #censored observations
  
  #Generate arrival times
  arrival_times <- runif(n, 0, recruit_end)
  
  #Order to help with large samples
  gen_ordered <- order(T_val+arrival_times)
  
  #All survival data
  surv.data <- data.frame(patient=c(1:n),
                          treatment=Z[gen_ordered],
                          surv_times = T_val[gen_ordered],
                          status = status[gen_ordered],
                          arrival_times = arrival_times[gen_ordered])
  
  
  
  #Generate longitudinal data for each patient
  long.data <- data.frame(patient = NULL, time_obs = NULL, long_obs = NULL)
  for(i in 1:n){
    
    #generate a random number of longitudinal measurement times
    obs_time <- obs_times[which(obs_times<surv.data$surv_times[i])]
    num_obs <- length(obs_time)
    
    
    #at each measurement time find true value of longitudinal observation, then 
    #add measurement error
    long_obs <- c(b_mat[gen_ordered[i],]%*%t(matrix(c(rep(1, num_obs),obs_time), ncol = 2)))+
      rnorm(num_obs, 0, sqrt(sigma_sq))
    
    long.data <- rbind(long.data, data.frame(patient = rep(surv.data$patient[i], num_obs),
                                             obs_time,
                                             long_obs))}
  
  return(list(surv.data, long.data))}

#Function that takes a dataset of patients at truncates at time time_k. This is
#to be used when truncating data at an interim analysis of a group sequnetial trial
data_cs <- function(surv.data, long.data, time_k){
  
  
  # Survival observations
  surv.data_k <- surv.data[surv.data$arrival_times <= time_k,]
  #Censor unobserved observations and change time
  un_obs_k <- surv.data_k$surv_times+surv.data_k$arrival_times > time_k
  surv.data_k$surv_times[un_obs_k] <- time_k-surv.data_k$arrival_times[un_obs_k]
  surv.data_k$status[un_obs_k] <- 0
  
  #Longitudinal observations
  ld <- sapply(1:nrow(surv.data_k), function(i){
    t_arrive <- surv.data_k$arrival_times[i]
    long_all_i <- long.data[long.data$patient == surv.data_k$patient[i],]
    long_val <- long_all_i[long_all_i$obs_time+t_arrive<=time_k,]
    long_val$patient <- rep(surv.data_k$patient[i], nrow(long_val))
    long_val
  }, simplify = F)
  long.data_k <- do.call(rbind, ld)
  
  return(list(surv.data_k, long.data_k))
}
