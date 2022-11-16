#Choice of form for the error spending functions - taking input t which is the
#obtained fraction of the maximum information.
error_spend_f = function(t, a){
  return(min(a*t^2,a))}
error_spend_g = function(t, b){
  return(min(b*t^2,b))}

error_spending_func_1s <- function(
    K,           #number of groups
    alpha,       #required significance
    beta,        #required power  
    delta,       #power requirement
    I,         #Vector of observed information levels
    I_max,       #maximum information
    
    #Choice of form for the error spending functions - taking input t which is the
    #obtained fraction of the maximum information.
    error_spend_f,
    error_spend_g)

#This function returns the boundary values for an error spending test

{
  
  
  #Pass anaylses that go backwards (or very close) in information
  for(k in 1:(K-1)){
    if(I[k+1] < I[k]){I[k+1] <- I[k]}}
  
  
  #############################
  ### Numerical integration ###
  #############################
  #gridpoints
  r_ind <- c(sapply(2:K, function(k) (I[k]-I[k-1])/I[k]),
             abs((I_max-I)/I_max))
  r_ind[r_ind=="NaN"] <- Inf
  if(min(r_ind)>0.01){r <- 16}
  if(min(r_ind)<=0.01){r <- 64}
  if(min(r_ind)<=0.001){r <- 128}
  
  x_unbounded <- grid(r,0)
  
  
  
  
  #Create empty vectors to store information
  h_null <- list()   #density of continuation region under Null hypothesis
  h_delta <- list()  #density of continuation region under alternative hypothesis
  z <- list()        #list of grid points
  w <- list()        #list of weights
  a_root <- rep(0,K) #lower boundaries
  b_root <- rep(0,K) #upper boundaries
  
  
  
  #Desired amount of error to be spent  
  #Type 1 error
  pi_1k <- rep(0,K)
  pi_1k[1] <- error_spend_f(I[1]/I_max, alpha)
  pi_1k[2:K] <- sapply(2:K, function(k)
  {error_spend_f(I[k]/I_max, alpha)-error_spend_f(I[k-1]/I_max, alpha)})
  #Type 2 error
  pi_2k <- rep(0,K)
  pi_2k[1] <- error_spend_g(I[1]/I_max, beta)
  pi_2k[2:K] <- sapply(2:K, function(k)
  {error_spend_g(I[k]/I_max, beta)-error_spend_g(I[k-1]/I_max, beta)})
  
  
  #Avoid underrunning by spending all error at final anaylsis
  if(I[K] < I_max){
    pi_1k[K] <- alpha-error_spend_f(I[K-1]/I_max, alpha)} 
  
  
  ##################
  ### Analysis 1 ###
  ##################
  b_root[1] <- qnorm(1-pi_1k[1])
  a_root[1] <- delta*sqrt(I[1])+qnorm(pi_2k[1])
  
  
  x <- grid_bounded(x_unbounded, c(a_root[1],b_root[1]))
  z[[1]] <- grid_z(x)
  w[[1]] <- weights(z[[1]])
  h_null[[1]] <- w[[1]]*dnorm(z[[1]])
  h_delta[[1]] <- w[[1]]*dnorm(z[[1]]-delta*sqrt(I[1]))
  
  
  
  
  
  #########################
  ### DEnsity functions ###
  #########################
  
  f_vals <- function(za,zb,I_prev,I_current,theta){
    #Returns the probability of moving from za in analysis k-1 to
    #zb in analysis k.
    dI <- I_current-I_prev
    f_val <- sqrt(I_current)*dnorm((zb*sqrt(I_current)-za*sqrt(I_prev)-theta*dI)/sqrt(dI))/sqrt(dI)
    return(f_val)
  }
  
  
  #Function returning probability of being below "bound"
  spent_error <- function(bound,        #variable
                          I_current,    #Ik
                          I_prev,       #I k-1
                          h_prev,       #densities at previous analysis
                          z_prev,       #gridpoints at previous analysis
                          theta){       #0 for type 1, delta for type 2
    
    #Gridpoints for numerical integration
    x <- grid_bounded(x_unbounded, c(-Inf, bound))
    z_current <- grid_z(x)
    w_current <- weights(z_current)
    
    #Density 
    h_current <- sapply(1:length(z_current), function(i){
      fs <- f_vals(z_prev,z_current[i],I_prev,I_current,theta)
      return(sum(h_prev*fs*w_current[i]))})
    
    #Returns the probability of being below "bound"
    return(sum(h_current))}
  
  
  
  ############################
  ### Analyses k = 2,...,K ###
  ############################
  
  k<-2
  while(I[k-1] < I_max & k <= K){
    #Stop as soon as maximum information achieved
    
    
    #Information going backwards - skip an analysis
    if(I[k] <= I[k-1] | (I[k]-I[k-1])/I_max < 1e-3 ){
      b_root[k] <- b_root[k-1]
      a_root[k] <- a_root[k-1]
      z[[k]] <- z[[k-1]]
      w[[k]] <- w[[k-1]]
      h_null[[k]] <- h_null[[k-1]]
      h_delta[[k]] <- h_delta[[k-1]]
      I[k] <- I[k-1]
      
      
      
    }else{
      #Calculate where the boundary would lie if all error is to be spent at this analysis
      all_error <- try(uniroot(function(u) spent_error(u,I[k],I[k-1],h_null[[k-1]],z[[k-1]],0)+
                                 alpha-error_spend_f(I[k-1]/I_max, alpha)-sum(h_null[[k-1]]), 
                               lower = 0, upper = 3, extendInt="yes")$root, silent =TRUE)
      if(class(all_error)=="try-error"){break}
      
      #Upper boundary - value returning probability of rejection (under Null) as pi_1
      b_root[k] <- uniroot(function(u) spent_error(u,I[k],I[k-1],h_null[[k-1]],z[[k-1]],0)+
                             pi_1k[k]-sum(h_null[[k-1]]), 
                           lower = 0, upper = 3, extendInt="yes")$root
      
      
      
      
      #If the amount of error spent is very close to spending all alpha, stop analysis here
      if(b_root[k]-all_error < 0.01){
        b_root[k] <- all_error
        break
      }
      
      
      if((I_max-I[k])/I_max < 0.015){
        a_root[k] <- b_root[k]
        break
      }
      
      #Lower boundary - value returning probability of accepting (under alternative) as pi_2
      a_root[k] <- uniroot(function(l) spent_error(l,I[k],I[k-1],h_delta[[k-1]],z[[k-1]],delta)-
                             pi_2k[k], 
                           lower = -2, upper = 2, extendInt="yes")$root
      
      
      #If boundaries are upside down, terminate error spending
      if(a_root[k] >= b_root[k]){break}
      
      if(all_error-a_root[k] < 0.01){
        b_root[k] <- all_error
        break
      }
      
      
      #Update variables
      x <- grid_bounded(x_unbounded, c(a_root[k],b_root[k]))
      z[[k]] <- grid_z(x)
      if(length(z[[k]])<=3){
        z[[k]] <- append(z[[k]], z[[k]][1:(length(z[[k]])-1)]+diff(z[[k]])/2)
      }
      w[[k]] <- weights(z[[k]])
      h_null[[k]] <- sapply(1:length(z[[k]]), function(i){
        hs <- h_null[[k-1]] 
        fs <- f_vals(z[[k-1]],z[[k]][i],I[k-1],I[k],0)
        return(sum(hs*fs*w[[k]][i]))})
      h_delta[[k]] <- sapply(1:length(z[[k]]), function(i){
        hs <- h_delta[[k-1]] 
        fs <-  f_vals(z[[k-1]],z[[k]][i],I[k-1],I[k],delta)
        return(sum(hs*fs*w[[k]][i]))})
      
    }
    
    
    k <- k+1
  }
  
  
  #Set boundaries at final analysis equal to each other.
  K <- max(which(b_root!=0))
  a_root[K] <- b_root[K]
  
  
  return(list(a_root, b_root, pi_1k, pi_2k))}


###########################################################################################
one_sided <- function(
    K,             #number of groups
    theta,         #value in Null Hypothesis
    I,             #Observed information levels
    b,             #Upper boundaries
    a )            #lower boundaries
  
  #This function takes as inputs: the number of groups, the required significance
  #level, the power requirements, the boundary values and the observed information
  #levels. The function above should be used prior to find I, b and a. It also takes
  #the value in the Null Hypothesis - put theta= 0 for the significance and change
  #to non-zero for the power. The function returns the error for a one-sided error
  #spending test.
  
{
  
  #value to determine number of grid points. If the change in information is
  #small then the number of grid points should be increased for accuracy.
  if(sum(diff(I)/(I[-1])<0.01)==0){r <- 16}else{r <- 20}   
  
  
  #determine grid points and weights - call functions from "Numerical Integration.R"
  #return list of K vectors, one for each group.
  #points over the whole real line
  x_unbounded <- sapply(1:K, function(k) grid(r,0))
  #remove all points outside the desired range
  x <- sapply(1:K, function(k) grid_bounded(x_unbounded[,k], c(a[k],b[k])))
  if(class(x)=="matrix"){x <- as.list(data.frame(x))}
  #add in the midpoints
  z <- lapply(x, grid_z)
  #find weights for the midpoints - the last analysis should only have one z value
  #as the boundaries join at the final analysis
  w <- lapply(z[1:(K-1)], weights)
  
  
  
  #Probability of staying inside the continuation region per analysis - function h.
  #Analysis k is dependent on analysis k-1 e.g probability of staying within the continuation
  #region at k is conditional on being within the continuation region at k-1.
  
  
  #Analysis 1
  #f1 - density at the first analysis, normal distribution.
  f1 <- function(z){
    dnorm(z- theta*sqrt(I[1]))}
  #h1 values
  h <- list(w[[1]]*f1(z[[1]])) 
  
  #Analysis k=2,...,K
  #fk - density at the analysis k, this depends upon z_k-1
  fk <- function(za,zb,k){
    dIk <- I[k]-I[k-1]
    f_val <- sqrt(I[k])*dnorm((zb*sqrt(I[k])-za*sqrt(I[k-1])-theta*dIk)/sqrt(dIk))/sqrt(dIk)
    return(f_val)
  }
  #h2,...,hK 
  hk <- function(k, i){
    #Function returning h value for given input i at analysis k
    hs <- h[[k-1]]
    w <- w[[k]][i]
    fs <- fk(z[[k-1]],z[[k]][i], k)
    return(sum(hs*w*fs))}
  #loop through values of k from 2 to K-1
  for(k in 2:(K-1)){
    #Apply the function hk accross i=1,...,m_k
    h[[k]] <- sapply(1:length(z[[k]]), function(i) hk(k,i))}
  #Apply the funtion to the last analysis - this only has one point
  hs <- h[[K-1]]
  fs <- fk(z[[K-1]],b[K],K)
  h[[K]] <- sum(hs*fs)
  
  
  
  #Functions for the probability of crossing the boundaries. This only depends on the
  #previous analysis, k-1, and the boundary value ak or bk
  
  #eu - density for greater than upper bound
  eu_k <- function(z, k){
    #Function to return the probability of a point z at analysis k-1 crossing the upper
    #bound at analysis k
    dI_kplus1 <- I[k+1]-I[k]
    val <- (z*sqrt(I[k])+theta*dI_kplus1-b[k+1]*sqrt(I[k+1]))/sqrt(dI_kplus1)
    return(pnorm(val))
  }
  #eu values
  eu <- list()
  for(k in 1:(K-1)){
    eu[[k]] <- sapply(z[[k]], function(z) eu_k(z, k))}
  #el - density for less than lower bound
  el_k <-function(z,k){
    #Function to return the probability of a point z at analysis k-1 crossing the lower
    #bound at analysis k
    dI_kplus1 <- I[k+1]-I[k]
    val <- (-z*sqrt(I[k])-theta*dI_kplus1+a[k+1]*sqrt(I[k+1]))/sqrt(dI_kplus1)
    return(pnorm(val))
  }
  #el values
  el <- list()
  for(k in 1:(K-1)){
    el[[k]] <- sapply(z[[k]], function(z) el_k(z, k))}
  
  
  #final function - probability of exiting at analysis k
  #probability of crossing upper boundary
  psi_k <- function(k){
    #Function returning probability of exiting at analysis k
    sum(h[[k-1]]*eu[[k-1]])
  }
  #probability of crossing lower boundary
  eta_k <- function(k){
    sum(h[[k-1]]*el[[k-1]])
  }
  
  
  #probabilities of exiting the trial
  psik <- c(1-pnorm(b[1]-theta*sqrt(I[1])),sapply(2:K, psi_k))
  etak <- c(pnorm(a[1]-theta*sqrt(I[1])),sapply(2:K, eta_k))
  
  ####overall probability of rejection
  power <- sum(psik)
  
  return(list(power, psik, etak))}



error_spending_constants_1s <- function(
    K1,      #number of groups
    alpha1,  #required significance level
    beta1,   #required power
    delta1,  #minimum clinically significant difference
    rho1)    #rho value for the error spending function
  
  #This function takes the above inputs and returns the sample size ratio for a test
  #with observed information levels that uses an error spending function. The test is 
  #H0:theta=0, for a test H0:theta=theta_0, standardise the test statistic
  
{
  
  
  power_r_variable <- function(r){
    
    #scaled by sample size ratio.
    I_f1 <- (qnorm(1-alpha1)+qnorm(1-beta1))^2/delta1^2 #fixed sample information
    I_max <- I_f1*r    #maximum information
    I1 <- (1:K1)*I_max/K1 #Observed information fractions 
    #scaled by sample size ratio.
    
    
    vals_r <- error_spending_func_1s(
      K = K1,           #number of groups
      alpha = alpha1,  #required significance
      beta = beta1,    #required power
      delta = delta1,
      I = I1,         #Vector of observed information levels
      I_max = I_max,
      #Choice of form for the error spending functions - taking input t which is the
      #obtained fraction of the maximum information.
      error_spend_f = function(t, a){
        return(min(a*t^rho1,a))},
      error_spend_g = function(t, b){
        return(min(b*t^rho1,b))})
    
    #if the loop was broken early, e.g if bk<ak for some analysis then we must define a
    #new maximum number of analyses - the length of the vectors a and b.
    K1 <- sum(vals_r[[1]]!=0)    #maximum number of analyses
    a1 <- vals_r[[1]][1:K1]      #shortened vector for lower boundaries
    b1 <- vals_r[[2]][1:K1]      #shortened vector for upper boundaries
    
    r_error_power <- one_sided(
      K= K1,             #number of groups
      theta = delta1,    #value in Null Hypothesis
      I = I1,            #Observed information levels times sample size ratio
      b = b1,            #Upper boundaries
      a = a1)[[1]]       #lower boundaries
    
    return(r_error_power)}
  
  #The value of R is the value that returns power as we desire - hence use a root finding
  #function to find this value.
  r_root <- uniroot(function(r) power_r_variable(r)-1+beta1,
                    lower = 1, upper = 1.2, extendInt = "yes")$root
  
  
  return(r_root)}