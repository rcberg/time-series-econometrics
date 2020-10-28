library(rfred)
library(tidyverse)
library(lubridate)

# this function constrains variance parameter (for MLE) to be positive 
cons_fn_nobreak = function(uncons){
  cons = uncons
  
  constraint1 = uncons[4] / (1+abs(uncons[4]))
  constraint2 = uncons[5] / (1+abs(uncons[5]))
  
  cons[4] = constraint1 + constraint2
  cons[5] = -(constraint1*constraint2)
  
  cons[2] = exp(uncons[2])
  cons[3] = exp(uncons[3])
  
  return(cons)
} 

cons_fn = function(uncons){
    cons = uncons
    
    constraint1 = uncons[5] / (1+abs(uncons[5]))
    constraint2 = uncons[6] / (1+abs(uncons[6]))
    
    cons[5] = constraint1 + constraint2
    cons[6] = -(constraint1*constraint2)
    
    cons[3] = exp(uncons[3])
    cons[4] = exp(uncons[4])
    
    return(cons)
  } 

# parts adapted from jeremy piger's kalman filter coded on matlab

log_lik_nobreak <- 
  function(theta){
    T = length(y)
    
    theta = cons_fn_nobreak(theta) # in case the parameters weren't going to be positive...
    
    delta = theta[1]
    sig_v = theta[2]
    sig_eta = theta[3]
    phi1 = theta[4]
    phi2 = theta[5]
    
    # setting up the state space form for the filter
    
    alpha = c(0 , 0)
    F = matrix( c(phi1, 1, phi2 , 0), nrow = 2)
    Q = matrix( c(sig_eta^2, 0, 0, 0), nrow = 2)
    H = matrix(c(1, -1), nrow=1)
    R = sig_v^2 
    lambda = delta
    
    # making initial values based off class notes
    
    x_init = solve( diag(sqrt(length(F))) - F )%*%alpha
    v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
    p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
    
    x_lag = x_init 
    p_lag = p_init
    LL = 0 
    
    # the kalman filter (see class notes)
    for( t in 1:T ){
      # these are called the 'prediction equations' of the Kalman filter
      x_pred = alpha + F%*%x_lag
      p_pred = F%*%p_lag%*%t(F) + Q
      
      pred_err = y[t] - lambda - H%*%x_pred
      var_pred_err = H%*%p_pred%*%t(H) + R
      
      # these are called the 'updating equations' of the Kalman filter
      x_upd = x_pred + p_pred%*%t(H)%*%solve(var_pred_err)%*%pred_err
      p_upd = p_pred - p_pred%*%t(H)%*%solve(var_pred_err)%*%H%*%p_pred
      
      LL = LL -0.5*log(2*pi) -0.5*log(var_pred_err) - 0.5*(pred_err^2/var_pred_err)
      
      x_lag = x_upd
      p_lag = p_upd
      
    }
    -LL
  }

log_lik_break <- 
  function(theta, B){
    T = length(y)
    
    theta = cons_fn(theta) # in case the parameters weren't going to be positive...
    
    delta1 = theta[1]
    delta2 = theta[2]
    sig_v = theta[3]
    sig_eta = theta[4]
    phi1 = theta[5]
    phi2 = theta[6]
    
    # setting up the state space form for the filter
    
    alpha = c(0 , 0)
    F = matrix( c(phi1, 1, phi2 , 0), nrow = 2)
    Q = matrix( c(sig_eta^2, 0, 0, 0), nrow = 2)
    H = matrix(c(1, -1), nrow=1)
    R = sig_v^2 
    lambda1 = delta1
    lambda2 = delta2
    
    # making initial values based off class notes
    
    x_init = solve( diag(sqrt(length(F))) - F )%*%alpha
    v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
    p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
    
    x_lag = x_init 
    p_lag = p_init
    LL = 0 
    
    # the kalman filter (see class notes)
    for( t in 1:T ){
      if( t >= B ){
        d_t = 1
      }else{ 
        d_t = 0
      }
      # these are called the 'prediction equations' of the Kalman filter
      x_pred = alpha + F%*%x_lag
      p_pred = F%*%p_lag%*%t(F) + Q
      
      pred_err = y[t] - (1-d_t)*lambda1 - d_t*lambda2 - H%*%x_pred
      var_pred_err = H%*%p_pred%*%t(H) + R
      
      # these are called the 'updating equations' of the Kalman filter
      x_upd = x_pred + p_pred%*%t(H)%*%solve(var_pred_err)%*%pred_err
      p_upd = p_pred - p_pred%*%t(H)%*%solve(var_pred_err)%*%H%*%p_pred
      
      LL = LL -0.5*log(2*pi) -0.5*log(var_pred_err) - 0.5*(pred_err^2/var_pred_err)
      
      x_lag = x_upd
      p_lag = p_upd
      
    }
    -LL
  }

rgdp <- 
  r_fred( series = "GDPC1" ,
          freq = "q" ,
          start_date = "1952-01-01",
          last_date = "2019-07-01" 
  ) %>%
  mutate(
    date = ymd(date) , 
    value = as.numeric(value) , 
    level = log(value) , 
    growth = c(NA,diff(level)*100)
  )

y = rgdp$growth[-1] # data

loglik_df = rep(0,48)
solution_df = list()

for(index in 1:48 ){
  tau = 37+4*index
  
  theta_in = c(mean(y[1:tau]) , mean(y[(tau+1):length(y)]) , -1 , -1 , 1, 1 )
  
  solution = optim(par=theta_in , fn=log_lik_break, B=tau) # maximizes likelihood function for optimal parameters
  loglik_df[index] = solution$value # value of the maximized exact log likelihood
  solution_df[[index]] = cons_fn(solution$par)
} 

# this will compute the kalman filter. it takes an initial parameter vector as argument.

kfest_nobreak = function(theta){
  T = length(y)
  
  theta = cons_fn_nobreak(theta) # in case the parameters weren't going to be positive...
  
  delta = theta[1]
  sig_v = theta[2]
  sig_eta = theta[3]
  phi1 = theta[4]
  phi2 = theta[5]
  
  # setting up the state space form for the filter
  
  alpha = c(0 , 0)
  F = matrix( c(phi1, 1, phi2 , 0), nrow = 2)
  Q = matrix( c(sig_eta^2, 0, 0, 0), nrow = 2)
  H = matrix(c(1, -1), nrow=1)
  R = sig_v^2 
  lambda = delta
  
  # making initial values based off class notes
  
  x_init = solve( diag(sqrt(length(F))) - F )%*%alpha
  v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
  p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
  
  x_lag = x_init 
  p_lag = p_init
  
  filter_out=data.frame(filter_index = 0 ,
                        cycle = 0, 
                        kalman_cyc = 0, 
                        var_predict = 0, 
                        var_update = 0, 
                        k_gain = 0)
  p_pred_list <- list()
  x_pred_list <- list()
  p_upd_list <- list()
  x_upd_list <- list()
  # the kalman filter (see class notes)
  for( t in 1:T ){
    
    x_pred = alpha + F%*%x_lag
    p_pred = F%*%p_lag%*%t(F) + Q
    
    pred_err = y[t] - lambda - H%*%x_pred
    var_pred_err = H%*%p_pred%*%t(H) + R
    kalman_gain = p_pred%*%t(H)%*%solve(var_pred_err)
    
    x_upd = x_pred + p_pred%*%t(H)%*%solve(var_pred_err)%*%pred_err
    p_upd = p_pred - p_pred%*%t(H)%*%solve(var_pred_err)%*%H%*%p_pred
    
    x_lag = x_upd
    p_lag = p_upd
    filter_out = bind_rows(filter_out, 
                           data.frame(st_predict = x_pred[1], 
                                      st_update = x_upd[1], 
                                      var_predict = p_pred[1,1], 
                                      var_update = p_upd[1,1], 
                                      k_gain = kalman_gain[1,1]))
    
    #save model objects for the smoother step
    p_pred_list[[t]] = p_pred
    x_pred_list[[t]] = x_pred
    p_upd_list[[t]] = p_upd
    x_upd_list[[t]] = x_upd
    
  }
  # initialize the smoother
  x_smooth_last = x_upd_list[[T]]
  p_smooth_last = p_upd_list[[T]]
  smooth_out <- data.frame(smooth_index = T, 
                           k_smooth = filter_out$kalman_cyc[T], 
                           var_smooth = filter_out$var_update[T])
  #kalman smoother
  for( t in (T-1):1){
    smooth_var_pred = p_pred_list[[t+1]]
    smooth_var_upd = p_upd_list[[t]]
    smooth_est_pred = x_pred_list[[t+1]]
    smooth_est_upd = x_upd_list[[t]]
    
    c_k = smooth_var_upd%*%t(F)%*%solve(smooth_var_pred)
    x_smooth = smooth_est_upd+c_k%*%(x_smooth_last - smooth_est_pred)
    p_smooth = smooth_var_upd + c_k%*%(p_smooth_last - smooth_var_pred)%*%t(c_k)
    
    smooth_out <- bind_rows(smooth_out,
                            data.frame(smooth_index = t,
                                       k_smooth = x_smooth[1],
                                       var_smooth = p_smooth[1,1])
    )
    x_smooth_last = x_smooth
    p_smooth_last = p_smooth
  }
  data_out <- bind_cols(filter_out[-1,], arrange(smooth_out,smooth_index))
  return(data_out)
}

kfest_break = function(theta, B){
  T = length(y)
  
  theta = cons_fn(theta) # in case the parameters weren't going to be positive...
  
  delta1 = theta[1]
  delta2 = theta[2]
  sig_v = theta[3]
  sig_eta = theta[4]
  phi1 = theta[5]
  phi2 = theta[6]
  
  # setting up the state space form for the filter
  
  alpha = c(0 , 0)
  F = matrix( c(phi1, 1, phi2 , 0), nrow = 2)
  Q = matrix( c(sig_eta^2, 0, 0, 0), nrow = 2)
  H = matrix(c(1, 0), nrow=1)
  R = sig_v^2 
  lambda1 = delta1
  lambda2 = delta2
  
  # making initial values based off class notes
  
  x_init = solve( diag(sqrt(length(F))) - F )%*%alpha
  v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
  p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
  
  x_lag = x_init 
  p_lag = p_init
  
  filter_out=data.frame(filter_index = 0 ,
                        cycle = 0, 
                        kalman_cyc = 0, 
                        var_predict = 0, 
                        var_update = 0, 
                        k_gain = 0)
  p_pred_list <- list()
  x_pred_list <- list()
  p_upd_list <- list()
  x_upd_list <- list()
  # the kalman filter (see class notes)
  for( t in 1:T ){
    if( t >= B ){
      d_t = 1
    }else{ 
      d_t = 0
    }
    x_pred = alpha + F%*%x_lag
    p_pred = F%*%p_lag%*%t(F) + Q
    
    pred_err = y[t] - (1-d_t)*lambda1 - d_t*lambda2 - H%*%x_pred
    var_pred_err = H%*%p_pred%*%t(H) + R
    kalman_gain = p_pred%*%t(H)%*%solve(var_pred_err)
    
    x_upd = x_pred + p_pred%*%t(H)%*%solve(var_pred_err)%*%pred_err
    p_upd = p_pred - p_pred%*%t(H)%*%solve(var_pred_err)%*%H%*%p_pred
    
    x_lag = x_upd
    p_lag = p_upd
    filter_out = bind_rows(filter_out, 
                           data.frame(cycle = x_pred[1], 
                                      kalman_cyc = x_upd[1], 
                                      var_predict = p_pred[1,1], 
                                      var_update = p_upd[1,1], 
                                      k_gain = kalman_gain[1,1],
                                      filter_index = t))
    #save model objects for the smoother step
    p_pred_list[[t]] = p_pred
    x_pred_list[[t]] = x_pred
    p_upd_list[[t]] = p_upd
    x_upd_list[[t]] = x_upd
  }
  # initialize the smoother
  x_smooth_last = x_upd_list[[T]]
  p_smooth_last = p_upd_list[[T]]
  smooth_out <- data.frame(smooth_index = T, 
                           k_smooth = filter_out$kalman_cyc[T], 
                           var_smooth = filter_out$var_update[T])
  #kalman smoother
  for( t in (T-1):1){
    smooth_var_pred = p_pred_list[[t+1]]
    smooth_var_upd = p_upd_list[[t]]
    smooth_est_pred = x_pred_list[[t+1]]
    smooth_est_upd = x_upd_list[[t]]
    
    c_k = smooth_var_upd%*%t(F)%*%solve(smooth_var_pred)
    x_smooth = smooth_est_upd+c_k%*%(x_smooth_last - smooth_est_pred)
    p_smooth = smooth_var_upd + c_k%*%(p_smooth_last - smooth_var_pred)%*%t(c_k)
    
    smooth_out <- bind_rows(smooth_out,
                            data.frame(smooth_index = t,
                                       k_smooth = x_smooth[1],
                                       var_smooth = p_smooth[1,1])
    )
    x_smooth_last = x_smooth
    p_smooth_last = p_smooth
  }
  data_out <- bind_cols(filter_out[-1,], arrange(smooth_out,smooth_index))
  return(data_out)
}

tau_lik <- 37+4*which(loglik_df == min(loglik_df))

theta_mle <- solution_df[[which(loglik_df == min(loglik_df))]]

theta_in_nobreak = c(mean(y) , -1 , -1 , 1, 1 )
theta_mle_nobreak = optim(par=theta_in_nobreak , fn=log_lik_nobreak)$par

kf_uc_nobreak <- 
  kfest_nobreak(theta_mle_nobreak) %>%
  bind_cols( growth = y , 
             date = lubridate::ymd(rgdp$date[2:length(rgdp$date)])) 

kf_optim_break <- 
  kfest_break(theta_mle, B=tau_lik) %>%
  bind_cols( growth = y , 
             date = lubridate::ymd(rgdp$date[2:length(rgdp$date)])) %>%
  mutate( d_t = ifelse( date < rgdp$date[tau_lik], 0, 1) , 
          trend = (1-d_t)*theta_mle[1] + d_t*theta_mle[2],
          kalman_est = kalman_cyc + trend , 
          uc_est = cycle + trend ,
          kalman_smooth = k_smooth+trend)

library(hrbrthemes)

ggplot(data = kf_uc_nobreak  ) + 
  geom_line(aes(x = date , y = st_predict ), color = 'red' , size = 1) + 
  geom_point(aes(x = date , y = growth)) + 
  labs( x = "Date" , 
        y = "Percent Change (Annualized)" , 
        title = "Kalman filter estimate of real GDP growth cycle",
        subtitle = "(No structural break)")+
  theme_ipsum_rc()

ggplot(data = kf_optim_break  ) + 
  geom_line(aes(x = date , y = uc_est ), color = 'red' , size = 1) + 
  geom_point(aes(x = date , y = growth)) + 
  labs( x = "Date" , 
        y = "Percent Change (Annualized)" , 
        title = "Kalman filter estimate of real GDP growth cycle",
        subtitle = "(Structural break at 2005-01-01)")+
  theme_ipsum_rc()
