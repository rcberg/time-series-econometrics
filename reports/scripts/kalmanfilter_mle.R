library(rfred)
library(tidyverse)

# adapted from jeremy piger's kalman filter coded on matlab

# this function constrains variance parameter (for MLE) to be positive 
cons_fn = function(uncons){
  cons = uncons
  cons[4] = exp(uncons[4])
  cons[5] = exp(uncons[5])
  cons
} 

# first need to get those optimal MLE parameters before we can do any kalman filtering.
# this function takes an initial parameter vector as argument

log_lik = function(theta){
  T = length(y)
  
  theta = cons_fn(theta) # in case the parameters weren't going to be positive...
  
  C = theta[1]
  phi1 = theta[2]
  phi2 = theta[3]
  sig_eta = theta[4]
  sig_v = theta[5]
  
  # setting up the state space form for the filter
  
  alpha = c(C, 0)
  F = matrix( c(phi1, 1, phi2, 0), nrow = 2)
  Q = matrix( c(sig_eta^2, 0, 0, 0), nrow = 2)
  H = matrix(c(1,0),nrow=1)
  R = sig_v^2 
  lambda = 0
  
  # making initial values based off class notes
  
  x_init = solve(diag(sqrt(length(F))) - F )%*%alpha
  v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
  p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
  
  x_lag = x_init 
  p_lag = p_init
  LL = 0 
  
  # the kalman filter (see class notes)
  
  for( t in 1:T ){
    x_pred = alpha + F%*%x_lag
    p_pred = F%*%p_lag%*%t(F) + Q
    
    pred_err = y[t] - lambda - H%*%x_pred
    var_pred_err = H%*%p_pred%*%t(H) + R
    
    x_upd = x_pred + p_pred%*%t(H)%*%solve(var_pred_err)%*%pred_err
    p_upd = p_pred - p_pred%*%t(H)%*%solve(var_pred_err)%*%H%*%p_pred
    
    LL = LL -0.5*log(2*pi) -0.5*log(var_pred_err) - 0.5*(pred_err^2/var_pred_err)
    
    x_lag = x_upd
    p_lag = p_upd
    
  }
  -LL
}

rgdp = r_fred( series = "GDPC1" ,
               freq = "q" ,
               start_date = "1947-01-01" ,
               last_date = "2019-07-01" ) %>%
  mutate(
    level = log(value) , 
    growth = c(NA,diff(level)*100)
  )

y = rgdp$growth[-1] # data
theta_in = c(mean(y) , 0 , 0 , log(sd(y)),1)

solution = optim(theta_in , log_lik ) # maximizes likelihood function for optimal parameters

theta_mle = cons_fn(solution$par) # MLE maximized parameters
exact_ll = solution$value # value of the maximized exact log likelihood

# this will compute the kalman filter. it takes an initial parameter vector as argument.

kfest = function(theta){
  T = length(y)
  
  C = theta[1]
  phi1 = theta[2]
  phi2 = theta[3]
  sig_eta = theta[4]
  sig_v = theta[5]
  
  alpha = c(C , 0)
  F = matrix(c(phi1, 1, phi2, 0), nrow = 2)
  Q = matrix(c(sig_eta^2, 0, 0, 0), nrow = 2)
  H = matrix(c(1,0),nrow=1)
  R = sig_v^2
  lambda = 0
  
  x_init = solve(diag(sqrt(length(F))) - F )%*%alpha
  v_p_init = solve(diag(length( F ))  - kronecker(F,F) )%*%matrix(Q , nrow = length(Q) , ncol = 1)
  p_init = matrix(v_p_init , nrow = sqrt(length(F)) , ncol = sqrt(length(F)) )
  
  x_lag = x_init 
  p_lag = p_init
  
  kf_est = rep(0,T)
  filter_out = data.frame(st_predict = 0, st_update = 0, var_predict = 0, var_update = 0)
  
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
  }
  return(filter_out[-1,])
}

kf_gdp <- 
  kfest(theta_mle) %>%
  bind_cols(growth = y , 
            date = lubridate::ymd(rgdp$date[2:length(rgdp$date)]))%>% # we needed the theta_mle optimized parameters to do this part
  mutate( ci_up = st_update + 1.96*var_predict , 
          ci_low = st_update - 1.96*var_predict )

library(hrbrthemes)

ggplot(data = kf_gdp  ) + 
  geom_line(aes(x = date , y = st_update ), color = 'red' , size = 1) + 
  geom_line(aes(x = date , y = growth),linetype=2) + 
  geom_point(aes(x = date , y = growth)) + 
  labs( x = "Date" , 
        y = "Percent Change (Annualized)" , 
        title = "Kalman filter estimate of real GDP growth trend")+
  theme_ipsum_rc()
