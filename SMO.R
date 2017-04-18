## SMO
g <- function(x, data_, y, alpha, b, sigma) {
  res = 0
  for(i in 1:nrow(data_)) {
    res = res + y[i]*alpha[i]*K(x, data_[i, ], sigma)
  ##  print(res)  
  }
  res + b
}

K <- function(x1, x2, sigma) {
  exp(-sqrt(sum((x1-x2)^2))/(2*sigma))
}

SMO <- function(data_,
                y_,
                C,
                sigma = 1,
                trace = T) {
  require(foreach)
  
  alpha = rep(0, nrow(data_))
  b = 0
  E = rep(0, nrow(data_))
  
  for(i in 1:1000) {
    if(trace & i%%100==0)
      print(i)
    ## Update E
    invisible(foreach(i=1:nrow(data_)) %dopar% {
      E[i] = g(data_[i, ], data_, y, alpha, b, sigma)-y[i]})
    
    ## search first alpha index
    alpha1_ = -1
    for(i in 1:nrow(data_)) {
      if(alpha[i] > 0 & alpha[i] < C & 
         y[i]*g(data_[i, ], data_, y, alpha, b, sigma)!=1)
        alpha1_ = i
    }
    
    if(alpha1_ == -1) {
      for(i in 1:nrow(data_)) {
        if(alpha[i] == C & y[i]*g(data_[i, ], data_, y, alpha, b, sigma)>1) {
          alpha1_ = i
          break;
        }
        
        if(alpha[i] == 0 & y[i]*g(data_[i, ], data_, y, alpha, b, sigma)<1) {
          alpha1_ = i
          break;
        }
      }
    }
  
    ## search second alpha index
    candidates = which(E==E[-alpha1_][
      if(E[alpha1_] > 0)  which.min(E[-alpha1_]) else which.max(E[-alpha1_])])
    alpha2_ = candidates[which(candidates!=alpha1_)][1]
    
    ## update alpha
    if(y[alpha1_]*y[alpha2_] == -1) {
      L = if(alpha[alpha2_]-alpha[alpha1_] < 0) 0 else  alpha[alpha2_]-alpha[alpha1_]
      H = if(C+alpha[alpha2_]-alpha[alpha1_] > C) C else C+alpha[alpha2_]-alpha[alpha1_]
    } else {
      L = if(alpha[alpha2_]+alpha[alpha1_]-C > 0) alpha[alpha2_]+alpha[alpha1_]-C else 0
      H = if(C < alpha[alpha2_]+alpha[alpha1_]) C else alpha[alpha2_]+alpha[alpha1_]
    }
    
    inta = K(data_[alpha1_, ], data_[alpha1_, ], sigma) + K(data_[alpha2_, ], data_[alpha2_, ], sigma) - 
      2*K(data_[alpha1_, ], data_[alpha2_, ], sigma)
    tmp2 = alpha[alpha2_] + (y[alpha2_]*(E[alpha1_]-E[alpha2_]))/inta
    
    tmp2 = if(tmp2 > H) H else if(tmp2 < L) L else tmp2
    tmp1 = alpha[alpha1_] + y[alpha1_]*y[alpha2_]*(alpha[alpha2_]-tmp2)
    
    if(tmp2 != alpha[alpha2_]) {
      ## update b
      K11 = K(data_[alpha1_, ], data_[alpha1_, ], sigma)
      K12 = K(data_[alpha1_, ], data_[alpha2_, ], sigma)
      K22 = K(data_[alpha2_, ], data_[alpha2_, ], sigma)
      b1 = -E[alpha1_] - y[alpha1_]*K12*(tmp1-alpha[alpha1_]) - y[alpha2_]*K12*(tmp2-alpha[alpha2_]) + b
      b2 = -E[alpha2_] - y[alpha1_]*K12*(tmp1-alpha[alpha1_]) - y[alpha2_]*K22*(tmp2-alpha[alpha2_]) + b
      
      b = (b1 + b2)/2
      alpha[alpha1_] = tmp1
      alpha[alpha2_] = tmp2
    } 
  }
  
  list(alpha = alpha, b = b)
}