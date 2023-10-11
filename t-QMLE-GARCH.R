#==================================================
# 这是一个可以用于估计各种模型的参数的函数。
# 这种方法的关键思想是根据数据本身估计参数,并假设误差的分布
# 定义GARCH(q,p) q冲击 p波动持久
#==================================================

# Adaptive Quasi Maximum Likelihood Estimation
A_tQMLE <- function(series, order = c(1,1)){
  # param series The original time-series that need to be fitted as GARCH model
  # param order The GARCH model orders. Includes GARCH terms and ARCH terms
  q <- order[2]; p <- order[1] # p意味着冲击性，q意味着波动持久性
  n <- length(series); max.iter <- 50
  Bdf.t <- c()
  
  # Estm包含参数est.df、ETA和SSE（总和平方误差）
  Estm <- matrix(NA, max.iter, (p+q+1+3))
  Est.model <- QMLE(series=series, LogLFunc = "LogL_GARCH_Norm", order = order) # QMLE function
  para <- Est.model$QMLE.N
  e.t <- Est.model$e
  Estm <- c(Est.model$QMLE.N, NA, NA, NA)
  names(Estm) <- c("alpha0", "alpha1", "beta","est.df", "eta", "diff_Para")
  
  # The scale parameter eta_f is the boundary. Equal to 1 is the condition.
  old.df <- 100; new.df <- Estdf(e.t)
  new.YITA <- YITAtQMLE(e=e.t,dfest=new.df) # YITAtQMLE function
  iter <-  1; diff <- 1
  while ((iter<max.iter) & (new.df != old.df)) {
    iter <-  iter+1
    Est.model <- QMLE(series = series, LogLFunc = "LogL_GARCH_t", order = order, dfest = new.df)
    e.t <- Est.model$e
    new.para <- Est.model$QMLE.t
    new.diff <- sum((para-new.para)^2)/(p+q+1)
    Estm <- rbind(Estm, c(Est.model$QMLE.t, new.df, new.YITA, new.diff))
    old.df <- new.df; new.df <- Estdf(e.t)
    new.YITA <- YITAtQMLE(e=e.t,dfest=new.df) # This is the crux! Using the assumed df to estimate yita!
    diff <- new.diff; para <- new.para
  }
  return(Estm)
}

#一种内部函数，用于计算GARCH模型的残差，以便进一步的自适应过程可以向前推进
com.residue <- function(alpha, beta, series){
  n <- length(series)
  sig2 <- numeric(n)
  q <- length(alpha)-1
  p <- length(beta)
  d <- max(q,p)
  if (alpha[1]==0) sig2[1:d] <- abs(series[1:d])
  else sig2[1:d] <- alpha[1]/(1.0-sum(alpha[2:(q+1)])-sum(beta))
  for (t in (d+1):n){
    sig2[t] <-  sum(alpha * c(1, series[t - (1:q)]^2)) + sum(beta * sig2[t - (1:p)])
    # It is changed based on my NOTES in red GARCH(P,Q)
  }
  e <- series/(sqrt(sig2))
  list(sig.sq = sig2, e=e)
}

#该函数是全自适应拟最大似然估计的关键。
#它在误差分布族的特定假设下，根据模型残差估计自由度。所以输入只是模型残差。
Estdf <- function(e){
  n <- length(e)
  like <- function(e){
    df.tL <- function(par){
      df <- par[1]
      #tQMLE.e=e/yi
      yi <- 1
      if (df>0)
      {
        #f=gamma((df+1)/2)/((pi*df)^0.5*gamma(df/2))*((df/(df-2))^0.5)*(1+tQMLE.e^2/(df-2))^(-(df+1)/2)
        f <- gamma((df+1)/2)/((pi*df)^0.5*gamma(df/2))*(1+(e/yi)^2/(df))^(-(df+1)/2)
        sum(log(yi)-log(f))/n
      }
      else Inf
    }
    df.tL
  }
  df.qmle <- nlminb(c(0.01),like(e),lower=c(0),upper=Inf)
  return(round(df.qmle$par))
}

# Example to generate GARCH(1,1) with normal innovation.
GARCH1_1 <- function(n, a, b){
  # param n Number of series
  # a Alpha vector, which includes intercept and alpha1
  # b Beta vector
  
  e <- rnorm(n)
  x <- double(n)
  sig2t <- x
  x[1:1] <- rnorm(1, sd = sqrt(a[1]/(1.0-a[2]-b[1])))
  
  x[1]
  sig2t[1] <- a[1]/(1.0-a[2]-b[1])
  sig2t[1]
  for(i in 2:n) {
    # Generate GARCH(1,1) process
    sig2t[i] = sum(a * c(1, x[i-1]^2)) + sum(b[1] * sig2t[i - 1])
    x[i] = e[i]*sqrt(sig2t[i])
  }
  list(series = x, e = e, sig2t = sig2t)
}

# Return Log-likelihood function of GARCH model with Pearson's Type IV (PIV) distributed error
LogL_GARCH_PIV <- function(series, p, q, lambda, a, nu, m){
  GARCH_e_PIV <- function(para){
    n <- length(series)
    sig2 <- numeric(n)
    alpha <- para[1:(q+1)]
    beta <- para[(q+2):(q+1+p)]
    d <- max(q,p)
    sig2[1:(d)] <- alpha[1]/(1.0-sum(alpha[2:(q+1)])-sum(beta))
    for (t in (d+1):n) {
      sig2[t] <-  sum(alpha * c(1, series[t - (1:q)]^2)) + sum(beta * sig2[t - (1:p)])
    }
    if (c(alpha,beta) > rep(0,(q+p+1)) && (sum(alpha[2:(q+1)])+sum(beta)) < 1 && sig2>=0)
    {
      e <- series / (sqrt(sig2))
      # Here is an issue!!!
      # cmplx <- complex(real = m, imaginary = nu/2)
      # K <- 2^(2*m-2)*abs(gammaz(cmplx))^2/(a*pi*gammaz(2*m-1))
      # f <- K * (1 + ((e - lambda)/a)^2)^(-m) * exp(-nu * atan((e - lambda)/a))
      # return(sum(log(sqrt(sig2))-log(f)))
      
      # Perfect! Pay attention: the sqrt(sig2) is the a_t in the K, which is conditional standard deviation.
      L <- sum(log(sqrt(sig2)) + m * log(1+(e^2)) +nu * atan(e))
      return(L)
    }
    else return(Inf)
  }
  GARCH_e_PIV
}

# The Log-likelihood function of GARCH(p,q) model.
LogL_GARCH_Norm <- function(series, p, q) {
  GARCH_Norm <- function(para){
    n <- length(series)
    sig2 <- numeric(n); # tem1 <- numeric(1); tem2 <- numeric(1)
    alpha <- para[1:(q+1)] # w is the alpha(1)
    beta <- para[(q+2):(q+1+p)]
    d <- max(q,p)
    sig2[1:(d)] <- alpha[1]/(1.0-sum(alpha[2:(q+1)])-sum(beta))
    for (t in (d+1):n){
      sig2[t] = sum(alpha * c(1, series[t - (1:q)]^2)) + sum(beta * sig2[t - (1:p)])
      # It is changed based on my NOTES in red GARCH(P,Q)
    }
    if (all(c(alpha,beta) > rep(0,(q+p+1))) && (sum(alpha[2:(q+1)])+sum(beta)) < 1 && all(sig2>=0)){
      return(sum(2*log(sqrt(sig2)) + series^2/(sig2)))
    } else return(Inf)
  }
  GARCH_Norm
}

# The Log likelihood function of GARCH model. Here is the GARCH(p,q) with student's t error
LogL_GARCH_t <- function(series, p, q, df){
  # df The specified degree of freedome of student's t innovation.
  GARCH_e_t <- function(para){
    n <- length(series)
    sig2 <- numeric(n)
    alpha <- para[1:(q+1)]
    beta <- para[(q+2):(q+1+p)]
    d <- max(q,p)
    sig2[1:(d)] <- alpha[1]/(1.0-sum(alpha[2:(q+1)])-sum(beta))
    for (t in (d+1):n) {
      sig2[t] <-  sum(alpha * c(1, series[t - (1:q)]^2)) + sum(beta * sig2[t - (1:p)])
    }
    if (all(alpha > 0) && all(beta > 0) && (sum(alpha[2:(q+1)]) + sum(beta)) < 1 && all(sig2 >= 0))
    {
      tQMLE.e <- series / (sqrt(sig2))
      g <- (gamma((df + 1)/2) / (pi * df)^0.5 / gamma(df/2))*(1+tQMLE.e^2/(df))^(-(df+1)/2)
      return(sum(log(sqrt(sig2))-log(g)))
    }
    else return(Inf)
  }
  GARCH_e_t
}

# The Log-likelihood function of linear model. Here is the y=a + b*x. Need to be changed.
LogL_Linear_Norm <- function(X, y){
  Linear_Norm <- function(para){
    n <- nrow(X)
    m <- ncol(X)+1
    X <- cbind(rep(1,times = n), X)
    alpha <- para[1]
    beta <- para[2]
    s <- sd(y, na.rm = FALSE)
    e <- y - X %*% as.matrix(para, m, 1)
    return(n*log(s) + 1/(2*s^2) * sum(e^2))
  }
  Linear_Norm
}

# Simple Maximum Likelihood Estimation based on normal residual assumption.
MLE <- function(y, X, LogLFunc = c("LogL_GARCH_Norm", "LogL_Linear_Norm"), order = c(1,1)){
  ## normal distribution innovation likelihood
  # y The dependent variable or time series.
  # X The single independent variable for single linear regression. It could be missing if use this function to estimate GARCH(p,q) model.
  # LogLFunc The log-likelihood function of models.
  if (!missing(LogLFunc)) 
    LogLFunc <- match.arg(LogLFunc)
  else LogLFunc <- "LogL_GARCH_Norm"
  
  if (missing(X) || LogLFunc == "LogL_GARCH_Norm"){
    q <- order[1]; p <- order[2]
    ini.para <- rep(0.01, p+q+1)
    low.cons <- rep(0, p+q+1)
    up.cons <- c(Inf, rep(1, p+q))
    LogLFunc <- LogL_GARCH_Norm(y,p,q)
    
    MLE.N <- nlminb(ini.para, LogLFunc, lower=low.cons, upper=up.cons)
    list(MLE.N = MLE.N$par)
  } else if (!missing(X) && LogLFunc == "LogL_Linear_Norm") {
    LogLFunc <- LogL_Linear_Norm
    MLE.N <- nlminb(c(0.01,0.01), LogLFunc(X, y))
  }
}

# tQMLE is a function that can be used to estimate parameters of GARCH(p,q) with student's t innovation by
# specified log-likelihood estimation as "LogL_GARCH_t" and "dfest"(degree of freedom). It also includes
# QMLE other
QMLE <- function(series, LogLFunc = c("LogL_GARCH_Norm", "LogL_GARCH_t", "LogL_GARCH_PIV"), order = c(1,1), dfest, params.PIV){
  # sig2 in the formula is the sig^2!!!!!!!!!!!!!!pay attention.
  if (missing(series)) {
    stop("No input data for 'series'!")
  }
  q <- order[1]; p <- order[2]
  if (!missing(LogLFunc)) 
    LogLFunc <- match.arg(LogLFunc)
  else LogLFunc <- "LogL_GARCH_Norm"
  
  #LogLFunc is the log likelihood function of GARCH(p,q) with student's t innovation. The model is setted as
  # \eqn{\sigma^{2}_{t|t-1}=\omega+\sum_{i=1}^{q}\alpha_{i}u^{2}_{t-i}+\sum_{j=1}^{p}\beta_{j}\sigma^{2}_{t-j}}
  # This input can be changed according to assumption. Default setting is normal distributed innovation.
  
  if (LogLFunc == "LogL_GARCH_Norm") {
    # || (missing(dfest)) || (missing(params.PIV))) 
    LogLFunc <- LogL_GARCH_Norm
    QMLE.N <- MLE(y = series, LogLFunc(series), order = order)$MLE.N # MLE function
    pred <- com.residue(alpha = QMLE.N[1:(q+1)], beta = QMLE.N[(q+2):(p+q+1)], series=series)
    print("Estimated as normal error!")
    list(QMLE.N = QMLE.N, sigma.sq = pred$sig.sq, e = pred$e)
  } else if ((LogLFunc == "LogL_GARCH_t") && (!missing(dfest)))
  {
    df <- dfest
    ini.para <- rep(0.01, p+q+1)
    low.cons <- rep(0, p+q+1)
    up.cons <- c(Inf, rep(1, p+q))
    
    LogLFunc <- LogL_GARCH_t(series,p,q,df)
    QMLE.t <- nlminb(ini.para, LogLFunc, lower=low.cons, upper=up.cons)
    QMLE.t <- QMLE.t$par
    pred <- com.residue(alpha = QMLE.t[1:(q+1)], beta = QMLE.t[(q+2):(p+q+1)], series=series)
    print("Estimated as student's t error!")
    list(QMLE.t = QMLE.t, sigma.sq = pred$sig.sq, e = pred$e)
  } else if ((LogLFunc == "LogL_GARCH_PIV") && (!missing(params.PIV)))
  {
    LogLFunc <- LogL_GARCH_PIV
    ini.para <- rep(0.01, p+q+1)
    low.cons <- rep(0, p+q+1)
    up.cons <- c(Inf, rep(1, p+q))
    m <- params.PIV$m; nu <- params.PIV$nu; lambda <- params.PIV$location; a <- params.PIV$scale
    QMLE.PIV <- nlminb(ini.para, LogLFunc(series, p, q, lambda, a, nu, m), lower=low.cons, upper=up.cons)
    QMLE.PIV <- QMLE.PIV$par
    pred <- com.residue(alpha = QMLE.PIV[1:(q+1)], beta = QMLE.PIV[(q+2):(p+q+1)], series=series)
    print("Estimated as Pearson Type IV error!")
    list(QMLE.PIV = QMLE.PIV, sigma.sq = pred$sig.sq, e = pred$e)
  }
}

# This function is used to estimate scale parameter \eqn{\eta_f} from the model residues under specified assumption of distribution of error. 
# So the input includes degree of freedom and model residues.
YITAtQMLE <- function(e,dfest){
  n <- length(e)
  like <- function(e,df){
    yitatL <- function(par){
      yi <- par[1]
      #tQMLE.e=e/yi
      if (yi>0)
      {
        #f=gamma((df+1)/2)/((pi*df)^0.5*gamma(df/2))*((df/(df-2))^0.5)*(1+tQMLE.e^2/(df-2))^(-(df+1)/2)
        f <- gamma((df+1)/2)/((pi*df)^0.5*gamma(df/2))*(1+(e/yi)^2/(df))^(-(df+1)/2)
        sum(log(yi)-log(f))/n
      }
      else Inf
    }
    yitatL
  }
  yitqmle <- nlminb(c(0.01),like(e,df=dfest),lower=c(0),upper=Inf)
  return(yitqmle$par)
}

# 模板
data(EuStockMarkets)
dax <- diff(log(EuStockMarkets))[,"DAX"]
dax.garch <- garch(dax) # Fit a GARCH(1,1) to DAX returns
summary(dax.garch) # ARCH effects are filtered. However,

library(xts)
test <- dflnSP500[-1,]
h <- length(test)
plot(test)

alpha <- c(0.1, 0.5); beta <- c(0.2) # GARCH(1,1) coefficients
test.arch <- garch(test, order = c(1,1)) 
est1 <- MLE(y = test, LogLFunc = "LogL_GARCH_Norm", order = c(1,1))
est1

est2 <- QMLE(series = test, LogLFunc = "LogL_GARCH_t", order = c(1,1), dfest = 20)
est2$sigma.sq #输出条件方差
plot(est2$sigma.sq, type='l')

est3 <- QMLE(series = test, order = c(1,1))
summary(est3)
est3$QMLE.N



