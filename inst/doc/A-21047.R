## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(cars) + 
  geom_point(aes(x = speed, y = dist), size = 1, shape = 1) + 
  geom_smooth(aes(x = speed, y = dist), method = "lm", formula = "y ~ x") + 
  labs(title = "Stopping Distance versus Speed", x = "speed", y = "dist")

## -----------------------------------------------------------------------------
data(cars)
speed <- cars$speed
dist  <- cars$dist
carlm <- lm(dist ~ speed)
summary(carlm)$coef

## -----------------------------------------------------------------------------
summary(cars)
plot(speed, dist, type = "p", pch = 1, cex = 0.8,
     main = "Stopping Distance versus Speed",
     xlab = "speed", ylab = "dist") 
abline(carlm, col = "blue", lwd = 2, lty = 1)

## -----------------------------------------------------------------------------
#par(mfrow = c(2,2))
plot(carlm)

## -----------------------------------------------------------------------------
library(xtable)
xtable::xtable(summary(cars))

## ----results = 'asis'---------------------------------------------------------
print(xtable(summary(cars)), type = "html", comment = FALSE)

## ---- fig.width=8,fig.height=5------------------------------------------------
cRayleigh = function(sigma=1){
  n <- 1000
  u <- runif(n)
  x <- sqrt(-2*sigma^2*log(1-u)) # F(x)= 1 - e^{-x^2/(2*sigma^2)}, 0<=x$
  hist(x, prob = TRUE, main = expression(f(x)== frac(x,sigma^2)*e^{-frac(x^2, 2*sigma^2)}),sub=paste("sigma = ", sigma))
  y <- seq(0, 100, .01)
  lines(y, y/sigma^2*exp(-y^2/(2*sigma^2)))  
}
#par(mfrow = c(2,2))
cRayleigh(0.1)
cRayleigh(1)
cRayleigh(6)
cRayleigh(10)

## ---- fig.width=8,fig.height=5------------------------------------------------
cNorm = function(p1=.75){
  n <- 1e3;
  j <- k <- 0;
  y <- numeric(n)

  while (k < n) {
    u <- runif(1) # Acceptance R.V.
    j <- j + 1
    x <- runif(1,-10,10) #random variate from g(.) 
    
    if (p1*dnorm(x,0,1) + (1-p1)*dnorm(x,3,1) > u) {
      #we accept x
      k <- k + 1
      y[k] <- x 
    }
  }
  hist(y, prob = TRUE, main = paste("p1 =", p1))
  y <- seq(-100, 100, .01)
  lines(y, p1*dnorm(y,0,1) + (1-p1)*dnorm(y,3,1))
}
#par(mfrow = c(2,3))
cNorm(0.0)
cNorm(0.1)
cNorm(0.25)
cNorm(0.5)
cNorm(0.75)
cNorm(1.0)

## -----------------------------------------------------------------------------
cGammaPois = function(n = 1e3, t = 10, lambda = 1,alpha = 1, beta = 1){
  # n: the total number
  # t: time
  # lambda: for poisson lambda
  # alpha, beta: for gamma param
  x1 <- 0
  for(i in seq(1,n)){
    N <- sum(rpois(t,lambda)) # Generate N(t) for every observation.
    data <- rgamma(N,shape = alpha,rate = beta)
    x1[i] <- sum(data)
  }
  return(x1)
  # Generate Gamma R.V.
  # data <- array(rgamma(n*N,shape = alpha, rate = beta),dim = c(N,n))
}
x <- cGammaPois(1e4,10,10,0.3,3)
hist(x, prob = TRUE)

## -----------------------------------------------------------------------------
n      <- 1e3 # 1e2, 1e3, 1e4
t      <- 10
lambda <- 1   # 0.1, 1, 10
alpha  <- 1   # 0.1, 1, 10
beta   <- 1   # 0.1, 1, 10

cCheckMeanVar = function(n = 1e3, t = 10, lambda = 1, alpha = 1, beta = 1){
  x <- cGammaPois(n,t,lambda,alpha,beta)
  mean_act <- mean(x)
  mean_thm <- lambda * t * alpha / beta
  var_act  <- var(x)
  var_thm  <- lambda * t * (alpha + alpha^2 ) / (beta^2)
  comD     <- data.frame(Mean = c(mean_act,mean_thm), Var  = c(var_act ,var_thm), row.names = c("Est", "Thm"))
  knitr::kable(comD,caption = paste("n=",n,"t=",t,"lambda=",lambda,"alpha=",alpha,"beta=",beta))
}

## -----------------------------------------------------------------------------
cCheckMeanVar(1e2,10,1,1,1)
cCheckMeanVar(1e3,10,1,1,1)
cCheckMeanVar(1e4,10,1,1,1)

## -----------------------------------------------------------------------------
cCheckMeanVar(1e3,10,0.1,1,1)
cCheckMeanVar(1e3,10,1,1,1)
cCheckMeanVar(1e3,10,10,1,1)

## -----------------------------------------------------------------------------
cCheckMeanVar(1e3,10,1,0.1,1)
cCheckMeanVar(1e3,10,1,1,1)
cCheckMeanVar(1e3,10,1,10,1)

## -----------------------------------------------------------------------------
cCheckMeanVar(1e3,10,1,1,0.1)
cCheckMeanVar(1e3,10,1,1,1)
cCheckMeanVar(1e3,10,1,1,10)

## -----------------------------------------------------------------------------
MC.beta <- function(n = 10000, start = 0, end = 1){
  x <- runif(n,start,end)
  hx <- 30*x^2 *(1-x)^2 *(end - start)
  return(mean(hx))
}

# MC.beta(end = 1) # for all x, namely F(1) I put them in the for loop

bef <- array(0,dim = c(11,3))
for (le in seq(1:11)) {
  xi <- (le - 1) / 10
  bef[le,1] <- xi 
  bef[le,2] <- MC.beta(end = xi)
  bef[le,3] <- pbeta(xi, 3, 3)
}
knitr::kable(bef, col.names = c("x","MC.beta","pbeta"), caption = paste("alpha = 3","beta = 3"))

## -----------------------------------------------------------------------------
# Generate hat theta with m = n
antiRayleigh <- function(x, sigma = 1, n = 1e3, anti = TRUE){
  u <- runif(n/2,0,x)
  if(anti){ 
    v = 1 - u
  }
  else {
    v <- runif(n/2,0,x)
  }
  u <- c(u,v)
  cdf <- numeric(length(x))

  for (i in 1:length(x)) {
    g <- x[i]/sigma^2*u*exp(-u^2/2/sigma^2)
    cdf[i] <- mean(g)
  }
  return(cdf)
}

# m samples for Rayleigh(x)
m <- 1000
x <- 1.95
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
  MC1[i] <- antiRayleigh(x, sigma = 2, anti = FALSE)
  MC2[i] <- antiRayleigh(x, sigma = 2)
}

knitr::kable(
  data.frame(
    sdMC = c(sd(MC1)),
    sdMCanti = c(sd(MC2)),
    MCantiverseMC = c(sd(MC2)/sd(MC1)))  
)

## -----------------------------------------------------------------------------
# compare the reduction with n = 2
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
  MC1[i] <- antiRayleigh(x, sigma = 2, n = 2,anti = FALSE)
  MC2[i] <- antiRayleigh(x, sigma = 2, n = 2)
}

knitr::kable(
  data.frame(
    sdMC = c(sd(MC1)),
    sdMCanti = c(sd(MC2)),
    MCantiverseMC = c(sd(MC2)/sd(MC1)))  
)

## -----------------------------------------------------------------------------
Eg <- 1/sqrt(2* pi * exp(1)) + pnorm(-1)
Eg1square <- pnorm(-1) * (3 * pnorm(-1)+4/sqrt(2 *pi *exp(1)))
Eg2square <- 3/2/pi/exp(1)
Varg1 <- Eg1square - Eg ** 2
Varg2 <- Eg2square - Eg ** 2
knitr::kable(
  data.frame(
    Eg = c(Eg),
    Varg1 = c(Varg1),
    Varg2 = c(Varg2)))

## -----------------------------------------------------------------------------
# generate f1 and f2 distribution using acceptance method
f1 <- function(x){
  y <- dnorm(x)/pnorm(-1)
  return(y)}
f2 <- function(x){
  y <- sqrt(exp(1))*x*exp(-x**2/2)
  return(y)}

# n samples from f1 and f2
fx <- function(n = 1e3){
  x1 <- x2 <- numeric(n)
  res <- matrix(ncol = n,nrow = 2)
  i  <- j  <- 1
  while (i <= n || j <= n) {
    y  <- runif(1)
    x <- runif(1,1,n)
    if(f1(x) > y && i<= n ){
      # we accept x
      x1[i] <- x
      i <- i + 1
    }
    if(f2(x) >y && j<= n){
      x2[j] <- x
      j = j + 1
    }
  }
  res[1,] <- x1
  res[2,] <- x2
  return(res)
}

# calculate one g1 and g2 with n samples
g <- function(n = 1e3){
  x <- fx(n)
  
  y1  <- x[1,]^2*pnorm(-1)
  me1 <- mean(y1)
  sd1 <- mean(y1^2) - me1^2
  
  y2  <- x[2,]/sqrt(2*pi*exp(1))
  me2 <- mean(y2)
  sd2 <- mean(y2^2) - me2^2

  res = c(me1,sd1,me2,sd2)
  return(res)
}

# calculate m gs with m groups
m <- 2 # change it to a larger one if you want to get concise data
#mafor5.14 <- matrix(ncol = 4,nrow = m)
#for (j in seq(1:m)) {
#  mafor5.14[j,] <- g()
#}
load("../data/mafor5.14.rda")
knitr::kable(
  data.frame(
    g1_thm = c(Eg, Varg1/m),
    g1_act = c(mean(mafor5.14[,1]),mean(mafor5.14[,2])/m),
    g2_thm = c(Eg, Varg2/m),
    g2_act = c(mean(mafor5.14[,3]),mean(mafor5.14[,4])/m),
    row.names = c("Mean","Var")
    ))

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
m = 1e3
LCL <- UCL <- numeric(m)
for (i in seq(1:m)){
  x <- rnorm(n,mean = 0,sd=2)
  LCL[i] <- mean(x) - sqrt(var(x)/n)*qt(1-alpha/2, n-1)
  UCL[i] <- mean(x) - sqrt(var(x)/n)*qt(alpha/2, n-1)  
}
# count the number of intervals that contain mu=0
sum( LCL < 0 & UCL > 0 )
# or compute the mean to get the confidence level
mean( LCL < 0 & UCL > 0 )

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
m = 1e3
LCL <- UCL <- numeric(m)
for (i in seq(1:m)){
  x <- rchisq(n,df=2)
  LCL[i] <- mean(x) - sqrt(var(x)/n)*qt(1-alpha/2, n-1)
  UCL[i] <- mean(x) - sqrt(var(x)/n)*qt(alpha/2, n-1)  
}
# count the number of intervals that contain mu=0
sum( LCL < 2 & UCL > 2 )
# or compute the mean to get the confidence level
mean( LCL < 2 & UCL > 2 )

## -----------------------------------------------------------------------------
ex6acp <- function(m =1e4, method = "norm"){
  n     <- 20
  alpha <- .05
  mu0   <- 1
  p     <- numeric(m) #storage for p-values
  for (j in 1:m) {
    x <- switch (method,
                 norm  = rnorm(n, mean = 1, sd = 2), # N(1,2)
                 chiq1 = rchisq(n, df = 1), # Chi(1)
                 unif  = runif(n, min = 0, max = 2), # U(0,2)
                 expo1 = rexp(n, rate = 1) # Exp(1)
                 )
    ttest <- t.test(x, alternative = "two.sided", mu = mu0)
    p[j]  <- ttest$p.value
    }
  p.hat  <- mean(p < alpha)
  se.hat <- sqrt(p.hat * (1 - p.hat) / m)
  print(c(p.hat, se.hat))
}
Norm  <- ex6acp(method = "norm") # samples from N(1,2)
Chiq1 <- ex6acp(method = "chiq1") # samples from Chiq(1)
Unif  <- ex6acp(method = "unif") # samples from Uni(0,2)
Expo1 <- ex6acp(method = "expo1") # samples from Expo(1)
res   <- data.frame(Item = c("p.hat","se.hat"),Norm, Chiq1, Unif, Expo1)
knitr::kable(res)

## -----------------------------------------------------------------------------
# generate samples from N(0,1) to form matrix as sm
sm_ge <- function(p = 1, n = 10, e = 0){
  sm <- matrix(nrow = p, ncol = n)
  
  sigma <- sample(c(1, 10), replace = TRUE, size = n, prob = c(1-e, e))
  for (pj in 1:p) {
      sm[pj,] <- rnorm(n, 0, sigma)
    }
  return(sm)
}

# calculate the sk of sm
sk_cal <- function(sm){
  p <- dim(sm)[1]
  n <- dim(sm)[2]
  xbarvec <- rowMeans(sm)
  
  # calculate the SigmaHat of p * p
  sh <- matrix(0,ncol = p, nrow = p)
  for (i in 1:n){
    sh <- sh + (sm[,i] - xbarvec) %*% t(sm[,i] - xbarvec)
    #we couldn't use the following formula because cov() uses n-1.
    # sigmahat[i,j] <- cov(sm[i,],sm[j,])
    }
  sigmahat <- solve(sh / n)

  sum <- 0
  for (i in 1:n){
    for (j in 1:n) {
      sum <- sum + (t(sm[,i]-xbarvec) %*% sigmahat %*% (sm[,j]-xbarvec))^3
    }
  }
  # sum
  sk <- sum / n^2
  return(sk)
}

check_sim <- function(p =1, n = 10, m = 1e3, alpha = 0.05, e = 0){
  cv <- qchisq(1-alpha, df = p*(p+1)*(p+2)/6 )

  # this takes quite long time pls wait for minites.
  sktests <- numeric(m)
  for (k in 1:m) {
    #test decision is 1 (reject) or 0
    sktests[k] <- as.integer(n/6*sk_cal(sm_ge(p, n, e)) >= cv)
    }
  p.reject <- mean(sktests) #proportion rejected
  return(p.reject)
}

## -----------------------------------------------------------------------------
n <- c(10,20,30,50,100) # it takes too long so I abandon 500
# resfor6_8 <- numeric(length(n))
#for (nj in 1:length(n)) {
#  resfor6_8[nj] <- check_sim(1,n[nj])
#}
load("../data/resfor6_8.rda")
knitr::kable(data.frame(n,resfor6_8))

## -----------------------------------------------------------------------------
epsilon <- c(seq(0,0.1,0.02), seq(.2, 1, .1))
#epsilon <- c(seq(0,1,0.1))
N <- length(epsilon)
# power <- numeric(N)
#for (j in 1:N) {
#  #for these epsilons 
#  power[j] <- check_sim(p =1, n=30, m=1e3, alpha = 0.1, epsilon[j])
#  #print(power)
#}
load("../data/powerfor6_8.rdata")
print(powerfor6_8)
#plot power vs epsilon 
plot(epsilon, powerfor6_8, type = "b",
     xlab = bquote(epsilon), 
     ylab = "power",
     xlim = c(0,1),
     ylim = c(0,1))
abline(h = 0.1, lty = 2.5)
se <- sqrt(powerfor6_8 * (1-powerfor6_8) / 1e3)
#standard errors 
lines(epsilon, powerfor6_8+se, lty = 3)
lines(epsilon, powerfor6_8-se, lty = 3)

## -----------------------------------------------------------------------------
load("../data/classscores.rda")
data <- data.frame(Mechanics_C=classscores[,1],
                   Vectors_C=classscores[,2],
                   Algebra_O=classscores[,3],
                   Analysis_O=classscores[,4],
                   Statistics_O=classscores[,5])
# I read the [84] An Introduction to the Bootstrap by Bradley Efron, Robert J. Tibshirani, p63, and I find he used the bias the estimator of cov with denominator n. However, in r package cov uses denominator n - 1 to ensure unbiased.
cp_cov <- function(data, bycol = TRUE){
  if(bycol!=TRUE){data <- t(data)}
  n <- dim(data)[2] # col
  m <- dim(data)[1] # row

  covc <- matrix(ncol = n,nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
        covc[i,j] <- (m-1)/m*cov(data[,i],data[,j])
    }
  }
  return(covc)
}

# We developed a function to calculate the eigenvalues.
theta_func <- function(data,demo_by_n = FALSE){
  data_cov <- cov(data) # calculate the covariance matrix
  if(demo_by_n == TRUE){data_cov <- cp_cov(data)}
  ev       <- eigen(data_cov)$values # calculate the eigenvalues of cov matrix
  theta    <- ev[1]/sum(ev) # return the eigenvalues sorted in decreasing order
  return(c(theta,ev))
}

round(theta_func(data,1)[2:6],1) # denominated by n
# we first calculate the same values as reported.

theta_func(data,1)[1] # denominated by n
theta_func(data,0)[1] # denominated by n-1
# we find the estimator of theta is equal, haha

## -----------------------------------------------------------------------------
# data already defined as data.frame above
estimate_without_packages <- function(data,seed=54321,N=1e3){
  set.seed(seed)
  thetastar <- numeric(N)
  theta <- theta_func(data)[1]
  for (n in 1:N) {
    datastar <- data[sample(nrow(data),replace =TRUE),]
    thetastar[n] <- theta_func(datastar)[1]
  }
  ewop <- round(c(ori=theta,bias_bs=mean(thetastar) - theta, se_bs=sd(thetastar)),3)
  return(ewop)
}

library(boot)
library(MASS)
only_theta_func <- function(data, i){
  data <- data[i,]
  data_cov <- cov(data) # calculate the covariance matrix
  ev       <- eigen(data_cov)$values # calculate the eigenvalues of cov matrix
  theta    <- ev[1]/sum(ev) # return the eigenvalues sorted in decreasing order
  return(c(theta))
}

estimate_with_packages <- function(data,seed=54321,N=1e3){
  set.seed(seed)
  res <- boot(data,statistic = only_theta_func,R = N)
  ewp <- round(c(ori=res$t0,bias_bs=mean(res$t)-res$t0,se_bs=sd(res$t)),3)
  return(ewp)
}

estimate_without_packages(data)
estimate_with_packages(data)

## -----------------------------------------------------------------------------
cp_JK <- function(data){
  n <- nrow(data)
  theta_hat <- only_theta_func(data,1:n)
  theta_jk <- numeric(n)
  for (i in 1:n) {
    theta_jk[i] <- only_theta_func(data,(1:n)[-i])
    }
  bias_jk <- (n-1)*(mean(theta_jk)-theta_hat)
  se_jk   <- sqrt((n-1)*mean((theta_jk-theta_hat)^2))
  ewpres <- estimate_with_packages(data)
  res <- round(c(ori=theta_hat,bias_jk=bias_jk,se_jk=se_jk,ewpres[3]),3)
  return(res)
}
cp_JK(data)

## -----------------------------------------------------------------------------
bs_BCa <- function(data, theta, statistic, conf = .95) {
  data <- as.matrix(data) 
  n <- nrow(data);N <- 1:n 
  alpha <- (1 + c(-conf, conf))/2;zalpha <- qnorm(alpha)
  # the bias correction factor 
  theta0 <- theta_func(data)[1];  z0 <- qnorm(sum(theta < theta0) / length(theta))
  # the acceleration factor (jackknife est.) 
  theta_jk <- numeric(n)
  for (i in N) {
    J <- N[1:(n-1)]
    theta_jk[i] <- statistic(data[-i, ], J) 
  }
  Loss <- mean(theta_jk) - theta_jk
  adj <- sum(Loss^3)/(6 * sum(Loss^2)^1.5)

  # BCa conf. limits 
  alpha_adj <- pnorm(z0 + (z0+zalpha)/(1-adj*(z0+zalpha))) 
  Q <- quantile(theta, alpha_adj, type=6)
  res <- as.matrix(Q)
  return(res)
}

ewp_data <- function(data,seed=54321,N=1e3){
  set.seed(seed)
  thetastar <- numeric(N)
  for (n in 1:N) {
    datastar <- data[sample(nrow(data),replace =TRUE),]
    thetastar[n] <- theta_func(datastar)[1]
  }
  return(thetastar)
}
N <- 1e3
theta_boot <- ewp_data(data,N)
res <- bs_BCa(data,theta_boot,stat=only_theta_func)
# 1st col percentile values, 2nd col the BCa
res
alpha_res <- sum(res[1,1] <= theta_boot & theta_boot <= res[2,1])/N
alpha_res

## -----------------------------------------------------------------------------
# define sk and sk_boot to calculate the value out.
sk <- function(x,justi=FALSE){
  n <- length(x)
  if (justi == TRUE){
    skew <- 1/n*sum((x-mean(x))^3)/(sd(x)*sqrt((n-1)/n))^3
  } else {
    skew <- 1/n*sum((x-mean(x))^3)/sd(x)^3
  }
  return(skew)
}

boot_skew <- function(x,i){
  skew <- sk(x[i],TRUE)
  return(skew)
}

## -----------------------------------------------------------------------------
library(boot)
library(progress)

check_bs_ci <- function(sample_size=500, mc = 1e3, func="norm"){
  bs_ci_norm <- bs_ci_basi <- bs_ci_perc <- matrix(NA, mc, 2)
  pb <- progress_bar$new(total = mc)
  msk <- switch (func,
    norm  = 0,
    chisq = 2*sqrt(2/5)
  )
  for (m in 1:mc){
    x <- switch (func,
      norm  = rnorm(sample_size),
      chisq = rchisq(sample_size,df=4) # n = 5, df = 4
    )
    bs    <- boot(x,statistic = boot_skew, R=sample_size) # 1 sampling with sample_size
    bs_ci <- boot.ci(bs,type = c("norm","basic", "perc"))
    bs_ci_norm[m,] <- bs_ci$norm[2:3]
    bs_ci_basi[m,] <- bs_ci$basic[4:5]
    bs_ci_perc[m,] <- bs_ci$percent[4:5]
    pb$tick() #if (m %% (mc/10) == 0){print(cat("loading =", (m/mc)*100,"%"))}
  }
  CI_norm <- c(mean(bs_ci_norm[,1]),mean(bs_ci_norm[,2]))
  CI_norm_alpha_L <- mean(bs_ci_norm[,1] > msk)
  CI_norm_alpha_R <- mean(bs_ci_norm[,2] < msk)
  CI_norm_alpha   <- mean(bs_ci_norm[,1] <= msk & bs_ci_norm[,2] >= msk)
  
  CI_basi <- c(mean(bs_ci_basi[,1]),mean(bs_ci_basi[,2]))
  CI_basi_alpha_L <- mean(bs_ci_basi[,1] > msk)
  CI_basi_alpha_R <- mean(bs_ci_basi[,2] < msk)
  CI_basi_alpha   <- mean(bs_ci_basi[,1] <= msk & bs_ci_basi[,2] >= msk)
  
  CI_perc <- c(mean(bs_ci_perc[,1]),mean(bs_ci_perc[,2]))
  CI_perc_alpha_L <- mean(bs_ci_perc[,1] > msk)
  CI_perc_alpha_R <- mean(bs_ci_perc[,2] < msk)
  CI_perc_alpha   <- mean(bs_ci_perc[,1] <= msk & bs_ci_perc[,2] >= msk)
  
  df <- data.frame(ID=c("Norm","Basi","Perc"),
             Alpha_L = c(CI_norm_alpha_L,CI_basi_alpha_L,CI_perc_alpha_L),
             Alpha_M = c(CI_norm_alpha  ,CI_basi_alpha  ,CI_perc_alpha  ),
             Alpha_R = c(CI_norm_alpha_R,CI_basi_alpha_R,CI_perc_alpha_R),
             CI_L    = c(CI_norm[1]     ,CI_basi[1]     ,CI_perc[1]     ),
             CI_R    = c(CI_norm[2]     ,CI_basi[2]     ,CI_perc[2]     ))
  print(df)
  return(df)
}
#resfor7Bchisq <- check_bs_ci(func = "chisq")
load("../data/resfor7Bchisq.rda")
#resfor7Bnorm <- check_bs_ci(func = "norm")
load("../data/resfor7Bnorm.rda")

## -----------------------------------------------------------------------------
# the cor() function can only get statistic, so we use cor.test() to get more parameters.
cp_check_indep <- function(x,y){
  Re  <- 1e3 - 1 # notice the denominator is Re + 1
  z   <- c(x,y)
  NAl <- 1:length(z)
  Nx  <- length(x)
  reps_s <- reps_d <- numeric(Re) 
  # d for default, s for spearman
  res_d <- cor.test(x,y)
  res_s <- cor.test(x,y,method = "spearman")
  
  for (i in 1:Re) {
    k <- sample(NAl, size = Nx, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k]
    reps_d[i] <- cor.test(x1,y1)$statistic
    reps_s[i] <- cor.test(x1,y1,method = "spearman")$statistic
  }
  # Notice the nominator is 1 + #{theta_hat >= theta}
  p_d <- mean(abs(c(res_d$statistic,reps_d)) >= abs(res_d$statistic))
  p_s <- mean(abs(c(res_s$statistic,reps_s)) >= abs(res_s$statistic))
  print(round(c(p_d,res_d$p.value,p_s,res_s$p.value),3))
  print(res_s$statistic)
}

n <- 1e2
set.seed(5431)
All <- rnorm(2*n); x <- All[1:n]; y <- All[-(1:n)]

res1 <- cp_check_indep(x,y)
res2 <- cp_check_indep(sort(x),sort(y,TRUE))
res3 <- cp_check_indep(sort(x),sort(y))

## -----------------------------------------------------------------------------
library(boot)

# NN test function
library(RANN)
Ran <- function(z, ix, sizes, k) {
  if(is.vector(z)) {
    z <- data.frame(z,0)
  }
  ss <- sizes;
  z <- z[ix, ];
  NN <- nn2(data = z, k = k+1)
  b1 <- NN$nn.idx[1:ss[1],-1] 
  b2 <- NN$nn.idx[(ss[1]+1):(ss[1]+ss[2]),-1]
  res <- (sum(b1 < ss[1] + 0.5) + sum(b2 > ss[2] + 0.5)) /
    (k * (ss[1] + ss[2]))
}
nn_test <- function(z,Nxy,Re = 9999,k=3){
  nn_obj <- boot(data = z, statistic = Ran, R = Re,
                   sim = "permutation", sizes = Nxy, k=k)
  res <- c(nn_obj$t0,nn_obj$t)
  p_val_nn <- mean(res >= res[1])
  list(statistic = res[1],p.value = p_val_nn)
  return(p_val_nn)
}

# energy test function
library(energy)
energy_test <- function(z,Nxy,Re=9999){
  energy_obs <- eqdist.etest(z, sizes=Nxy, R = Re)
  p_val_energy <- energy_obs$p.value
  return(p_val_energy)
}

# Ball test function
# install.packages("Ball")
library(Ball)
ball_test <- function(x,y,Re=9999, seed = 12345){
  p_val_ball <- bd.test(x,y, num.permutations = Re,
                        seed = seed)$p.value
  return(p_val_ball)
}

## -----------------------------------------------------------------------------
library(progress)
set.seed(12345)
k <- 3
n1 <- n2 <- 50
n  <- n1+n2
Nxy <- c(n1,n2)

# implement an function to generate x,y,z for every i in m. This part is just for testing, we will change it later.
p <- 2
xyz <- function(){
  x <- matrix(rnorm(n1*p),ncol=p)
  y <- matrix(rnorm(n2*p),ncol=p)
  z <- rbind(x,y)
  return(z)
}

# the final function
final_test <- function(m = 1e3, Re = 999){
  p.val <- matrix(NA,m,3)
  pb <- progress_bar$new(total = m)
  
  for(i in 1:m){
    z <- xyz()
    p.val[i,1] <- nn_test(z,Nxy,Re,k)
    p.val[i,2] <- energy_test(z,Nxy,R=Re)
    p.val[i,3] <- ball_test(z[1:n1,],z[-(1:n1),],Re = Re,seed = i*12345)
    pb$tick()
  }
  alpha <- 0.1;
  power <- colMeans(p.val < alpha)
  print(power)
  return(power)
}

## -----------------------------------------------------------------------------
p <- 1 # you can change p to 2 if you want
xyz <- function(){
  x <- matrix(rnorm(n1*p,0,1), ncol = p);
  y <- matrix(rnorm(n2*p,0,1.6), ncol = p);
  z <- rbind(x,y)
  return(z)
}
#r1 <- final_test(m = 1e2, Re = 999) 
load("../data/r1forQs.rdata")
r1forQs

## -----------------------------------------------------------------------------
p <- 1 # you can change p to 2 if you want
xyz <- function(){
  x <- matrix(rnorm(n1*p,0,1), ncol = p);
  y <- matrix(rnorm(n2*p,0.2,1.6), ncol = p);
  z <- rbind(x,y)
  return(z)
}
#r2 <- final_test(m = 1e2, Re = 999) 
# You can change 1e2 to 1e3, if you want more precise result.
load("../data/r2forQs.rdata")
r2forQs
# for N(0.6,1.6^2)
xyz <- function(){
  x <- matrix(rnorm(n1*p,0,1), ncol = p);
  y <- matrix(rnorm(n2*p,1,1.6), ncol = p);
  z <- rbind(x,y)
  return(z)
}
#r2.1 <- final_test(m=1e2)
load("../data/r2.1forQs.rdata")
r2.1forQs

## -----------------------------------------------------------------------------
p <- 1
rmix_norm <- function(n,alpha=0.5,mu1=0,mu2=0,sd1=1,sd2=1){
  # clot 0 then Normal 1 else Normal 2
  # alpha * Normal1 + (1 - alpha)* Normal2
  clot <- sample(c(0,1),size=n, replace = TRUE,prob = c(alpha,1-alpha))
  x <- rnorm(n,mu1,sd1)
  y <- rnorm(n,mu2,sd2)
  ret <- x
  for (i in 1:n) {if(clot[i] == 1){ret[i] <- y[i]}}
  return(ret)
}
xyz <- function(){
  x <- matrix(rt(n1*p,1), ncol = p);
  y <- matrix(rmix_norm(n2*p,alpha = 0.5, 0, 0, 1, 2), ncol = p);
  z <- rbind(x,y)
  return(z)
}
#r3forQs <- final_test(m=1e2)
load("../data/r3forQs.rdata")
r3forQs
xyz <- function(){
  x <- matrix(rt(n1*p,1), ncol = p);
  y <- matrix(rmix_norm(n2*p,alpha = 0.5,0, 1, 1, 2), ncol = p);
  z <- rbind(x,y)
  return(z)
}
#r3.1forQs <- final_test(m=1e2)
load("../data/r3.1forQs.rdata")
r3.1forQs

## -----------------------------------------------------------------------------
# since the y has 10 obs and paired with x
# there are other methods to perform this experiment
n1 <- 10
Nxy <- c(n1,10*n1)
xyz <- function(){
  # we need to generate x,y that are not indepent
  x <- matrix(rnorm(n1,10,1), ncol = 1);
  y <- matrix(x, nrow = n1, ncol = 10) + # baseline
    matrix(0.9,nrow = n1,ncol = 10) + # fixed effect
    matrix(rnorm(n1*10,0,1), ncol =10) # random effect
  z <- as.matrix(c(x,y))
  return(z)
}

final_test_modified <- function(m = 1e3, Re = 999){
  p.val <- matrix(NA,m,3)
  pb <- progress_bar$new(total = m)
  
  for(i in 1:m){
    z <- xyz()
    p.val[i,1] <- nn_test(z,Nxy,Re,k)
    p.val[i,2] <- energy_test(z,Nxy,R=Re)
    p.val[i,3] <- ball_test(z[1:n1,],z[-(1:n1),],Re = Re,seed = i*12345)
    pb$tick()
  }
  alpha <- 0.1;
  power <- colMeans(p.val < alpha)
  print(power)
  return(power)
}

#r4forQs <- final_test_modified(m = 1e2, Re = 999)
load("../data/r4forQs.rdata")
r4forQs

## -----------------------------------------------------------------------------
# return the cauchy density
cp_cauchy <- function(x,theta = 1, eta = 0){
  stopifnot(theta >0)
  return(1/(theta * pi * (1 + ((x-eta)/theta)^2)))
}

# version one with chisq
cp_mc_random_v1 <- function(length_of_chains = 1e4, x0 = NA,from = 1001){
  to <- length_of_chains + from - 1
  x <- numeric(to)
  xI <- x
  # Notice x is positive so we need xI to generate the negative value
  if (is.na(NA)){
    xI[1] <- x[1] <- rchisq(1, df = 100)
  } else {
    xI[1] <- x[1] <- x0
  }
  k <- 0
  u <- runif(to)
  
  for (i in 2:to){
    xt <- x[i-1]
    xIt <- xI[i-1]
    y  <- rchisq(1, df = xt)
    num <- cp_cauchy(y ) * dchisq(xt, df = y)
    den <- cp_cauchy(xt) * dchisq(y , df = xt)
    if (u[i] <= num / den ){
      x[i] <- y
      xI[i] <- y * sample(c(-1,1),size = 1) # generate negative or positive values
      } else {
        x[i] <- xt
        xI[i] <- xIt
        k <- k + 1
      }
  }
  return(c(k,xI[from:to]))
}
cp_mc_random_v2 <- function(length_of_chains = 1e4, x0 = NA,from = 1001){
  to = length_of_chains + from - 1
  x <- numeric(to)
  if (is.na(x0)){
    x[1] <- rt(1, df = 100)  
  } else {
    x[1] <- x0
  }
  
  k <- 0
  u <- runif(to)
  
  for (i in 2:to){
    xt <- x[i-1]
    y  <- rt(1, df = abs(xt))
    num <- cp_cauchy(y ) * dt(xt, df = abs(y))
    den <- cp_cauchy(xt) * dt(y , df = abs(xt))
    if (u[i] * den <= num ){
      x[i] <- y
      } else {
        x[i] <- xt
        k <- k + 1
      }
  }
  return(c(k,x[from:to]))
}
k1 <- cp_mc_random_v1()[1]
x1 <- cp_mc_random_v1()[-1]

k2 <- cp_mc_random_v2()[1]
x2 <- cp_mc_random_v2()[-1]

cp_show_chain <- function(x,from_p = 1, to_p = length(x)){
  idx <- from_p:to_p
  y1 <- x[idx]
  plot(idx,y1,type="l",ylab="x")
}

#par(mfrow=c(2,2))
cp_show_chain(x1)
cp_show_chain(x1,8500)
cp_show_chain(x2)
cp_show_chain(x2,8500)

cp_show_quantile <- function(x,show_all = FALSE){
  a <- ppoints(100)
  QR <- qcauchy(a)
  Q <- quantile(x,a)
  qqplot(QR,Q,xlab = "Std Cauchy Quantiles",ylab = "Sample Quantiles")

  if (show_all == FALSE){
    x_show <- x[x> min(QR) & x < max(QR)]
  } else {
    x_show <- x
  }
  hist(x_show,breaks = "Scott",freq = FALSE)
  lines(QR,cp_cauchy(QR))
}
#par(mfrow=c(2,2))
cp_show_quantile(x1)
cp_show_quantile(x2)

## -----------------------------------------------------------------------------
# this func is to check my calculation, not useful for this solution. Ignore it.
cp_ck <- function(n,a,b){
  sum1 <- 0
  for (i in 0:n){
    sum1 <- sum1 + choose(n+a+b-2,i+a-1)*choose(n,i)
  }
  sum1 <- sum1 * (n+a+b-1)
  sum2 <- choose(2*n+a+b-2,n+b-1) * (n+a+b-1)
  sum3 <- gamma(2*n+a+b-1)/(gamma(n+a)*gamma(n+b))* (n+a+b-1)
  return(c(sum1,sum2,sum3))
}
# cp_ck(4,1,1)

## -----------------------------------------------------------------------------
cp_gibbs <- function(length_of_chains = 1e4,from = 1001, a = 1, b = 1,x_range = 10, mu_xy = c(NA,NA)){
  # E(Y) = a / (a+b)
  # E(X) = n * y
  if (is.na(mu_xy[1]) & is.na(mu_xy[2])){ # x0 is not defined
    muy <- a/(a+b)
    mux <- x_range*muy
  } else { # x0 is defined
    muy <- mu_xy[2]
    mux <- mu_xy[1]
  }
  to = length_of_chains + from - 1
  Z <- matrix(0,nrow = to,ncol = 2)
  ### Start
  Z[1,] <- c(mux,muy) # Z initialize
  for (i in 2:to){
    # update X
    Zy <- Z[i-1,2]
    Z[i,1] <- rbinom(1,x_range,Zy)
    # update Y
    Zx <- Z[i,1]
    Z[i,2] <- rbeta(1,Zx+a,x_range-Zx+b)
  }
  return(Z[from:to,])
}

cp_show_gibbs_stat <- function(Z){
  Z_show <- Z
  col_mean_Z <- colMeans(Z_show)
  cov_Z <- cov(Z_show)
  cor_Z <- cor(Z_show)
  print(col_mean_Z)
  print(cov_Z)
  print(cor_Z)
  plot(Z_show,main = "",cex=0.5,xlab="x",ylab="y",ylim=range(Z_show[,2]))
}

Z <- cp_gibbs()
cp_show_gibbs_stat(Z)
#par(mfrow=c(2,2))
cp_show_chain(Z[,1],9500)
cp_show_chain(Z[,1])
cp_show_chain(Z[,2],9500)
cp_show_chain(Z[,2])

## -----------------------------------------------------------------------------
cp_GeRu_func <- function(psy){
  psy <- as.matrix(psy)
  num_col <- ncol(psy)
  num_row <- nrow(psy)
  
  mean_of_psy <- rowMeans(psy)
  var_bewteen_chains <- num_col * var(mean_of_psy)
  vars_within_chains <- apply(psy, 1, "var")
  var_estimate_within_chains <- mean(vars_within_chains)
  var_estimate_upper <- var_estimate_within_chains * (num_col-1)/num_col + var_bewteen_chains/num_col
  GeRu_stat <- var_estimate_upper/var_estimate_within_chains
  return(GeRu_stat)
}

# all generation for chains
cp_chains_generaton <- function(num_of_chains=4,length_of_chains=1e4, method="v1",x0=NA){
  x1 <- matrix(0,nrow = num_of_chains, ncol = length_of_chains)
  for (n_o_c in 1:num_of_chains){
    x1[n_o_c,] <- switch (method,
      "v1" = cp_mc_random_v1(length_of_chains,x0 = x0[n_o_c])[-1],
      "v2" = cp_mc_random_v2(length_of_chains,x0 = x0[n_o_c])[-1],
      "v3" = cp_gibbs(length_of_chains,mu_xy = c(x0[n_o_c,]))[,1],
      "v4" = cp_gibbs(length_of_chains,mu_xy = c(x0[n_o_c,]))[,2]
    )
  }
  return(x1)
}

#x1 <- cp_chains_generaton()

cp_psy_calculation <- function(x){
  # compute diagnostic statistics for x1 and x2
  psy <- t(apply(x1, 1, cumsum))
  for (i in 1:nrow(psy)){
    psy[i,] <- psy[i,] / (1:ncol(psy))
  }
  return(psy)
}

#psy <- cp_psy_calculation(x1)[1]

cp_plot_chains <- function(psy,length_of_chains=ncol(psy),num_of_chains=4){
  #par(mfrow=c(2,2))
  for (i in 1:4){
    plot(psy[i,],type="l",xlab=i,ylab=bquote(psy))
  }
  #par(mfrow=c(1,1))
}

# cp_plot_chains(psy)

cp_plot_r_hat <- function(psy,length_of_chains=ncol(psy)){
  rhat <- rep(0,length_of_chains)
  for (j in 1:length_of_chains){
    rhat[j] <- cp_GeRu_func(psy[,1:j])
  }
  plot(rhat[1:length_of_chains],type = "l",xlab="",ylab="R")
  abline(h=1.1,lty=2)
}

# Finally We draw this pic altogether
method = c("v1","v2","v3","v4")[-2] # ignore v2
GeRu_of_convergiance <- chain_length_of_convergiance <- numeric(length(method))

for (me in 1:length(method)){
  method_use = method[me]
  x0 <- switch (method_use,
    'v1' = c(NA,1,2,3),
    # 'v2' = c(NA,-1,1,10),
    'v3' = matrix(c(1,0.5,1,0.5,1,0.5,1,0.5),ncol = 2,nrow = 4),
    'v4' = matrix(c(1,0.5,1,0.5,1,0.5,1,0.5),ncol = 2,nrow = 4)
  )
  length_of_chains <- 1000

  for (turn in 1:10){
    x1 <- cp_chains_generaton(length_of_chains = length_of_chains, method = method_use,x0 = x0)
    psy <- cp_psy_calculation(x1)
    cp_dias <- cp_GeRu_func(psy)

    if (cp_dias < 1.2){
      chain_length_of_convergiance[me] <- length_of_chains
      print(length_of_chains)
      GeRu_of_convergiance[me] <- cp_dias
      print(cp_dias)
      cp_plot_chains(psy)
      cp_plot_r_hat(psy)
      break
    } else {
      length_of_chains <- length_of_chains + 1e3
    }
  }
}

## -----------------------------------------------------------------------------
# the kth value needs to be calculated seperately because each of them was too large. It can be calculated easily by twins.
L0 <- function(a,k){
  norm_a <- sqrt(t(a) %*% a)
  #norm_a <- norm(a)
  d <- length(a)
  denominator <- (-1)^k * norm_a ^{2 *k +2} * gamma((d+1)/2) * gamma(k+3/2)
  nominator <- factorial(k) * 2^k * (2*k +1) * (2*k+2) * gamma(k+d/2+1)
  L = denominator / nominator
  return(L)
}

L0(matrix(22,1),10)
L0(100,100)

## -----------------------------------------------------------------------------
# first calculate the left small number
L_k_of_left <- function(a,k){
  norm_a <- sqrt(t(a)%*%a)
  Lkol <- (2*k+2)*log(norm_a) - log(factorial(k)) - k*log(2)-log(2*k+1)-log(2*k+2)
  return(exp(Lkol))
}

#Then calculate the whole part
L_a_k <- function(a,k){
  d <- length(a)
  lkol <- L_k_of_left(a,k)
  lak <- lkol*(-1)^k * beta((d+1)/2,k+3/2) * (k+d/2+1)
  return(lak)
}

L_a_k(matrix(22,1),10)
L_a_k(100,100)

## -----------------------------------------------------------------------------
sum_of_L_a_k <- function(a, echo = FALSE){
  epsilon <- 10^{-8}
  k = 0
  sum_of_lak = 0
  while (abs(L_a_k(a,k)) > epsilon) {
    sum_of_lak <- sum_of_lak + L_a_k(a,k)
    k <- k + 1
  }
  if(echo){
    print(k)
  }
  return(sum_of_lak)
}

sum_of_L_a_k(50)

## -----------------------------------------------------------------------------
a <- matrix(c(1,2),ncol = 1)
sum_of_L_a_k(a)

## -----------------------------------------------------------------------------
# the ck value
c_k <- function(a,k){
  res <- sqrt(a^2 * k / (k+1-a^2))
  return(res)
}

# the parameters in front of the function
param <- function(k){
  res <- 2 * exp(lgamma((k+1)/2) - lgamma(k/2) - log(pi*k)/2)
  return(res)
}

# the integrated function
inte_func <- function(a,k){
  res <- exp(- (k+1)/2 * log(1+a^2/k))
  return(res)
}

# part of the function with k
half_func <- function(a,k){
  res <- param(k) * integrate(inte_func, 0, c_k(a,k), rel.tol=.Machine$double.eps^0.25, k = k)$value
  return(res)
}

solve_func_without_math <- function(k = 4){
  Final <- function(a){
    half_func(a,k) - half_func(a,k-1)
  }
  y <- 1:999/1000
  y <- sqrt(k)*y
  sj <- numeric(length(y))
  for(i in 1:length(y)){
    sj[i] <- Final(y[i])
  }
  for (i in 1:(length(sj)-1)) {
      if(sj[i] * sj[i+1] < 0){
      return(uniroot(Final,c(y[i], y[i+1]))$root )
      }
    
  }
}

solve_func_without_math(k=10)

## -----------------------------------------------------------------------------
P_k <- function(a,k){
  p <- 1 - pt(c_k(a,k),k)
  return(p)
}

S_func <- function(a,k){
  s <- P_k(a,k) - P_k(a,k-1)
  return(s)
}

solve_func <- function(k = 2){
  ssfunc <- function(a){
    ssf <- P_k(a,k) - P_k(a,k-1)
    return(ssf)
  }
  y <- 1:999/1000
  y <- sqrt(k)*y
  sj <- ssfunc(y)
  epsilon <- 10^{-8}
  
  for (i in 1:length(sj)) {
    if(i != length(sj)){
      if(sj[i] * sj[i+1] < 0){
      return(uniroot(ssfunc,c(y[i],y[i+1]))$root )
    }
  }
} 
}

## -----------------------------------------------------------------------------
ks <- c(4:25,100,500,1000)
res1 <- res2 <- numeric(length(ks))
for(i in 1:length(ks)){
  res1[i] <- solve_func(k=ks[i])
  res2[i] <- solve_func_without_math(k=ks[i])
}
matrix(c(ks,res1,res2),ncol = 3)

## -----------------------------------------------------------------------------
Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
# we choose the frequency of Y >= 1 to estimate lambda as lambda0
# which means 0.3 = e^{-\lambda}
lambda0 <- - log(0.3)

N <- 1000
n <- 10
k <- 7
for(j in 1:N){
  lambda1 <- n/(sum(Y) + (n-k)/lambda0)
  if(abs(lambda1 - lambda0) < 10^{-8}){
    iter <- j
    print(c(j,lambda1))
    break
  } else {
    lambda0 <- lambda1
  }
}

lambda2 <- 1/lambda1

## ----eval=FALSE---------------------------------------------------------------
#  trims <- c(0, 0.1, 0.2, 0.5)
#  x <- rcauchy(100)
#  lapply(trims, function(trim) mean(x, trim = trim))
#  lapply(trims, mean, x = x)

## ----eval=FALSE---------------------------------------------------------------
#  rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
trims <- c(1, 2, 3, 4)
x <- rcauchy(100)
lapply(trims, function(trim) sum(x, trim = trim))
lapply(trims, sum, x)
# here the sum function takes trims as the optional parameter in the end

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp), 
  mpg ~ disp+wt,
  mpg ~ I(1/disp)+wt
)
# method 1 for loop
models1 <- list()
for(i in 1:4){
  z <- lm(formula = formulas[[i]],mtcars)
  models1[i] <- list(z)
}

# method 2: lapply
models2 <- lapply(formulas,lm,mtcars)

rsq <- function(mod) summary(mod)$r.squared

# R^2 for method 1
lapply(models1, rsq)
# R^2 for method 2
lapply(models2, rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
         rows <- sample(1:nrow(mtcars), rep = TRUE)
         mtcars[rows, ]
})

# method 1 for loop
models3 <- list()
for(i in 1:10){
  z <- lm(formula = mpg ~ disp, bootstraps[[i]])
  models3[i] <- list(z)
}

# method 2 lapply
models4 <- lapply(bootstraps,lm, formula = mpg ~ disp)

rsq <- function(mod) summary(mod)$r.squared

# R^2 for method 1
lapply(models3, rsq)
# R^2 for method 2
lapply(models4, rsq)

## -----------------------------------------------------------------------------
testdf1 <- data.frame(idx = c(1:12),
                      normrandom = rnorm(12),
                      trandom = rt(12,1), 
                      chisqrandom = rchisq(12,1))
testdf1

vapply(testdf1, sd, FUN.VALUE = numeric(1))

## -----------------------------------------------------------------------------
testdf2 <- data.frame(idx = c(1:12), 
                      normrandom = rnorm(12), 
                      cnidx = c("一","二","三","四","五","六","七","八","九","十","十一","十二"), 
                      trandom = rt(12,1), 
                      chisqrandom = rchisq(12,1))
testdf2
# the class of this dataframe is
cla <- vapply(testdf2, class, FUN.VALUE = character(1))
cla

# extract the numeric data from dataframe
testdf2num <- testdf2[cla == "numeric"]
testdf2num

# compute the sd of every column
vapply(testdf2num, sd, FUN.VALUE = numeric(1))

# you can wirte it in one line code or in a function
vapply(testdf2[vapply(testdf2, class, FUN.VALUE = character(1)) == "numeric"],
       sd, FUN.VALUE = numeric(1))

compute_numeric <- function(df,FUN){
  vapply(df[vapply(df, class, FUN.VALUE = character(1)) == "numeric"],
       sd, FUN.VALUE = numeric(1))
}
compute_numeric(testdf2,sd)

## -----------------------------------------------------------------------------
library(parallel)
# A multicore version of sapply
mcsapply <- function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, 
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), 
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL,
                      simplify = TRUE, USE.NAMES = TRUE) {

    FUN <- match.fun(FUN)
    answer <- mclapply(X, FUN, ..., 
                       mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
                       mc.silent = mc.silent, mc.cores = mc.cores, mc.cleanup = mc.cleanup,
                       mc.allow.recursive = mc.allow.recursive, affinity.list = affinity.list)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!isFALSE(simplify)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

# simple example where the multicore version costs more time
b <- list(x = 1:10, y = matrix(1:12, 3, 4))
system.time(sapply(b, sum))
system.time(mcsapply(b,sum))

# complex example where the multicore version saves more time
jack_df <- function(x){x[-(sample(nrow(x),1)),]}
r2      <- function(model){summary(model)$r.square}
jack_lm <- function(i){r2(lm(mpg ~ I(1/disp)+wt, data = jack_df(mtcars)))}
system.time(sapply(1:5000, jack_lm))
system.time(mcsapply(1:5000, jack_lm))

## ----eval=FALSE---------------------------------------------------------------
#  ###########################
#  sapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
#  {
#      FUN <- match.fun(FUN)
#      answer <- lapply(X = X, FUN = FUN, ...)
#      if (USE.NAMES && is.character(X) && is.null(names(answer)))
#          names(answer) <- X
#      if (!isFALSE(simplify))
#          simplify2array(answer, higher = (simplify == "array"))
#      else answer
#  }
#  
#  vapply <- function (X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE)
#  {
#      FUN <- match.fun(FUN)
#      if (!is.vector(X) || is.object(X))
#          X <- as.list(X)
#      .Internal(vapply(X, FUN, FUN.VALUE, USE.NAMES))
#  }
#  
#  lapply <- function (X, FUN, ...)
#  {
#      FUN <- match.fun(FUN)
#      if (!is.vector(X) || is.object(X))
#          X <- as.list(X)
#      .Internal(lapply(X, FUN))
#  }
#  #########################

## -----------------------------------------------------------------------------
# The following is the R version
cp_gibbsR <- function(length_of_chains = 1e4,from_point = 1001, a = 1, b = 1, x_range = 10, mu_x = NA, mu_y = NA){
  ## length_of_chains: the whole chain length
  ## from_point: from which point to output, notice we actually generate 1e4 + 1000, because we have to dump the first 1000 points
  ## a, b: they are in the beta parameters
  ## x_range: the n parameter in Binomial(n,y)
  ## mu_x.mu_y: the beginning value of this chain
  to_point = length_of_chains + from_point - 1
  Z <- matrix(0,nrow = to_point,ncol = 2)
  # Z initialize
  if (is.na(mu_x) | is.na(mu_y)){ 
    # condition 1: x0 is not defined
    Z[1,2] <- a/(a+b)
    Z[1,1] <- x_range * Z[1,2]
  } else { 
    # condition 2: x0 is defined
    Z[1,] <- c(mu_x,mu_y)
  }
  for (i in 2:to_point){
    # update X
    Zy <- Z[i-1,2]
    Z[i,1] <- rbinom(1,x_range,Zy)
    # update Y
    Zx <- Z[i,1]
    Z[i,2] <- rbeta(1,Zx+a,x_range-Zx+b)
  }
  return(Z[from_point:to_point,])
}

## -----------------------------------------------------------------------------
library(Rcpp)
# the following is the C version
# notice: C use 0, R use 1 at the beginning, so something looks different from R version.
cppFunction('
NumericMatrix cp_gibbsC(int length_of_chains=1e4,int from_point=1000,int a =1,int b = 1,int x_range=10,float mu_x = -1,float mu_y= -1) {
  // x,y can not be negative, and the NA is not accepted, so I choose -1 to avoid no initial values. 
  int to_point = length_of_chains + from_point;
  NumericMatrix Z(to_point, 2);
  // Initial Value
  if (mu_x == -1 || mu_y == -1){
    Z(0,1) = a/(a+b);
    Z(0,0) = x_range * Z(0,1);
  } else {
    Z(0,0) = mu_x;
    Z(0,1) = mu_y;
  }
  for(int i = 1; i < to_point; i++) {
      // Update x
      float Zy = Z(i-1,1);
      Z(i,0) = rbinom(1,x_range,Zy)[0];
      // Update y
      float Zx = Z(i,0);
      Z(i,1) = rbeta(1,Zx+a,x_range-Zx+b)[0];
  }
  return(Z(Range(from_point,to_point-1),_));
  // 1000 - 10999: 10000 in total.
}
')

## -----------------------------------------------------------------------------
ZR <- cp_gibbsR()
ZC <- cp_gibbsC()

## -----------------------------------------------------------------------------
cp_show_chain <- function(x,from_p = 1, to_p = length(x)){
  idx <- from_p:to_p
  y1 <- x[idx]
  plot(idx,y1,type="l",ylab="x")
}

cp_show_gibbs_stat <- function(Z){
  Z_show <- Z
  col_mean_Z <- colMeans(Z_show)
  cov_Z <- cov(Z_show)
  cor_Z <- cor(Z_show)
  print(col_mean_Z)
  print(cov_Z)
  print(cor_Z)
  plot(Z_show,main = "",cex=0.5,xlab="x",ylab="y",ylim=range(Z_show[,2]))
}

## -----------------------------------------------------------------------------
cp_show_gibbs_stat(ZR)
cp_show_gibbs_stat(ZC)

# Result for R version
#par(mfrow=c(2,2))
cp_show_chain(ZR[,1],9500)
cp_show_chain(ZR[,1])
cp_show_chain(ZR[,2],9500)
cp_show_chain(ZR[,2])

# Result for C version
cp_show_chain(ZC[,1],9500)
cp_show_chain(ZC[,1])
cp_show_chain(ZC[,2],9500)
cp_show_chain(ZC[,2])

## -----------------------------------------------------------------------------
qqplot(ZR[,1],ZC[,1],xlab = "ZR",ylab = "ZC")
qqplot(ZR[,2],ZC[,2],xlab = "ZR",ylab = "ZC")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(ZR=cp_gibbsR(),ZC=cp_gibbsC())
summary(ts)[,c(1,3,5,6)]

