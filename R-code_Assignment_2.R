load("~/Desktop/Edinburgh/Semester 1/Incomplete Data Analysis/Assignment 2/dataex2.Rdata")
load("~/Desktop/Edinburgh/Semester 1/Incomplete Data Analysis/Assignment 2/dataex4.Rdata")
load("~/Desktop/Edinburgh/Semester 1/Incomplete Data Analysis/Assignment 2/dataex5.Rdata")

########## Question 2 ##########
############ Task b ############

#normally distributed with sigma^2 known and mu unknown
library(maxLik)

x <- dataex2$X

logLikdataex2 <- function(param){
  mu <- param[1]
  sigma = 1.5
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}
mle <- maxLik(logLik = logLikdataex2, start = c(mu = 5))
summary(mle)

########## Question 4 ##########

x_all= dataex4$X
x_mis = dataex4[which(is.na(dataex4$Y)),]$X
y_obs = dataex4[which(is.na(dataex4$Y) == FALSE),]$Y
x_obs = dataex4[which(is.na(dataex4$Y)== FALSE),]$X

log_like_dataex4 <- function(param, data){
  x <- data[,1]; y <- data[,2]
  beta0 <- param[1]
  beta1 <- param[2]
  sum(y_obs*(beta0+x_obs*beta1))-(beta0+sum(x_mis*beta1))^2-log(1+exp(beta0+sum(x_all*beta1)))
}
mle <- maxLik(logLik = log_like_dataex4, data = dataex4, start = c(beta0 = 1, beta1 = 1))
summary(mle)



####tried another way and this did not work
expex <- function(dataex4, Beta, eps){
  diff <- 1
  beta <- Beta
  beta0 <- Beta[1]
  beta1 <- Beta[2]
  while(diff > eps){
    beta.old <- beta
    beta0 <- sum(y_obs)-2*(beta0+sum(x_mis*beta1))#-((exp(beta0+sum(x_all*beta1)))/(1+exp(beta0+sum(x_all*beta1))))
    beta1 <- sum(y_obs*x_obs)-2*sum(x_mis*(beta0+sum(x_mis*beta1)))#-((sum(x_all*exp(beta0+sum(x_all*beta1))))/(1+exp(beta0+sum(x_all*beta1))))
    beta <- c(beta0, beta1)
    print(beta0)
    print(beta1)
    diff <- sum(abs(beta-beta.old))
  }
  return(theta)
}
resultat <- expex(dataex4, Beta = c(0, 0),eps=  0.00001)
beta0 <- resultat[1]; beta1 <- resultat[2];


########## Question 5 ##########
############ Task b ############

dataex5 <- data.frame(dataex5)
colnames(dataex5) = 'y'
em.mixture <- function(y, theta0, eps){
  n <- length(y)
  theta <- theta0
  p <- theta[1]; lambda <- theta[2]; mu <- theta[3]; sigma2 <- theta[4]
  diff <- 1
  while(diff > eps){
    theta.old <- theta
    #E-step
    ptilde1 <- p*dlnorm(y, meanlog = mu, sdlog = sigma2)
    ptilde2 <- (1-p)*dexp(y, rate = lambda)
    ptilde <- ptilde1/(ptilde1 + ptilde2)
    #M-step
    p <- (1/n)*sum(ptilde)
    lambda <- sum(1-ptilde)/sum(y*(1-ptilde))
    mu <- sum(ptilde*(log(y)))/sum(ptilde)
    sigma2 <- 3*sum(ptilde*(log(y)-mu)**2)/sum(ptilde)
    theta <- c(p, lambda,mu, sigma2)
    diff <- sum(abs(theta-theta.old))
    print(diff)
  }
  return(theta)
}

res <- em.mixture(y = dataex5$y, theta0 = c(0.1, 0.15, 1, 20), eps = 0.00001)
pest <- res[1]; lambdaest <- res[2]; muest <- res[3]; sigma2est <- res[4]
pest; lambdaest; muest; sigma2est

x = dataex5$y
hist(dataex5$y, main = "Estimated Density",
     xlab = "x", freq = FALSE, 
     ylim =c(0, 0.14))
curve(pest*dlnorm(x, meanlog = muest, sdlog = sigma2est) + (1-pest)*dexp(x, rate = lambdaest),
      add = TRUE, col = "blue2")
