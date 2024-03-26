

#Packages
library(fitdistrplus)
library(eeptools)
library(actuar)
library(pracma)
library(Rlab)
library(QRM)

#Number of simulations.
R <- 150000

#Vector to store summed variables S_{n-1} = X_1 + X_2 +...+ X_n-1.
vectS<-c(0)

#Size of portfolio, 5960 loans.
PortSize <- length(hmeq$LOAN)

#Exposure at default taken as 1.
EAD <- rep(1,PortSize)

#for loop to generate vector of S_{n-1} = X_1 + X_2 + ... + X_n-1.
for (i in 1:R) {
  
  #Default probability vector to implement Bernoulli mixture model with a beta mixing distribution.
  PD<- rbinomial.mixture(PortSize, m = 1, model = c("beta"), shape1 = 0.05638, shape2 = 0.22566)
  
  #Loss given default modeled with beta(0.2,0.95) distribution.
  LGD <- rbeta(length(hmeq$LOAN),0.2,0.95)
  
  #Loss portfolio.
  LossVariable <- EAD * LGD * PD
  
  #Vector to store sum of n-1 losses.
  vectS[i] <- sum(LossVariable[1:length(LossVariable)-1])
}

#Sample used to fit X1,X2,...,Xn to a distribution.
PD<- rbinomial.mixture(PortSize, m = 1, model = c("beta"), shape1 = 0.05638, shape2 = 0.22566)
LGD <- rbeta(length(hmeq$LOAN),0.2,0.95)
EAD <- rep(1,PortSize)
sample <- PD*EAD*LGD
summary(sample)
hist(sample)

#Scaling may be necessary to better fit distributions. Take values greater than 0.
adjSAMPLE <- sample[sample>0.0]

#Plot histogram of adjusted sample and the distribution of losses.
par(mfrow=c(1,2))
hist(adjSAMPLE,main = c("Adjusted Sample"), xlab = "Losses", breaks = 50)
hist(vectS, main = c("Loss Distribution"), xlab = "Losses", breaks = 50)
par(mfrow=c(1,1))

#Fit adjSAMPLE to four distributions.
fit_beta <- fitdist(adjSAMPLE,"beta", start = list(shape1 = 1, shape2 = 1))
fit_lognorm <- fitdist(adjSAMPLE, distr = "lnorm")
fit_weibull <- fitdist(adjSAMPLE,distr = "weibull", start = list(shape = 0.05, scale = 0.1))
fit_gamma <- fitdist(adjSAMPLE, "gamma", start = list(shape = 1, rate = 20))

#Summaries of the four fits.
summary(fit_beta)
summary(fit_weibull)
summary(fit_lognorm)
summary(fit_gamma)

#Goodness-of-fit tests.
gofstat(list(fit_weibull,fit_lognorm,fit_gamma,fit_beta))

#Might need to scale vectS for evaluation purposes and scaled by 1000 to plot.
summary(vectS)
ScaledTotalLosses <- vectS

#This is the conditional MC estimate for each S_(n-1) in ScaledTotalLosses.
Fncond <- function(x,ScaledTotalLosses) {
  s<-c(0) 
  for (i in 1:R) {
    s[i] <- pbeta(x-ScaledTotalLosses[i], fit_beta$estimate[1],fit_beta$estimate[2])
  }
  return(mean(s))
}

#Quantile.
alpha <- 0.95

#Solve Fncond(1/(1-q)-1,ScaledTotalLosses) = 0.95 for q, store as qu.
qu <- uniroot(f = function(q) alpha - Fncond(1/ (1 - q) - 1,ScaledTotalLosses), interval = c(0,1))$root

#VaR estimate.
VaR <- 1 /(1 - qu) - 1

#Multiply back by 10
actualVaR <- 10*VaR

#Estimate PDF using conditional MC.
fncond <- function(x,ScaledTotalLosses) {
  s<-c(0)
  for (i in 1:R) {
    s[i] <- dbeta(x-ScaledTotalLosses[i],fit_beta$estimate[1],fit_beta$estimate[2]) #this s will hold all evaluations from F(x-S_n)
  }
  return(mean(s))
}

###### Expected Shortfall ######

#Compute inverse of the estimated Fncond, using the same ScaledTotalLosses as above.
Fninvscond <- function(x) {
  sINV<-c(0) 
  for (i in 1:R) {
    sINV[i] <- pbeta(x-ScaledTotalLosses[i],fit_beta$estimate[1],fit_beta$estimate[2])
  }
  return(1-mean(sINV))
}

#Vectorize the inverse distribution function.
FninvscondVect <- Vectorize(Fninvscond, "x")

#Integrate the inverse distribution function from the VaR estimate to infinity to get the estimate for m.
m <- integral(FninvscondVect,VaR,Inf)

#Compute expected shortfall. ES should be greater than or equal to VaR.
ES <- VaR + (m)/(1-alpha)
actualES <- 10*ES

#Density function plot.
fncondVect <- Vectorize(fncond,"x")
p <- seq(0.15,0.37,by = 0.0053)
plot(p,fncondVect(p,ScaledTotalLosses),type = "l",xlab = "Losses",ylab = "Density",main = "Loss Density")
legend("topright",leg = paste0(c("99.99%-VaR","99.99%-ES")),lty = 1, col = c("darkgreen","red"))
abline(v=.24601,col = "darkgreen")
abline(v = 0.24982, col = "red")
