setwd("~/Desktop/ICA")

# ICA TASK 1 R SCRIPT:

D <- read.csv("group16.csv", stringsAsFactors=FALSE, header=FALSE)
n <- dim(D)[1]
n

plot(1:n, D[,], xlab="t", ylab=expression('x'[t]))


unscaled_posterior <- function(k) {
  a <- 1+ k/2
  b <- 1+ (n-k)/2
  return (
    (1/(n-1))*((1/((2*pi)^0.5))^n)*
      (gamma(a)/((1+(1/2)*sum(D[1:k,]^2))^a))*
      (gamma(b)/((1+(1/2)*sum(D[(k+1):n,]^2))^b))
  )
}


norm_constant <- sum(sapply(1:(n-1), unscaled_posterior))

true_posterior <- function(k) {
  (unscaled_posterior(k)/norm_constant)
}


plot(sapply(1:(n-1), true_posterior), xlab=expression(paste(tau)), ylab="Posterior density", type='l')
tau <- which.max(sapply(1:(n-1), true_posterior))

tau
true_posterior(tau)

##############################################################################################

# ICA TASK 2 R SCRIPT:

#relevant packages
install.packages("VineCopula")
library(VineCopula)
install.packages("fGarch")
library(fGarch)
install.packages("KScorrect")
library(KScorrect)
install.packages("ADGofTest")
library(ADGofTest)
library(stats)

setwd("~/Desktop/STAT0011/ICA")
CAC40 <- read.csv("CAC40.csv")
DAX <- read.csv("DAX.csv")
#For submission, the above data files were combined into one csv file, called "Task2_CAC40&DAX.csv".

DP <- DAX$Adj.Close #Let DP be the weekly price of DAX.
CP <- CAC40$Adj.Close #Let CP be the weekly price of CAC40.
LDP <- log(DP) #Take the logarithm of the weekly price for DAX.
LCP <- log(CP) #Take the logarithm of the weekly price for CAC40
LRDP <- diff(LDP) #Calculate the weekly log-returns for DAX.
LRCP <- diff(LCP) #Calculate the weekly log-returns for CAC40.
ret1 <- c(0,LRDP) #Because the log-return on the first day is unknown, then replaced by 0.
ret2 <- c(0,LRCP) #Because the log-return on the first day is unknown, then replaced by 0.
date <- DAX$Date
dataset <- cbind(date,ret1,ret2)
View(dataset)
write.table(dataset, file="DesktopICA_Task2.csv", sep=",", row.names = FALSE)
?write.table


jarqueberaTest(ret1) #p Value: < 2.2e-16
jarqueberaTest(ret2) #p Value: < 2.2e-16
#p value is quite small, against the null hypothesis that ret1 and ret2 come from normal distribution.

plot(ret1, type="l")
plot(ret2, type="l")


# Building AR models: the Box - Jenkins approach
# Step 1: Identification
# returns 1
par(mfrow=c(2,2))
acf(ret1, col="green", lwd=2) 
pacf(ret1, col="green", lwd=2) 
acf(ret1^2, col="red", lwd=2) 
par(mfrow=c(1,1))

# returns 2
par(mfrow=c(2,2))
acf(ret2, col="green", lwd=2) 
pacf(ret2, col="green", lwd=2)
acf(ret2^2, col="red", lwd=2) 
par(mfrow=c(1,1))


# Step 2: Estimation
model1=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="sstd")
model2=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="sstd")

#To select conditional distribution for GARCH(2,3) model:
model1b=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="norm")
model1c=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="std")
model1d=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="snorm")
model1e=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="ged")
model1f=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="sged")
model1g=garchFit(formula=~garch(2,3),data=ret1,trace=F,cond.dist="QMLE")

ic=rbind(model1@fit$ics,model1b@fit$ics,model1c@fit$ics,model1d@fit$ics, model1e@fit$ics,model1f@fit$ics,model1g@fit$ics)
ic

#To select conditional distribution for GARCH(2,1) model:
model2b=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="norm")
model2c=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="std")
model2d=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="snorm")
model2e=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="ged")
model2f=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="sged")
model2g=garchFit(formula=~garch(2,1),data=ret2,trace=F,cond.dist="QMLE")

ic2=rbind(model2@fit$ics,model2b@fit$ics,model2c@fit$ics,model2d@fit$ics, model2e@fit$ics,model2f@fit$ics,model2g@fit$ics)
ic2

# Step 3: Model checking
# returns 1
res1 <- residuals(model1, standardize=TRUE)
par(mfrow=c(1,2))
acf(res1, col="green", lwd=2) 
acf(res1^2, col="red", lwd=2) 
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1) ##p-value = 0.9374
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) ##p-value = 0.5461
model1@fit$ics  ## to obtain AIC & BIC
u1<-psstd(res1, mean=0, sd=1, nu=coef(model1)[9], xi=coef(model1)[8]) 
coef(model1)
hist(u1)

# Further distributional checks
#Corrected Kolmogorov-Smirnov test
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value  ##p value: 0.27
#Anderson-Darling test
ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value  ##p value: 0.1116544 


# returns 2
res2 <- residuals(model2, standardize=TRUE)
par(mfrow=c(1,2))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2) 
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1) ##p-value = 0.9736
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) ##p-value = 0.2209
model2@fit$ics 
u2<-psstd(res2, mean=0, sd=1, nu=coef(model2)[7], xi=coef(model2)[6])
coef(model2)
hist(u2)

# Further distributional checks
#Kolmogorov-Smirnov test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value   ##p value: 0.9688
#Anderson-Darling test
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value   ##p value: 0.9888776


# Step 4: Copula modelling
model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05, se = TRUE)
model  ##t (par = 0.91, par2 = 8, tau = 0.73) 

N=10000
u_sim=BiCopSim(N, family=model$family, model$par, model$par2)

res1_sim=qsstd(u_sim[,1], mean = 0, sd = 1, nu=coef(model1)[9], xi=coef(model1)[8]) 
res2_sim=qsstd(u_sim[,2], mean = 0, sd = 1, nu=coef(model2)[7], xi=coef(model2)[6]) 

plot(u_sim)  ##data is not independent
hist(u_sim[,1]) ##looks like standard uniform 


# So, the next step is to re-introduce autocorrelation and GARCH effects observed in data
head(u_sim)

y1<-psstd(ret1, mean=0, sd=1, nu=coef(model1)[9], xi=coef(model1)[8])
y2<-psstd(ret2, mean=0, sd=1, nu=coef(model2)[7], xi=coef(model2)[6])

model=BiCopSelect(y1, y2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05, se = TRUE)
model  ##t (par = 1, par2 = 2.39, tau = 0.99) 

N=10000
y_sim=BiCopSim(N, family=model$family, model$par, model$par2)

y1simulated<-qsstd(y_sim[,1], mean=0, sd=1, nu=coef(model1)[9], xi=coef(model1)[8])
y2simulated<-qsstd(y_sim[,2], mean=0, sd=1, nu=coef(model2)[7], xi=coef(model2)[6])

par(mfrow=c(1,1))
plot(dataset, ylab="ret2", xlab="ret1", col="blue", main="Original log-returns")
plot(data.frame(y1simulated, y2simulated), ylab="ret2_sim", xlab="ret1_sim",
     col="blue", main="Simulated log-returns")
par(mfrow=c(1,1))

portsim <- matrix(0, nrow = N, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

# Equally weighted portfolio, w1=w2=0.5
portsim=log(1+((exp(y1simulated)-1)+(exp(y2simulated)-1))*(1/2))
varsim=quantile(portsim,c(0.01,0.05))
varsim

## 1%        5% 
##-2.827649 -1.758483 



