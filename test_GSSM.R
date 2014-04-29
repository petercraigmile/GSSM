
## ======================================================================
## Some test code for GSSM
## pfc@stat.osu.edu, April 2014.
## ======================================================================

source("GSSM.R")



## An autoregressive (AR) process

ar.sdf <- function (freqs, phi, sigma2 = 1, delta.t = 1) 
{
    ws <- -2 * pi * delta.t * freqs
    js <- seq(length(phi))
    reals <- sapply(ws, function(w, js, phi) 1 - sum(phi * cos(js * 
        w)), js = js, phi = phi)
    imags <- sapply(ws, function(w, js, phi) sum(phi * sin(js * 
        w)), js = js, phi = phi)
    (sigma2 * delta.t)/(reals * reals + imags * imags)
}

x <- GSSM(1024, ar.sdf, phi=0.6)

par(mfrow=c(2,3), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")

plot(x, ylab="AR(1) process", type="l")
acf(x)
pacf(x)



## A fractionally differenced (FD) process

fd.sdf <- function (fs, d, sigma2=1, deltat=1, deriv=0)
{
  tsy  <- abs(2 * sin(pi * fs * deltat))
  deltat * sigma2 * tsy^(-2*d) * (-2*log(tsy))^deriv 
}

fd.acvs <- function (max.lag, d, sigma2=1)
{
  acvs0 <- sigma2 * exp(lgamma(1-2*d)-2*lgamma(1-d))
  if (max.lag>0)
  {
    ks <- 1:max.lag
    cumprod(c(acvs0, (ks-1+d)/(ks-d)))
  }
  else acvs0
}


N <- 1024
d <- 0.4

the.acvs <- fd.acvs(N-1, d=d)

sum.acvs <- the.acvs[1] + 2*sum(the.acvs[-1])

y <- GSSM(N, fd.sdf, d=d, LM=TRUE, sum.acvs=sum.acvs)

plot(y, type="l", ylab="FD(d) process")
acf(y)
pacf(y)

