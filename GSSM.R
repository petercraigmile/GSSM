
## ======================================================================
## Copyright 2013--2014, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## File     : GSSM.R
## Contains : R code to simulate stationary Gaussian processes
##            using the Gaussian Spectral Synthesis Method (GSSM).
## Updated  : pfc@stat.osu.edu, April 2014.
##
## Reference: pp. 291-292, D. B. Percival and A. T. Walden. Wavelet
##            Methods for Time Series Analysis. Cambridge University
##            Press, Cambridge, 2000 (P&W).
## ======================================================================



GSSM <- function (N, the.sdf, L=4, LM=FALSE, sum.acvs, ...) {
  ## ======================================================================
  ## Using the Gaussian Spectral Synthesis Method (GSSM), approximately
  ## simulate a mean zero stationary Gaussian process of length 'N*L',
  ## truncated to length 'N', with a given spectral density function
  ## 'the.sdf', with arguments '...'.
  ##
  ## The approximation improves as L increases.
  ##
  ## For long memory processes, you need to use the adjustment that
  ## appears at the bottom of p.291 and top of p.292 of P&W.  To
  ## approximately simulate a long memory process you need to provide
  ## a 'sum.acvs' value which is equal to
  ##
  ##  \sum_{\tau = -(N-1)}^{N-1} s_{X, \tau}
  ##  = 2 \int_0^{1/2} \frac{\sin([2N_1] \pi f)}{\sin (\pi f)} S_X(f) df,
  ##
  ## where s_{X,\tau} is the autocovariance sequence (ACVS) at lag
  ## \tau and S_X(f) is the spectral density function.  P&W warn that
  ## the integral calculation (if no ACVS is available) needs to be
  ## numerically evaluated 'with sufficient care'.
  ##
  ## Peter Craigmile, pfc@stat.osu.edu, Apr 2014
  ## ======================================================================
  
  M <- N*L
  K <- M/2

  if (missing(sum.acvs)) {

    ss <- the.sdf((0:K)/M, ...)

  } else {

    wk <- pi * (1:(M-1))/M

    ss <- c(NA, the.sdf((1:K)/M, ...))
    
    ss[1] <- M * (sum.acvs -
                  sum(ss[c(2:(K+1), K:2)] * sin((2*N-1) * wk) / sin(wk))/M) / (2*N-1)
  }

  w <- sqrt(0.5*ss[2:K]) * (rnorm(K-1) + 1i * rnorm(K-1))

  U <- c(sqrt(ss[1]) * rnorm(1),
         w,
         sqrt(ss[K+1]) * rnorm(1),
         rev(Conj(w)))
  
  Re(fft(U, inverse=TRUE))[1:N]/sqrt(M)
}


