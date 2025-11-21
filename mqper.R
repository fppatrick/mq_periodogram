#' Robust M-quantile periodogram
#' This function computes the univariate robust M-quantile periodogram using
#' M-quantile regression.
#' @param series univariate time series
#' @param tau quantile level (0 < tau < 1)
#' @return spec: Vector of robust level-crossings spectral estimates
#' @return freq: Vector of corresponding frequencies
#' @author Patrick Ferreira Patrocinio
#' @import MASS
library(MASS)
mqper <- function(series, tau) {
  n <- length(series)
  g <- n %/% 2
  per <- FFT <- rep(0,g)
  for (j in 1:(g-1)) {
    X1 <- X2 <- NULL
    w <- 2 * pi * j / n
    X1 <- cos(w * 1:n)
    X2 <- sin(w * 1:n)
    MX <- cbind(X1, X2)
    fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
    FFT[j] <- sqrt(n / (8 * pi)) * complex(real = fitrob$coefficients[2], imaginary = -fitrob$coefficients[3])
    
  }
  w <- 2 * pi * seq_len(g) / n
  X1 <- cos(w * 1:n)
  fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
  FFT[g] <- sqrt(n / (8 * pi)) * complex(real = fitrob$coefficients[2], imaginary = -0)
  per <- Mod(FFT)^2
  
  perior <- c(per[-1], rev(per))
  
  w <- 2 * pi * seq.int(g) / n
  return(list(spec = perior, freq = w))
}

