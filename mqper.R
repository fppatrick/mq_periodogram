#' Robust M-Quantile Periodogram
#'
#' Computes the robust M-quantile periodogram from univariate time series using
#' the M-quantile regression as introduced in Breckling & Chambers (1988).
#' This estimator provides robustness to outliers and heavy-tailed noise. Furthermore,
#' it offer a more rich view of the spectrum than the classical approach.
#'
#' @details
#' For each Fourier frequency \eqn{\lambda_j = 2 \pi j / n}, the function fits
#' the M-quantile regression:
#' \deqn{ X_t = \beta_{1,j}(\tau) \cos(\lambda_j t) +
#'        \beta_{2,j}(\tau) \sin(\lambda_j t) + \varepsilon_t, }
#'
#' using the `mqlm` estimator (\url{https://github.com/fppatrick/mq_regression}).
#'  
#' The robust M-quantile periodogram is then defined as
#' \deqn{
#'     I_\tau(\lambda_j)
#'       = \frac{n}{4} \left( \beta_{1,j}^2(\tau) + \beta_{2,j}^2(\tau) \right).
#' }
#'
#' @param timeseries A numeric vector representing a univariate time series.
#' @param tau Quantile level used in the M-quantile regression.
#'   Must satisfy \eqn{0 < \tau < 1}.
#'
#' @return A list with:
#' \item{perior}{Vector of robust spectral estimates.}
#' \item{freq}{Vector of corresponding Fourier frequencies.}
#'
#' @author Patrick Ferreira Patrocinio
#'
#' @import mquantile.R
#' @export

mqper <- function(timeseries, tau) {

  if (!is.numeric(timeseries))
    stop("'timeseries' must be numeric.")

  n <- length(timeseries)
  if (n < 4)
    stop("Time series must have length >= 4.")

  if (!is.numeric(tau) || tau <= 0 || tau >= 1)
    stop("'tau' must satisfy 0 < tau < 1.")

  g <- n %/% 2
  perior <- numeric(g)
  FFT    <- complex(length.out = g)

  t_index <- seq_len(n)

  for (j in 1:g) {

    w <- 2 * pi * j / n
    X1 <- cos(w * t_index)

    if (j == n / 2 && n %% 2 == 0) {
      MX <- cbind(X1)
      fit <- mqlm(MX, timeseries, maxit = 100, q = tau, k = 1.345)
      FFT[j] <- sqrt(n / 4) * complex(real = fit$coef[2], imaginary = 0)
    } else {
      X2 <- sin(w * t_index)
      MX <- cbind(X1, X2)
      fit <- mqlm(MX, timeseries, maxit = 100, q = tau, k = 1.345)
      FFT[j] <- sqrt(n / 4) * complex(real = fit$coef[2],
                                      imaginary = -fit$coef[3])
    }

    perior[j] <- Mod(FFT[j])^2
  }

  freq <- 2 * pi * seq_len(g) / n

  if (n %% 2 == 0) {
    return(list(perior = perior, freq = freq))
  } else {
    return(list(perior = perior[-g], freq = freq[-g]))
    
  }
}
