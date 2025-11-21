# Robust M-quantile periodogram (`mqper`)

An R function for estimate the so-called M-quantile spectrum of a univariate time series using M-quantile regression.

---

## Description
The `mqper` function calculates a robust periodogram using M-quantile regression, which is less sensitive to outliers compared to traditional Fourier-based methods. This approach is particularly useful for analyzing time series with heavy-tailed noise or contamination.

Key features:
- Robust spectral estimation via M-quantile regression
- Returns spectral estimates and corresponding frequencies
- Works with any quantile level (0 < tau < 1)

---

## Installation
Ensure you have the 'mquantile.R' file
```R
source('mquantile.R')
source("mqper.R")

# Generate sample time series
set.seed(123)
x <- rnorm(256) + 0.5*sin(2*pi*(1:256)/16)  # Sine wave with noise

# Compute median (τ = 0.5) periodogram
result <- mqper(x, tau = 0.5)

# Plot the periodogram
plot(result$freq, result$spec, type = "l", 
     main = "Robust M-Quantile Periodogram",
     xlab = "Frequency", ylab = "Spectral Density")

# Compute periodograms at different quantiles
spec_low <- mqper(x, tau = 0.25)  # Lower quartile
spec_med <- mqper(x, tau = 0.5)   # Median
spec_high <- mqper(x, tau = 0.75) # Upper quartile

# Plot comparison
plot(spec_med$freq, spec_med$spec, type = "l", col = "black",
     main = "M-quantile Periodogram Comparison",
     xlab = "Frequency", ylab = "Spectral Density")
lines(spec_low$freq, spec_low$spec, col = "blue")
lines(spec_high$freq, spec_high$spec, col = "red")
legend("topright", legend = c("τ=0.25", "τ=0.5", "τ=0.75"),
       col = c("blue", "black", "red"), lty = 1)
