# Load necessary libraries
library(coda)
library(lattice)
set.seed(42)
# Estimating Poisson Mean using ABC

# Set true lambda and sample size
trueLambda <- 2
sampleSize <- 5

# Example data
dataY <- c(1, 1, 2, 2, 3)
priorAlpha <- 1
priorBeta <- 1

numSimulations <- 100000
acceptedLambdas <- numeric(0)
tolerance <- 1  # Define a tolerance level for accepting lambda

for (i in 1:numSimulations) {
  sampledLambda <- rgamma(1, shape = priorAlpha, rate = priorBeta)
  simulatedData <- sort(rpois(sampleSize, lambda = sampledLambda))
  if (sum((simulatedData - dataY)^2) <= tolerance) {  # Less strict condition
    acceptedLambdas <- c(acceptedLambdas, sampledLambda)
  }
}

# Plotting the results
histogramObject <- hist(acceptedLambdas, 20, plot = FALSE)
plot(histogramObject$mids, histogramObject$density, ylim = c(0, 1.1), type = 'l', col = 2, xlab = 'Lambda', ylab = 'Posterior', lwd = 2)
lines(histogramObject$mids, dgamma(histogramObject$mids, shape = priorAlpha + sum(dataY), rate = priorBeta + sampleSize), lwd = 2)

# ABC with mean as statistic
numABC <- 10000
lambdaABC <- meanABC <- numeric(0)
meanDataY <- mean(dataY)
toleranceDelta <- 0.5
for (i in 1:numABC) {
  sampledLambda <- rgamma(1, shape = priorAlpha, rate = priorBeta)
  simulatedData <- sort(rpois(sampleSize, lambda = sampledLambda))
  if (abs(mean(simulatedData) - meanDataY) < toleranceDelta) {
    lambdaABC <- c(lambdaABC, sampledLambda)
    meanABC <- c(meanABC, mean(simulatedData))
  }
}
histogramABC <- hist(lambdaABC, 50, plot = FALSE)
lines(histogramABC$mids, histogramABC$density, type = 'l', col = 3, lwd = 2)

# Regression correction for ABC
adjustLambda <- function(lambda, simMeans, dataMean) {
  regressionCoef <- coef(lm(lambda ~ simMeans))[2]
  adjustedLambda <- lambda + regressionCoef * (dataMean - simMeans)
  return(adjustedLambda)
}

lambdaAdjusted <- adjustLambda(lambdaABC, meanABC, meanDataY)
histogramAdjusted <- hist(lambdaAdjusted, 50, plot = FALSE)
lines(histogramAdjusted$mids, histogramAdjusted$density, type = 'l', col = 3, lwd = 2, lty = 2)
legend("topright", 
       legend = c("Rejection Sampling", "True Posterior", "ABC Result", "ABC Adjusted"),
       col = c(2, 1, 3, 3),
       lwd = 2,
       lty = c(1, 1, 1, 2))
 