# Load the mgcv library
library(mgcv)

# Set a seed for reproducibility
set.seed(42)

# 1. Generate a more extreme dataset (smaller and noisier)
n <- 20 # Reduced sample size
noise_level <- 0.8 # Increased noise
x <- runif(n, 0, 1)
true_f <- sin(2 * pi * x)
y <- true_f + rnorm(n, 0, noise_level)
dat <- data.frame(x = x, y = y)

# 2. Fit a GAM to the data
b <- gam(y ~ s(x, k = 8, bs = "cr"), data = dat, method = "REML")

# Create a sequence of new data points for a smooth plot
new_data <- data.frame(x = seq(0, 1, length.out = 200))

# 3. Calculate Uncorrected Pointwise CIs
pred_uncorrected <- predict(b, newdata = new_data, se.fit = TRUE)
uncorrected_upper <- pred_uncorrected$fit + 1.96 * pred_uncorrected$se.fit
uncorrected_lower <- pred_uncorrected$fit - 1.96 * pred_uncorrected$se.fit

# 4. Calculate Corrected Pointwise CIs
Xp <- predict(b, newdata = new_data, type = "lpmatrix")
V_corrected <- vcov(b, unconditional = TRUE)
corrected_se <- sqrt(diag(Xp %*% V_corrected %*% t(Xp)))
corrected_upper <- pred_uncorrected$fit + 1.96 * corrected_se
corrected_lower <- pred_uncorrected$fit - 1.96 * corrected_se

# 5. Calculate Simultaneous CIs via Simulation
set.seed(42) # for reproducibility of the simulation
N_sim <- 10000
sim_coefs <- rmvn(N_sim, coef(b), V_corrected)
sim_fits <- Xp %*% t(sim_coefs)

# Find the critical distance for the simultaneous interval
# --- THIS IS THE CORRECTED LINE ---
abs_deviations <- abs(sweep(sim_fits, 1, pred_uncorrected$fit, "-"))
# --- END OF CORRECTION ---

max_abs_deviations <- apply(abs_deviations, 2, max) # Max deviation for each simulation
critical_distance <- quantile(max_abs_deviations, 0.95) # 95th percentile of max deviations

# Construct the simultaneous interval
simultaneous_upper <- pred_uncorrected$fit + critical_distance
simultaneous_lower <- pred_uncorrected$fit - critical_distance

# 6. Plot all results for comparison
plot(y ~ x, data = dat, pch = 19,
     main = "Pointwise vs. Simultaneous Confidence Intervals",
     xlab = "x", ylab = "y", ylim = c(-5,5))

# Add the confidence intervals as shaded polygons, from widest to narrowest
# Simultaneous CI (blue)
polygon(c(rev(new_data$x), new_data$x), c(rev(simultaneous_upper), simultaneous_lower),
        col = rgb(0, 0, 1, 0.2), border = NA)
# Corrected Pointwise CI (purple)
polygon(c(rev(new_data$x), new_data$x), c(rev(corrected_upper), corrected_lower),
        col = rgb(0.5, 0, 0.5, 0.2), border = NA)
# Uncorrected Pointwise CI (red)
polygon(c(rev(new_data$x), new_data$x), c(rev(uncorrected_upper), uncorrected_lower),
        col = rgb(1, 0, 0, 0.2), border = NA)

# Add the fitted line on top
lines(new_data$x, pred_uncorrected$fit, col = "black", lwd = 2)

# Add legend
legend("topright",
       legend = c("Fitted Smooth", "Simultaneous 95% CI", "Corrected Pointwise CI", "Uncorrected Pointwise CI"),
       col = c("black", "blue", "purple", "red"),
       lty = c(1, NA, NA, NA),
       fill = c(NA, rgb(0, 0, 1, 0.2), rgb(0.5, 0, 0.5, 0.2), rgb(1, 0, 0, 0.2)),
       border = "white",
       bty = "n")

