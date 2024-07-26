# Purpose: To compare expectation and variance of across and within approach
# Creator: Matthew LH. Cheng
library(here)
library(tidyverse)

# Extract a vector of CAA from Experiment 3
om_path <- here("output", "Experiment 3")
om_scenario <- here(om_path, "Fem60_Mal40")
load(here(om_scenario, paste("Fem60_Mal40",".RData",sep = "")))
list2env(oms,globalenv()) # output into global environment
caa <- oms$CAA[30,,,,1] # Get CAA

n_sims <- 1e3 # number of simulations
f <- caa[,1] # females
m <- caa[,2] # males
iss <- 1e5 # input sample size

# par(mfrow = c(4,2))
# Simulation across -------------------------------------------------------
across <- c(f,m) / sum(f,m) # Compute Across proportions
across_sim <- rmultinom(n_sims, iss, across) # Simulate values for comparison
cor_across_sim <- cov2cor(cov(t(across_sim))) # get correlation matrix

# compare expected mean
plot(rowMeans(across_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Expectation", main = "Proportions Across (Mean)")
lines(iss * across, type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 10000, legend = c("Simulated", "Analytical"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# compare expected variance here
across_var_theor <- iss * across * ( 1 - across) # calculate theoretical
plot(matrixStats::rowVars(across_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Variance", main = "Proportions Across (Variance)")
lines(across_var_theor, type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 8000, legend = c("Simulated", "Analytical"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# Simulation Within -------------------------------------------------------
f_w <- f / sum(f) # Female proportions
m_w <- m / sum(m) # Male proportions
rat_f <- sum(f) / sum(f,m) # Female ratio
rat_m <- 1 - rat_f # Male ratio

# binomial
samps_f <- rbinom(n_sims, iss, rat_f) # female samples
samps_m <- iss - samps_f # male samples

# multinomial
within_f_sim <- rmultinom(n_sims, samps_f, f_w) # female multinomial samples
within_m_sim <- rmultinom(n_sims, samps_m, m_w) # male multinomial samples
within_sim <- rbind(within_f_sim, within_m_sim) # concatenate the two vectors
cor_within_sim <- cbind(cov2cor(cov(t(within_f_sim))), cov2cor(cov(t(within_m_sim)))) # get correlation matrix

# compare expected mean
plot(rowMeans(within_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Expectation", main = "Proportions Within (Mean)")
lines(c(f_w * iss * rat_f, m_w * iss * rat_m), type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 10000, legend = c("Simulated", "Analytical"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# compare expected variance here
within_var_theor <- c((iss * rat_f) * f_w * (1 - f_w), # females
                      (iss * rat_m) * m_w * (1 - m_w)) # males
plot(matrixStats::rowVars(within_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Variance", main = "Proportions Within (Variance)")
lines(within_var_theor, type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 7000, legend = c("Simulated", "Analytical"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# Compare empirical -------------------------------------------------------

# Empirical Means
plot(rowMeans(across_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Expectation", main = "Empirical Means")
lines(rowMeans(within_sim), type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 10000, legend = c("Across", "Within"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# Empirical Variance
plot(matrixStats::rowVars(across_sim), type = 'l', col = '#35577b', lwd = 3.5, xlab = "Index", ylab = "Variance", main = "Empirical Variances")
lines(matrixStats::rowVars(within_sim), type = 'l', lty = 2, lwd = 3.5, col = '#e69c4c')
legend(x = 45, y = 8000, legend = c("Across", "Within"), fil = c("#35577b", "#e69c4c"))
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.2))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.2))

# Compare analytical -----------------------------------------------------

# Analytical Means
png(here("figs", "ms_figs", "Acr_With_Comparison.png"), width = 850, height = 700)
par(mfrow = c(2,2), mar = c(4, 6, 1, 6))
plot(iss * across, type = 'l', col = '#35577b', lwd = 4, xlab = "Index", ylab = "Expectation", cex.lab = 1.5, cex.axis = 1.5)
lines(c(f_w * iss * rat_f, m_w * iss * rat_m), type = 'l', lty = 2, lwd = 4, col = '#e69c4c')
legend(x = 40, y = 10000, legend = c("Across", "Within"), fil = c("#35577b", "#e69c4c"), cex = 1.5)
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.075))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.075))

# Analytical Variance
plot(across_var_theor, type = 'l', col = '#35577b', lwd = 4, xlab = "Index", ylab = "Variance", cex.lab = 1.5, cex.axis = 1.5)
lines(within_var_theor, type = 'l', lty = 2, lwd = 4, col = '#e69c4c')
legend(x = 40, y = 8000, legend = c("Across", "Within"), fil = c("#35577b", "#e69c4c"), cex = 1.5)
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.075))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.075))

# Differences
plot(round(iss * across - c(f_w * iss * rat_f, m_w * iss * rat_m)), cex.lab = 1.5, cex.axis = 1.5,
     xlab = "Index", ylab = "Difference in Expectation (Across - Within)", lty = 1, lwd = 4, type = 'l')
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.075))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.075))

plot(across_var_theor - within_var_theor, cex.lab = 1.5, cex.axis = 1.5,
     xlab = "Index", ylab = "Difference in Variance (Across - Within)", lty = 1, lwd = 4, type = 'l')
rect(xleft = -3, xright = 30, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#DC3220", alpha = 0.075))
rect(xleft = 30, xright = 63, ybottom = -1e3, ytop = 1e5, col = adjustcolor("#005AB5", alpha = 0.075))
dev.off()

