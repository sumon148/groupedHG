# Application to the study of estimating levels of Salmonella bacteria on pig farms

# Application to Levels of Salmonella bacteria on pig farms: Case 1  --------------


# Define constants
N <- 80  # Population size
barN <- c(20, 15, 10, 5, 1)  # Different group sizes
b <- 1  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta <- function(barN, wp, rho_hat, C) {
  1 -  exp(-wp * C / barN * (1 - exp(-rho_hat / wp)))
}


calculate_delta_alt <- function(barN, wp, rho_hat, C) {
  1 -  (exp(-wp * C * (1 - exp(-rho_hat / wp))))^{1/barN}
}

# Calculate item-level sensitivity (delta)
item_level_sensitivities <- calculate_delta(barN, wp, rho_hat, C)
item_level_sensitivities_alt <- calculate_delta_alt(barN, wp, rho_hat, C)



# Generate FI for varying group sizes : AD and PMF based using HG Distribution ----------

item_level_sensitivities
# 0.1796423 0.2320416 0.3270132 0.5470888 0.9809425

# Different cases: Item level sensitivity

fisher_results_HG_AD <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


fisher_results_HG_PMF <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})






png("Figure/Fig 3 HG AD Case 1.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info,
                  fisher_results_HG_AD[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info,
                  fisher_results_HG_AD[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()






png("Figure/Fig 3 HG PMF Case 1.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info,
                  fisher_results_HG_PMF[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info,
                  fisher_results_HG_PMF[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



range_FI <- c(range_FI[1],0.005)

plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


dev.off()



# Generate FI for varying group sizes : AD and PMF based using BN Distribution ----------

item_level_sensitivities
# 0.1796423 0.2320416 0.3270132 0.5470888 0.9809425

# Different cases: Item level sensitivity

fisher_results_BN_AD <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


fisher_results_BN_PMF <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})






png("Figure/Fig 3 BN AD Case 1.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info,
                  fisher_results_BN_AD[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Binomial",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.0025))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info,
                  fisher_results_BN_AD[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.01, by = 0.001))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info,
                  fisher_results_BN_AD[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.005)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.005, by = 0.0001))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()






png("Figure/Fig 3 BN PMF Case 1.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info,
                  fisher_results_BN_PMF[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.0025))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info,
                  fisher_results_BN_PMF[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.01, by = 0.001))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



range_FI <- c(range_FI[1],0.005)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.005, by = 0.0005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


dev.off()




# Application to Levels of Salmonella bacteria on pig farms: Case 2  --------------

# Define constants
N <- 80  # Population size
barN <- c(20, 10, 5, 1)  # Different group sizes
b <- c(1,2,4,20)  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
# wn <- b*ws*barN <- 20  # Total weight of the pooled sample
# wp <- ws*barN <- 20  # Total weight of the pooled sample
ws <- 1  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta_ws <- function(barN, ws, rho_hat, C) {
  1 -  exp(-ws * C * (1 - exp(-rho_hat / (ws*barN))))
}


# Calculate item-level sensitivity (delta)
item_level_sensitivities_case2 <- calculate_delta_ws(barN, ws, rho_hat, C)

# Generate FI for varying group sizes : AD and PMF based using HG Distribution -------------

item_level_sensitivities_case2
# 0.1796423 0.3233888 0.5325257 0.9544191

fisher_results_HG_AD <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})



fisher_results_HG_PMF <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


png("Figure/Fig 3 HG AD Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info)

range_FI <- c(0,1)
plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()



png("Figure/Fig 3 HG PMF Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info)

#range_FI <- c(0,1)
plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



# Plot results
range_FI <- c(range_FI[1],0.03)

plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.03, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()





# ----------------------------------------------------#
# Generate FI for varying group sizes : AD and PMF based using BN Distribution -------------

item_level_sensitivities_case2
# 0.1796423 0.3233888 0.5325257 0.9544191

fisher_results_BN_AD <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})



fisher_results_BN_PMF <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


png("Figure/Fig 3 BN AD Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.03)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.03, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.005)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.005, by = 0.0005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

dev.off()



png("Figure/Fig 3 BN PMF Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info)

#range_FI <- c(0,1)
plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



# Plot results
range_FI <- c(range_FI[1],0.03)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.03, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.0005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


dev.off()





