
library(ggpubr)
library(ggplot2)
library(metR)

# Application to the study of estimating levels of Salmonella bacteria on pig farms

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta_ws <- function(barN, ws, rho_hat, C) {
  1 -  exp(-ws * C * (1 - exp(-rho_hat / (ws*barN))))
}

N <- 60 # Population size
max_group_size <- 20
max_groups <- 20
wn <- 20         # a total weight of faecal matter tested
lambda <- 1      # Specificity for item-level
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- c(9,15,51)  # Range of contaminated items



# Combination of b and barN based on the assumptions

ws <- c()
b_store <- c()
barN_store <- c()

# Loop through possible combinations
for (b in 1:max_groups) {

  for (N_bar in 1:max_group_size) {

    if (b * N_bar <= N) {
      # Store group size and number of groups
      b_store <- c(b_store, b)
      barN_store <- c(barN_store, N_bar)

      # Calculate ws and store
      ws <- c(ws, wn / (b * N_bar))
    }

  }
}

# Combine into a data frame
combinations <- data.frame(b = b_store, barN = barN_store, ws = ws)


# Calculate item-level sensitivity (delta)
combinations$delta <- calculate_delta_ws(combinations$barN, ws=combinations$ws, rho_hat, C)

# Calculation of Hellinger Information

fisher_results_HG_PMF <- lapply(1:dim(combinations)[1], function(i) {

  # Extract corresponding values
  N <- 60
  n <- combinations$barN[i]
  b_i <- combinations$b[i]  # Ensure correct mapping of b to barN
  ws_i <- combinations$ws[i]
  # Calculate delta based on barN
  delta <- combinations$delta[i]
  lambda <- 1
  # Compute Fisher Information for each Tx value
  fisher_info <- lapply(Tx_values, function(Tx) {
    info_hg_tx_imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF-HI", type = "item",verbose=FALSE
    )
  })
  HI <- sapply(fisher_info, function(x) x$FI_Tx)

  # Store results for this barN
  list(barN = n, b=b_i,ws=ws_i,delta=delta,lambda=lambda,HI = HI)
})

HI_HG_Tx_9 <- rep(NA, dim(combinations)[1])
HI_HG_Tx_15 <- rep(NA, dim(combinations)[1])
HI_HG_Tx_51 <- rep(NA, dim(combinations)[1])
for (i in 1: dim(combinations)[1]) {
  HI_HG_Tx_9[i] <- fisher_results_HG_PMF[[i]]$HI[1]
  HI_HG_Tx_15[i] <- fisher_results_HG_PMF[[i]]$HI[2]
  HI_HG_Tx_51[i] <- fisher_results_HG_PMF[[i]]$HI[3]
}

combinations$HI_HG_Tx_9 <- HI_HG_Tx_9
combinations$HI_HG_Tx_15 <- HI_HG_Tx_15
combinations$HI_HG_Tx_51 <- HI_HG_Tx_51



combinations_HG_Tx_9 <- combinations[,c("b","barN","HI_HG_Tx_9")]
combinations_HG_Tx_15 <- combinations[,c("b","barN","HI_HG_Tx_15")]
combinations_HG_Tx_51 <- combinations[,c("b","barN","HI_HG_Tx_51")]



fisher_results_BN_PMF <- lapply(1:dim(combinations)[1], function(i) {

  # Extract corresponding values
  N <- 60
  n <- combinations$barN[i]
  b_i <- combinations$b[i]  # Ensure correct mapping of b to barN
  ws_i <- combinations$ws[i]
  # Calculate delta based on barN
  delta <- combinations$delta[i]
  lambda <- 1
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    info_bn_tx_imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF-HI", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, b=b_i,ws=ws_i,delta=delta,lambda=lambda,HI = fisher_info)
})

HI_BN_Tx_9 <- rep(NA, dim(combinations)[1])
HI_BN_Tx_15 <- rep(NA, dim(combinations)[1])
HI_BN_Tx_51 <- rep(NA, dim(combinations)[1])
for (i in 1: dim(combinations)[1]) {
  HI_BN_Tx_9[i] <- fisher_results_BN_PMF[[i]]$HI[1]
  HI_BN_Tx_15[i] <- fisher_results_BN_PMF[[i]]$HI[2]
  HI_BN_Tx_51[i] <- fisher_results_BN_PMF[[i]]$HI[3]
}

combinations$HI_BN_Tx_9 <- HI_BN_Tx_9
combinations$HI_BN_Tx_15 <- HI_BN_Tx_15
combinations$HI_BN_Tx_51 <- HI_BN_Tx_51



combinations_BN_Tx_9 <- combinations[,c("b","barN","HI_BN_Tx_9")]
combinations_BN_Tx_15 <- combinations[,c("b","barN","HI_BN_Tx_15")]
combinations_BN_Tx_51 <- combinations[,c("b","barN","HI_BN_Tx_51")]

# Draw the Plots -----------------#

min_value <- min(combinations_HG_Tx_9$HI_HG_Tx_9, combinations_BN_Tx_9$HI_BN_Tx_9, na.rm = TRUE)
max_value <- max(combinations_HG_Tx_9$HI_HG_Tx_9, combinations_BN_Tx_9$HI_BN_Tx_9, na.rm = TRUE)


Tx_9_HG <- ggplot(combinations_HG_Tx_9, aes(x = barN, y = b, fill = HI_HG_Tx_9)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 9)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Hellinger information–HG, " ~ T[x] ~ "= 9")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


Tx_9_BN <- ggplot(combinations_BN_Tx_9, aes(x = barN, y = b, fill = HI_BN_Tx_9)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 9)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Hellinger information–BN, " ~ T[x] ~ "= 9")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


min_value <- min(combinations_HG_Tx_15$HI_HG_Tx_15, combinations_BN_Tx_15$HI_BN_Tx_15, na.rm = TRUE)
max_value <- max(combinations_HG_Tx_15$HI_HG_Tx_15, combinations_BN_Tx_15$HI_BN_Tx_15, na.rm = TRUE)

Tx_15_HG <- ggplot(combinations_HG_Tx_15, aes(x = barN, y = b, fill = HI_HG_Tx_15)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 15)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Hellinger information–HG, " ~ T[x] ~ "= 15")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


Tx_15_BN <- ggplot(combinations_BN_Tx_15, aes(x = barN, y = b, fill = HI_BN_Tx_15)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 15)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Hellinger information–BN, " ~ T[x] ~ "= 15")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


min_value <- min(combinations_HG_Tx_51$HI_HG_Tx_51, combinations_BN_Tx_51$HI_BN_Tx_51, na.rm = TRUE)
max_value <- max(combinations_HG_Tx_51$HI_HG_Tx_51, combinations_BN_Tx_51$HI_BN_Tx_51, na.rm = TRUE)


Tx_51_HG <- ggplot(combinations_HG_Tx_51, aes(x = barN, y = b, fill = HI_HG_Tx_51)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 51)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Hellinger information–HG, " ~ T[x] ~ "= 51")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))




Tx_51_BN <- ggplot(combinations_BN_Tx_51, aes(x = barN, y = b, fill = HI_BN_Tx_51)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, n.breaks = 10,limits = c(min_value, max_value),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 51)"
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Hellinger information–BN, " ~ T[x] ~ "= 51")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))



ggarrange(Tx_9_HG,Tx_9_BN,Tx_15_HG,Tx_15_BN,Tx_51_HG,Tx_51_BN,ncol = 2,nrow = 3)

ggsave("Paper Figures/Figure 3.png",width = 8,height = 12,units = "in",dpi = 300)


# Draw the Plots : In Standard Error Format -----------------#

combinations_HG_Tx_9$SE_HG_Tx_9 <- 1 / sqrt(combinations_HG_Tx_9$HI_HG_Tx_9)
combinations_BN_Tx_9$SE_BN_Tx_9 <- 1 / sqrt(combinations_BN_Tx_9$HI_BN_Tx_9)

step <- 1
min_value <- floor(min(combinations_HG_Tx_9$SE_HG_Tx_9, combinations_BN_Tx_9$SE_BN_Tx_9, na.rm = TRUE) / step) * step
max_value <- ceiling(max(combinations_HG_Tx_9$SE_HG_Tx_9, combinations_BN_Tx_9$SE_BN_Tx_9, na.rm = TRUE) / step) * step
Tx_9_quantile_levels <- c(3,4,5,6,8,10,12,14,16)


Tx_9_HG <- ggplot(combinations_HG_Tx_9, aes(x = barN, y = b, fill = SE_HG_Tx_9)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(2, 22, by = 2),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 9)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_HG_Tx_9,
               aes(x = barN, y = b, z = SE_HG_Tx_9),
               breaks = Tx_9_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_HG_Tx_9,
                    aes(x = barN, y = b, z = SE_HG_Tx_9),
                    breaks = Tx_9_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Standard Deviation – HG, " ~ T[x] ~ "= 9")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


Tx_9_BN <- ggplot(combinations_BN_Tx_9, aes(x = barN, y = b, fill = SE_BN_Tx_9)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(4, 22, by = 2),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 9)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_BN_Tx_9,
               aes(x = barN, y = b, z = SE_BN_Tx_9),
               breaks = Tx_9_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_BN_Tx_9,
                    aes(x = barN, y = b, z = SE_BN_Tx_9),
                    breaks = Tx_9_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Standard Deviation – BN, " ~ T[x] ~ "= 9")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


combinations_HG_Tx_15$SE_HG_Tx_15 <- 1 / sqrt(combinations_HG_Tx_15$HI_HG_Tx_15)
combinations_BN_Tx_15$SE_BN_Tx_15 <- 1 / sqrt(combinations_BN_Tx_15$HI_BN_Tx_15)

step <- 1
min_value <- floor(min(combinations_HG_Tx_15$SE_HG_Tx_15, combinations_BN_Tx_15$SE_BN_Tx_15, na.rm = TRUE) / step) * step
max_value <- ceiling(max(combinations_HG_Tx_15$SE_HG_Tx_15, combinations_BN_Tx_15$SE_BN_Tx_15, na.rm = TRUE) / step) * step
Tx_15_quantile_levels <- c(4,5,6,8,10,12,14,16,18,20)


Tx_15_HG <- ggplot(combinations_HG_Tx_15, aes(x = barN, y = b, fill = SE_HG_Tx_15)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(5, 25, by = 5),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 15)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_HG_Tx_15,
               aes(x = barN, y = b, z = SE_HG_Tx_15),
               breaks = Tx_15_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_HG_Tx_15,
                    aes(x = barN, y = b, z = SE_HG_Tx_15),
                    breaks = Tx_15_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Standard Deviation – HG, " ~ T[x] ~ "= 15")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


Tx_15_BN <- ggplot(combinations_BN_Tx_15, aes(x = barN, y = b, fill = SE_BN_Tx_15)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(5, 25, by = 5),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 15)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_BN_Tx_15,
               aes(x = barN, y = b, z = SE_BN_Tx_15),
               breaks = Tx_15_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_BN_Tx_15,
                    aes(x = barN, y = b, z = SE_BN_Tx_15),
                    breaks = Tx_15_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Standard Deviation – BN, " ~ T[x] ~ "= 15")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

combinations_HG_Tx_51$SE_HG_Tx_51 <- 1 / sqrt(combinations_HG_Tx_51$HI_HG_Tx_51)
combinations_BN_Tx_51$SE_BN_Tx_51 <- 1 / sqrt(combinations_BN_Tx_51$HI_BN_Tx_51)

step <- 1
min_value <- floor(min(combinations_HG_Tx_51$SE_HG_Tx_51, combinations_BN_Tx_51$SE_BN_Tx_51, na.rm = TRUE) / step) * step
max_value <- ceiling(max(combinations_HG_Tx_51$SE_HG_Tx_51, combinations_BN_Tx_51$SE_BN_Tx_51, na.rm = TRUE) / step) * step
Tx_51_quantile_levels <- c(5,6,8,10,12,14,16,18,20,25,30,40,50,60)


Tx_51_HG <- ggplot(combinations_HG_Tx_51, aes(x = barN, y = b, fill = SE_HG_Tx_51)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(10, 70, by = 10),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 51)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_HG_Tx_51,
               aes(x = barN, y = b, z = SE_HG_Tx_51),
               breaks = Tx_51_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_HG_Tx_51,
                    aes(x = barN, y = b, z = SE_HG_Tx_51),
                    breaks = Tx_51_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-hypergeometric model",
       subtitle = expression("Standard Deviation – HG, " ~ T[x] ~ "= 51")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


Tx_51_BN <- ggplot(combinations_BN_Tx_51, aes(x = barN, y = b, fill = SE_BN_Tx_51)) +
  geom_tile() +                         # Heatmap with tiles
  scale_fill_viridis_c(option = "mako", direction = -1,name = NULL, breaks = seq(10, 70, by = 10),limits = c(min_value, max_value),
                       oob = scales::squish,guide = guide_colourbar(barheight = unit(6, "cm"))) +  # , name = "H(Tx = 51)"
  # Superimpose contours: Constant Standard Error
  geom_contour(data = combinations_BN_Tx_51,
               aes(x = barN, y = b, z = SE_BN_Tx_51),
               breaks = Tx_51_quantile_levels,
               color = "black", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +

  geom_text_contour(data = combinations_BN_Tx_51,
                    aes(x = barN, y = b, z = SE_BN_Tx_51),
                    breaks = Tx_51_quantile_levels,
                    skip = 0,          # Do not skip any contours
                    size = 3,          # Smaller text helps fit
                    stroke = 0.1,      # Thin outline for readability
                    inherit.aes = FALSE) +
  labs(x = expression("Group size (" ~ bar(N) ~ ")"),
       y = expression("Number of groups (" ~ b ~ ")"),
       title = "Grouped-binomial model",
       subtitle = expression("Standard Deviation – BN, " ~ T[x] ~ "= 51")) +
  theme_minimal(base_size = 14) +       # Clean theme
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


ggarrange(Tx_9_HG,Tx_9_BN,Tx_15_HG,Tx_15_BN,Tx_51_HG,Tx_51_BN,ncol = 2,nrow = 3)

ggsave("Paper Figures/Figure 3 Updated.png",width = 8,height = 12,units = "in",dpi = 300)



