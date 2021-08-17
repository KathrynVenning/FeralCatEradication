# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# remove everything
rm(list = ls())

# libraries
library(ggplot2)


# Compensatory density feedback
# K = carry capacity
# Population rate of increase relative to carry capacity
# Larger distance between populationa and K = faster population growth
pop_found <- 1629
k_max <- 2 * pop_found
k_vec <- c(1, pop_found / 2, pop_found, 0.75 * k_max, k_max)
red_vec <- c(1, 0.965, 0.89, 0.79, 0.71)
k_red_dat <- data.frame(k_vec, red_vec)

ggplot(k_red_dat, aes(k_vec, red_vec)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  labs(x = "N", y = "Reduction factor")
ggsave("reports/figures/reduction_factor.jpg")
