# Create vectors
# Fertility
# KI cat birth rates matrix, data for female offsping produced each year.
# Data from Budke, C & Slater, M (2009)
fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)

# Fertility errors based on Budke & Slater
# Mean and standard deviations, juvenile fertility:
juv_m_sd <- mean(c(((0.745 / 3 - 0.352 / 3) / 2), ((1.58 / 3 - 0.745 / 3) / 2)))
fy_m_sd <- mean(c(((0.745 - 0.352) / 2), ((1.58 - 0.745) / 2))) # Mean and standard deviations, juvenile fertility
a_m_sd <- mean(c(((2.52 - 1.98) / 2), ((3.78 - 2.52) / 2))) # Mean and standard deviations, adult fertility
# Mean and standard deviations vector, juvenile and adult fertility:
std_fertility <- c(0.18 * fertility[1], 0.18 * fertility[2], a_m_sd, a_m_sd, a_m_sd, a_m_sd, a_m_sd)

# Survival
# KI cat survival
# probability of surviving from one year to the next. e.g surviving fourth year of life
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)

# survival errors based on Budke & Slater
y1_2_s_sd <- mean(c(((0.46 - 0.27) / 2), ((0.73 - 0.46) / 2))) # mean and standard deviations, juvenile survival
a_s_sd <- mean(c(((0.7 - 0.55) / 2), ((0.78 - 0.7) / 2))) # mean and standard deviations, adult survival
# Mean and standard deviations vector, juvenile and adult survival:
std_survival_probability <- c(y1_2_s_sd, y1_2_s_sd, a_s_sd, a_s_sd, a_s_sd, a_s_sd)