# libraries
library(ggplot2)

source("Cholera_params.R")

# Assumptions:
#   1. v1 = v2 = 0
# Parameters:
#   1. In overleaf, "xi_i_in" is "xi_i" in R;
#                   "xi_i_out" is "nu_i" in R
#   2. lambda_i = 1/(gamma_i + mu_i + delta_i + n_i)
#   3. sigma_i = n_i*lambda_i
#   4. tau = 1/(1-(sigma1*sigma2)) (called "nu" in Overleaf)
#   5. eta_i = 1/(nu_i + rho_i)
#   6. omega_i = rho_i/(nu_i + rho_i)
chol_r0 <- function(params, N0_1, N0_2) {
  with(as.list(c(params, IC)), {
    # calculate intermediate parameters
    lambda1 <- 1 / (gamma1 + mu1 + delta1 + n1)
    lambda2 <- 1 / (gamma2 + mu2 + delta2 + n2)
    sigma1 <- n1 * lambda1
    sigma2 <- n2 * lambda2
    tau <- 1 / (1 - (sigma1 * sigma2))
    eta1 <- 1 / (nu1 + rho1)
    eta2 <- 1 / (nu2 + rho2)
    omega1 <- rho1 / (nu1 + rho1)
    omega2 <- rho2 / (nu2 + rho2)
    # calculate RIi and RWi
    RI1 <- beta_I1 * N0_1 * lambda1
    RI2 <- beta_I2 * N0_2 * lambda2
    RW1 <- beta_W1 * N0_1 * xi1 * lambda1 * eta1
    RW2 <- beta_W2 * N0_2 * xi2 * lambda2 * eta2
    # calculate components of R0
    R11 <- tau * (RI1 + RW1)
    R12 <- sigma2 * tau * (RI1 + RW1)
    R21 <- sigma1 * tau * (RI2 + RW2) + tau * omega1 * ((xi1 * lambda1) / (xi2 * lambda2)) * RW2
    R22 <- tau * (RI2 + RW2) + tau * sigma2 * omega1 * ((xi1 * lambda1) / (xi2 * lambda2)) * RW2
    return((R11 + R22 + sqrt((R11 - R22)^2 + 4 * R12 * R21)) / 2)
  })
}

test_params <- expand.grid(
  beta_I1 = c(0, 2.64e-05),
  beta_I2 = c(0, 2.64e-05),
  beta_W1 = c(1.21E-5, 1.21E-4, 1.21E-3),
  beta_W2 = c(1.21E-5, 1.21E-4, 1.21E-3)
)

test_params$R0_vals <- apply(test_params, 1, function(i) {
  testp <- params
  testp["beta_I1"] <- i[1]
  testp["beta_I2"] <- i[2]
  testp["beta_W1"] <- i[3]
  testp["beta_W2"] <- i[4]
  return(chol_r0(testp, 10000, 10000))
})

ggplot(data = test_params, aes(x = as.factor(beta_W1), y = as.factor(beta_W2), fill = R0_vals)) +
  geom_tile(alpha = 0.9) +
  geom_text(aes(label = round(R0_vals, 2))) +
  facet_grid(rows = vars(beta_I1), cols = vars(beta_I2), labeller = "label_both") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")
ggsave("figures/parameter_sweeps/sampleR0.pdf", width = 5, height = 5, units = "in")
