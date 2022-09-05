### SETUP PLOTTING -------------------------------------------------------------
## colors
patch_colors <- c("#1b9e77", "#d95f02")
# cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))
control_type_labs <- c("uniform", "non-uniform") # << change name of unique here
names(control_type_labs) <- c("uniform", "unique")

# define labels
chol_control_labs <- c("Vaccination effort", "Sanitation effort")
names(chol_control_labs) <- c("v", "u")
# cholera states
chol_I_labs <- c("Infectives in Patch 1", "Infectives in Patch 2")
names(chol_I_labs) <- 1:2

# define labels
ebola_control_labs <- c("Vaccination effort", "Hospitalization effort")
names(ebola_control_labs) <- c("v", "u")
# ebola states
ebola_I_labs <- c("Infectives in Patch 1", "Infectives in Patch 2")
names(ebola_I_labs) <- 1:2
# linetypes 
ebola_ltys <- c("solid",  "23", "solid","28")
names(ebola_ltys) <- c("1-unique", "1-uniform",  "2-unique","2-uniform")
# colors
ebola_cols <- sort(rep(patch_colors,2))
names(ebola_cols) <- c("1-unique", "1-uniform", "2-unique", "2-uniform")
# control labels
ebola_gov_labs <- c( "Patch 1: uniform", "Patch 1: non-uniform","Patch 2: uniform", "Patch 2: non-uniform")
names(ebola_gov_labs) <- c("1-uniform", "1-unique", "2-uniform", "2-unique")


#### TIME SERIES PLOT FUNCTIONS ------------------------------------------------

# plot_df should have columns for time, value, patch, control_type
plot_timeseries <- function(plot_df, patch_colors, facet_type, facet_labs,
                            lty_lab, y_lab, leg_pos, patch_ltys) {
  # plot
  p <- ggplot(
    data = plot_df,
    aes(
      x = time, y = value,
      color = patch_control, linetype = patch_control
    )
  ) +
    geom_line() +
    scale_color_manual(values = patch_colors, labels = lty_lab, name = "Control:") +
    scale_linetype_manual(values = patch_ltys, labels = lty_lab, name = "Control:") +
    scale_x_continuous( expand = expansion(mult = c(0,0.05))) +
    scale_y_continuous( expand = expansion(mult = c(0,0.05))) +
    labs(
      x = "Days since start of control",
      y = y_lab,
      linetype = "Control type:",
      size = "Control type:"
    ) +
    theme_minimal_grid() +
    theme(
      legend.position = leg_pos,
      legend.key.width = unit(0.8, "cm"),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement = "outside"
    ) +
    panel_border()
  if (facet_type == "states") {
    p <- p +
      facet_wrap(vars(variable),
                 labeller = labeller(variable = facet_labs),
                 nrow = 1,
                 scales = "free", strip.position = "top"
      )
  } else if (facet_type == "Patch:") {
    p <- p +
      facet_wrap(vars(patch),
                 labeller = labeller(patch = facet_labs),
                 nrow = 2,
                 scales = "free",
                 strip.position = "top"
      )
  }
  return(p)
}

create_multipanel_ts_plot <- function(model_name, states, patch_colors,
                                      I_labs, control_labs, control_type_labs, patch_ltys) {
  patch_ltys2 = c(patch_ltys[1:3], patch_ltys[2])
  names(patch_ltys2) = names(patch_ltys)
  # plot infectious class
  p_Istates <- states %>%
    filter(
      model == model_name,
      variable == "I"
    ) %>%
    plot_timeseries(
      patch_colors = patch_colors,
      facet_labs = I_labs,
      facet_type = "Patch:",
      lty_lab = control_type_labs,
      y_lab = "",
      leg_pos = "bottom", 
      patch_ltys = patch_ltys2
    )
  # plot controls
  p_controls <- states %>%
    filter(
      model == model_name,
      variable %in% c("v", "u")
    ) %>%
    plot_timeseries(
      patch_colors = patch_colors,
      facet_labs = control_labs,
      facet_type = "states",
      lty_lab = control_type_labs,
      y_lab = "",
      leg_pos = "none", 
      patch_ltys = patch_ltys
    )
  # put together
  l <- get_legend(p_Istates)
  p <- plot_grid(p_Istates + 
                   theme(legend.position = "none"), 
                 p_controls,
                 rel_widths = c(0.3, 0.7),
                 align = "h",
                 axis = "b"
  )
  b <- ggplot + theme_nothing()
  p <- plot_grid(p, 
                 plot_grid(b,l,b, nrow = 1), 
                 rel_heights = c(0.9,0.1), ncol = 1)
  return(p)
}