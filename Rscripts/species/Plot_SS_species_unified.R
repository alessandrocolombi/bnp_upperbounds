# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc, wd_unicatt, wd_g100, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd, "bnp_upperbounds/Rscripts/species")
setwd(wd)

# Packages ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(grid)))

# Plot options -------------------------------------------------------------
save_img = TRUE
img_width = 18
img_height = 10
cex_plot = 2.3  # overall expansion factor for plot text/elements
cex_axis = 2    # size of axis tick labels
cex_lab  = 2    # size of axis titles
cex_xlab = 1.5    # size of x-axis titles only
mycol = c("Freq" = "darkorange",
          "PD" = "darkred",
          "FDP" = "darkblue",
          "Dir-Multi" = "lightblue")
lgd_names = c("Oracle", "Freq", "PD", "FDP", "Dir-Multi")
alpha <- 0.05

# Choose which grid to plot: "Mfix" or "nfix" -----------------------------
plot_mode = "nfix" # <--- modify here

if (plot_mode == "Mfix") {
  x_grid = seq(100, 5000, by = 500)
  x_breaks = c(100, 1000, 2000, 3000, 4000, 5000)
  x_lab = "n"
  save_name_base = "save/SS_species_MfixMAP_"
  img_name = "img/SS_species_unified_Mfix.pdf"
} else if (plot_mode == "nfix") {
  x_grid = seq(50, 1000, by = 50)
  x_breaks = c(100, 300, 500, 700, 900)
  x_lab = "M"
  save_name_base = "save/SS_species_nfixMAP_"
  img_name = "img/SS_species_unified_nfix.pdf"
} else {
  stop("plot_mode must be either 'Mfix' or 'nfix'.")
}

# The 3x3 layout is:
# row 1 -> Zips
# row 2 -> Geom
# row 3 -> Unif/Neg.Bin
panel_specs = data.frame(
  row_label = c(rep("Zips", 3),
                rep("Geom", 3),
                rep("Unif/Neg.Bin", 3)),
  col_id = rep(c("1", "2", "3"), 3),
  exp_name = c("Zipfs", "Zipfs", "Zipfs",
               "Geom", "Geom", "Geom",
               "Uniform", "NegBin", "NegBin"),
  trim_params = c("09", "102", "2",
                  "085", "09", "095",
                  "NA", "1_0003", "5_0003"),
  stringsAsFactors = FALSE
)

panel_specs$row_label = factor(
  panel_specs$row_label,
  levels = c("Zips", "Geom", "Unif/Neg.Bin")
)
panel_specs$col_id = factor(panel_specs$col_id, levels = c("1", "2", "3"))

build_panel_df = function(exp_name, trim_params, row_label, col_id) {
  load(paste0(save_name_base, exp_name, "_", trim_params, ".Rdat"))

  oracle = sapply(ExpRes_list, function(x) quantile(x[, 1], 1 - alpha))
  expres_qnt = lapply(
    ExpRes_list,
    function(x) apply(x[, 2:5], 2, quantile, probs = c(0.025, 0.5, 0.975))
  )
  expres_qnt = simplify2array(expres_qnt) # 3 x 4 x length(x_grid)

  oracle_df = data.frame(
    x = x_grid,
    y = 1000 * oracle,
    method = "Oracle",
    row_label = row_label,
    col_id = col_id,
    stringsAsFactors = FALSE
  )

  method_names = c("Freq", "PD", "FDP", "Dir-Multi")
  lines_df = do.call(rbind, lapply(seq_along(method_names), function(j) {
    data.frame(
      x = x_grid,
      y = 1000 * expres_qnt[2, j, ],
      ymin = 1000 * expres_qnt[1, j, ],
      ymax = 1000 * expres_qnt[3, j, ],
      method = method_names[j],
      row_label = row_label,
      col_id = col_id,
      stringsAsFactors = FALSE
    )
  }))

  list(oracle = oracle_df, methods = lines_df)
}

panel_data = lapply(seq_len(nrow(panel_specs)), function(i) {
  build_panel_df(
    exp_name = panel_specs$exp_name[i],
    trim_params = panel_specs$trim_params[i],
    row_label = panel_specs$row_label[i],
    col_id = panel_specs$col_id[i]
  )
})

oracle_df = do.call(rbind, lapply(panel_data, `[[`, "oracle"))
methods_df = do.call(rbind, lapply(panel_data, `[[`, "methods"))

oracle_df$method = factor(oracle_df$method, levels = lgd_names)
methods_df$method = factor(methods_df$method, levels = lgd_names)
oracle_df$row_label = factor(oracle_df$row_label, levels = levels(panel_specs$row_label))
methods_df$row_label = factor(methods_df$row_label, levels = levels(panel_specs$row_label))
oracle_df$col_id = factor(oracle_df$col_id, levels = levels(panel_specs$col_id))
methods_df$col_id = factor(methods_df$col_id, levels = levels(panel_specs$col_id))

gg = ggplot() +
  geom_ribbon(
    data = methods_df,
    aes(x = x, ymin = ymin, ymax = ymax, fill = method),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data = methods_df,
    aes(x = x, y = y, color = method),
    linewidth = 1.1
  ) +
  geom_line(
    data = oracle_df,
    aes(x = x, y = y, color = method),
    linewidth = 1.1
  ) +
  facet_grid(
    rows = vars(row_label),
    cols = vars(col_id)
  ) +
  scale_color_manual(
    values = c("Oracle" = "black", mycol),
    breaks = lgd_names,
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Oracle" = "black", mycol),
    breaks = lgd_names,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_breaks,
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    x = x_lab,
    y = "1000 * bound"
  ) +
  guides(
    fill = "none",
    color = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_bw(base_size = 13 * cex_plot) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_blank(),
    strip.text.y.right = element_text(angle = 270, face = "bold", size = 12 * cex_plot),
    axis.text.x = element_text(size = 11 * cex_axis),
    axis.text.y = element_text(size = 11 * cex_axis),
    axis.title.x = element_text(size = 12 * cex_xlab),
    axis.title.y = element_text(size = 12 * cex_lab),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11 * cex_plot)
  )

if (plot_mode == "Mfix") {
  gg = gg +
    scale_y_log10(labels = label_number())
}

if (save_img) {
  ggsave(
    filename = img_name,
    plot = gg,
    width = img_width,
    height = img_height
  )
}

print(gg)
