# =============================================================================
# DOES SOIL INTERACT WITH SYMBIONT STATUS TO INFLUENCE DEMOGRAPHY?
# Full Bayesian Analysis via rstan
# Author: Jacob Moutouama
# =============================================================================


# ggsave(
#   file.path(fig_dir, "PPC_models_soilendo.pdf"),
#   plot   = ppc_panel,
#   width  = 170,
#   height = 100,
#   units  = "mm",
#   dpi    = 600,
#   device = cairo_pdf
# )

# =============================================================================
# PART 3 — POSTERIOR PREDICTIONS
# =============================================================================

# ── 3.1  Extract posterior fitted means across precipitation gradient ─────────
extract_posterior_precip <- function(fit, sd_obj, link = "identity",
                                     n_grid = 50, ref_df = NULL) {
  draws_beta <- as.matrix(fit, pars = "beta")
  draws_upop <- as.matrix(fit, pars = "u_pop")
  d          <- sd_obj$d_used
  # Use full reference dataset for back-transform if supplied; fall back to d.
  ref        <- if (!is.null(ref_df)) ref_df else d
  n_iter     <- nrow(draws_beta)

  precip_seq <- seq(min(d$precip_std), max(d$precip_std), length.out = n_grid)

  map_dfr(endo_levels, function(sym) {
    sym_val <- if (sym == "S+") 1 else 0
    map_dfr(seq_along(precip_seq), function(i) {
      p <- precip_seq[i]
      # Design vector: intercept, Symbiont, precip_std, age_std=0, Symbiont:precip
      xvec <- c(1, sym_val, p, 0, sym_val * p)
      eta  <- draws_beta %*% xvec +
        rowMeans(draws_upop)   # marginalise over populations
      mu <- switch(link, log = exp(eta), identity = eta)
      tibble(
        Symbiont    = sym,
        precip_std  = p,
        # FIX: use ref$precip_mean (full dataset) not d$precip_mean (filtered)
        precip_mean = p * sd(ref$precip_mean) + mean(ref$precip_mean),
        draw        = seq_len(n_iter),
        mu          = as.numeric(mu)
      )
    })
  })
}

# Pass aghysoils as ref_df so back-transform uses the full dataset's mean/SD.
post_precip_biomass <- extract_posterior_precip(fit_biomass, sd_biomass, "identity",
                                                ref_df = aghysoils) %>%
  mutate(mu = exp(mu))
post_precip_inflo   <- extract_posterior_precip(fit_inflo,   sd_inflo,   "log",
                                                ref_df = aghysoils)
post_precip_totspi  <- extract_posterior_precip(fit_totspi,  sd_totspi,  "log",
                                                ref_df = aghysoils)
post_precip_avgspi  <- extract_posterior_precip(fit_avgspi,  sd_avgspi,  "log",
                                                ref_df = aghysoils)

# =============================================================================
# PART 4 — FIGURES
# =============================================================================
# Layout: 3-column combined figure.
#   Col 1 — Biomass:      upper (predicted) / lower (Δ contrast), 2:1 height
#   Col 2 — Inflo:        upper (predicted) / lower (Δ contrast), 2:1 height
#   Col 3 — Coefficients: caterpillar for focal betas, faceted by response
#
# Panel tags (A–E) placed OUTSIDE plots via plot_annotation(tag_levels).
# Each sub-plot carries tag = "" so patchwork assigns the letter externally.
# =============================================================================

# ── 4.1  Data-building helpers ────────────────────────────────────────────────

build_upper_data <- function(post_df) {
  post_df %>%
    group_by(Symbiont, precip_mean) %>%
    summarise(
      med  = median(mu),
      lo95 = quantile(mu, 0.025),
      hi95 = quantile(mu, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Symbiont = factor(Symbiont, levels = endo_levels))
}

build_lower_data <- function(post_df) {
  post_df %>%
    group_by(precip_mean, draw) %>%
    pivot_wider(names_from  = Symbiont,
                values_from = mu,
                id_cols     = c(precip_mean, draw)) %>%
    mutate(contrast = `S+` - `S-`) %>%
    group_by(precip_mean) %>%
    summarise(
      med  = median(contrast),
      lo95 = quantile(contrast, 0.025),
      hi95 = quantile(contrast, 0.975),
      .groups = "drop"
    )
}

build_obs_data <- function(raw_df, response_col) {
  raw_df %>%
    filter(!is.na(.data[[response_col]]), .data[[response_col]] > 0) %>%
    mutate(
      Symbiont = factor(Symbiont, levels = endo_levels),
      y_obs    = .data[[response_col]]
    ) %>%
    dplyr::select(precip_mean, Symbiont, y_obs)
}

# ── 4.2  Upper plot ───────────────────────────────────────────────────────────
make_upper_plot <- function(post_df, raw_df, response_col, y_label,
                            show_legend = TRUE) {
  pd  <- build_upper_data(post_df)
  obs <- build_obs_data(raw_df, response_col)

  ggplot() +
    geom_ribbon(
      data = pd,
      aes(x = precip_mean, ymin = lo95, ymax = hi95,
          fill = Symbiont, group = Symbiont),
      alpha = 0.2, color = NA
    ) +
    geom_line(
      data = pd,
      aes(x = precip_mean, y = med,
          color = Symbiont, group = Symbiont),
      linewidth = 0.6
    ) +
    geom_point(
      data = obs,
      aes(x = precip_mean, y = y_obs, color = Symbiont),
      alpha = 0.45, size = 2.5, shape = 16,
      show.legend = FALSE
    ) +
    scale_color_manual(values = ENDO_COLORS, breaks = endo_levels,
                       labels = ENDO_LABELS, name = "Symbiont") +
    scale_fill_manual(values  = ENDO_COLORS, breaks = endo_levels,
                      labels = ENDO_LABELS, name = "Symbiont") +
    labs(x = NULL, y = y_label) +
    vr_theme() +
    theme(
      legend.position = if (show_legend) c(0.2, 0.7) else "none",
      legend.key.size = unit(8, "pt"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    )
}

# ── 4.3  Lower plot ───────────────────────────────────────────────────────────
make_lower_plot <- function(post_df,
                            x_label = "Soil origin MAP (mm/yr)") {
  pd <- build_lower_data(post_df)

  ggplot(pd, aes(x = precip_mean)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95),
                fill = "#D9BFD6", alpha = 0.45, color = NA) +
    geom_line(aes(y = med), color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.4, color = "black") +
    labs(x = x_label,
         y = expression(Delta ~ (italic(S)^"+" - italic(S)^"\u2212"))) +
    vr_theme() +
    theme(legend.position = "none")
}

# ── 4.4  Reusable column builder (kept for Part 4b / spikelets figure) ───────
make_vr_column <- function(post_df, raw_df, response_col, y_label,
                           show_legend = TRUE) {
  upper <- make_upper_plot(post_df, raw_df, response_col, y_label,
                           show_legend = show_legend)
  lower <- make_lower_plot(post_df)
  upper / lower + plot_layout(heights = c(2, 1))
}

# ── 4.5  Build predicted/Δ leaf plots ─────────────────────────────────────────
# Kept as flat (non-nested) plot objects, rather than pre-combined columns, so
# the final 3-row assembly in 4.7 can control every row's height explicitly
# instead of two different internal height ratios (2:1 vs 1:1) fighting when
# nested inside one outer plot_layout().
upper_biomass <- make_upper_plot(post_precip_biomass, aghysoils,
                                 "calc_abg_mass_tot", "Aboveground biomass (g)",
                                 show_legend = TRUE)
lower_biomass <- make_lower_plot(post_precip_biomass)

upper_inflo <- make_upper_plot(post_precip_inflo, aghysoils,
                               "calc_total_inflo", "Total inflorescences",
                               show_legend = FALSE)
lower_inflo <- make_lower_plot(post_precip_inflo)

# ── 4.6  Coefficient caterpillar (third column) ───────────────────────────────
FOCAL_IDX    <- c(2L, 3L, 5L)
FOCAL_BREAKS <- paste0("beta[", FOCAL_IDX, "]")
FOCAL_LABS <- expression(
  Symbiont ~ (italic(S)^"+"),
  "Precipitation (std)",
  "Symbiont" ~ "\u00d7" ~ "Precipitation"
)
extract_focal_betas <- function(fit, model_label) {
  as.matrix(fit, pars = "beta")[, FOCAL_IDX, drop = FALSE] %>%
    as.data.frame() %>%
    setNames(paste0("beta[", FOCAL_IDX, "]")) %>%
    pivot_longer(everything(),
                 names_to  = "parameter",
                 values_to = "estimate") %>%
    mutate(model = model_label)
}

coef_long <- bind_rows(
  extract_focal_betas(fit_biomass, "Aboveground\nbiomass (g)"),
  extract_focal_betas(fit_inflo,   "Total\ninflorescences")
  # extract_focal_betas(fit_totspi,  "Total\nspikelets (projected)")
)

coef_summary <- coef_long %>%
  group_by(model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    prob_gt0   = mean(estimate > 0),
    prob_lt0   = mean(estimate < 0),
    .groups    = "drop"
  ) %>%
  mutate(
    strong_effect = prob_gt0 > 0.9 | prob_lt0 > 0.9,
    parameter = factor(parameter,
                       levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    model     = factor(model,
                       levels = c("Aboveground\nbiomass (g)",
                                  "Total\ninflorescences",
                                  "Total\nspikelets (projected)")),
    term_type = if_else(parameter == "beta[5]", "Interaction", "Main effect"),
    term_type = factor(term_type, levels = names(TERM_COLORS)),
    pt_fill   = if_else(strong_effect, as.character(term_type), "white")
  )

n_coef      <- nlevels(coef_summary$parameter)
shade_rects <- Filter(Negate(is.null), lapply(seq_len(n_coef), function(i) {
  if (i %% 2 == 1)
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = i - 0.49, ymax = i + 0.49,
             fill = "grey90", color = NA, alpha = 0.55)
}))

fig_coef <- ggplot(coef_summary, aes(y = parameter, color = term_type)) +
  shade_rects +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90),
                linewidth = 1.2, width = 0) +
  geom_point(aes(x = median_est, fill = pt_fill),
             shape = 21, size = 3.0, stroke = 0.9) +
  scale_fill_manual(values = c(TERM_COLORS, "white" = "white"),
                    guide  = "none") +
  scale_color_manual(
    values = TERM_COLORS,
    name   = NULL,
    guide  = guide_legend(
      override.aes = list(
        fill   = unname(TERM_COLORS),
        shape  = 21, size = 3, stroke = 0.9
      )
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  facet_wrap(~ model, ncol = 1, scales = "free_x",
             strip.position = "right") +
  scale_y_discrete(breaks = FOCAL_BREAKS, labels = FOCAL_LABS,
                   expand = expansion(add = 0.6))  +
  labs(x = "Posterior coefficient", y = NULL) +
  vr_theme() +
  theme(
    legend.position      = "none",
    legend.key.size      = unit(8, "pt"),
    panel.grid.major.y   = element_blank(),
    strip.text.y.right   = element_text(size = 10, color = "black", angle = 90,
                                        lineheight = 0.85),
    axis.text            = element_text(size = 8),
    axis.title           = element_text(size = 10)
  )

# ── 4.7  Assemble 3-row figure with external panel tags ──────────────────────
# Row 1 (upper, predicted):     upper_biomass | upper_inflo        (= a, b)
# Row 2 (lower, Δ contrast):    lower_biomass | lower_inflo        (= c, d)
# Row 3 (coefficient panel):    fig_coef, spanning the full width  (= e)
#
# All three rows are flat leaf-level plots assembled in ONE plot_layout(), so
# heights = c(2, 1, 2) is the single, explicit source of truth for row sizing
# (2:1 matches the original upper:lower ratio; the coef row gets a comparable
# weight to the top row rather than inheriting a mismatched internal ratio).
#
# plot_annotation(tag_levels = "a") walks every leaf plot in patchwork order
# and assigns a, b, c, d, e from left to right, top to bottom.
# plot.tag.position = "topleft" puts each letter outside and above the panel.
fig_greenhouse_combined <-
  (upper_biomass | upper_inflo) /
  (lower_biomass | lower_inflo) /
  ((plot_spacer() | fig_coef | plot_spacer()) +
     plot_layout(widths = c(-0.27, 2, -0.06))) +
  plot_layout(heights = c(2, 1, 3)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 8, face = "plain"),
    plot.tag.position = "topleft"
  )
print(fig_greenhouse_combined)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Combined.pdf"),
  plot   = fig_greenhouse_combined,
  width  = 170,
  height = 150,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Fig_Greenhouse_Combined.tiff"),
  plot        = fig_greenhouse_combined,
  width       = 170,
  height      = 150,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# PART 4b — STANDALONE SPIKELETS FIGURE
# =============================================================================
# Two-panel figure for total spikelets (projected):
#   Left  — col_spike: upper (posterior predicted lines + observed points) /
#                      lower (Δ contrast S+ − S−), stacked 2:1
#   Right — fig_coef_spike: caterpillar plot for the three focal betas from
#                            fit_totspi (beta[2], beta[3], beta[5])
#
# This is the spikelet-only analogue of Fig_Greenhouse_Combined and is saved
# separately as Fig_Spikelets_Combined.pdf.
# =============================================================================

# ── 4b.1  Spikelet predicted column ──────────────────────────────────────────
col_spike <- make_vr_column(
  post_df      = post_precip_totspi,
  raw_df       = aghysoils,
  response_col = "calc_total_spik",
  y_label      = "Total spikelets",
  show_legend  = TRUE
)

# ── 4b.2  Spikelet coefficient caterpillar ────────────────────────────────────
# Extract focal betas from fit_totspi only.
coef_spike <- extract_focal_betas(fit_totspi, "Total\nspikelets (projected)")

coef_summary_spike <- coef_spike %>%
  group_by(model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    prob_gt0   = mean(estimate > 0),
    prob_lt0   = mean(estimate < 0),
    .groups    = "drop"
  ) %>%
  mutate(
    strong_effect = prob_gt0 > 0.9 | prob_lt0 > 0.9,
    parameter = factor(parameter,
                       levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    term_type = if_else(parameter == "beta[5]", "Interaction", "Main effect"),
    term_type = factor(term_type, levels = names(TERM_COLORS)),
    pt_fill   = if_else(strong_effect, as.character(term_type), "white")
  )

n_coef_spike  <- nlevels(coef_summary_spike$parameter)
shade_rects_spike <- Filter(Negate(is.null), lapply(seq_len(n_coef_spike), function(i) {
  if (i %% 2 == 1)
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = i - 0.49, ymax = i + 0.49,
             fill = "grey93", color = NA, alpha = 0.55)
}))

fig_coef_spike <- ggplot(coef_summary_spike, aes(y = parameter, color = term_type)) +
  shade_rects_spike +
  geom_errorbar(aes(xmin = lo95, xmax = hi95), linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90), linewidth = 1.2,  width = 0) +
  geom_point(aes(x = median_est, fill = pt_fill),
             shape = 21, size = 3.0, stroke = 0.9) +
  scale_fill_manual(values = c(TERM_COLORS, "white" = "white"), guide = "none") +
  scale_color_manual(
    values = TERM_COLORS,
    name   = NULL,
    guide  = guide_legend(
      override.aes = list(fill = unname(TERM_COLORS),
                          shape = 21, size = 3, stroke = 0.9)
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  scale_y_discrete(labels = FOCAL_LABS, expand = expansion(add = 0.6)) +
  labs(x = "Posterior coefficient (spikelets)", y = NULL) +
  vr_theme() +
  theme(
    legend.position    = "right",
    legend.key.size    = unit(8, "pt"),
    panel.grid.major.y = element_blank()
  )

# ── 4b.3  Assemble and save spikelet figure ───────────────────────────────────
fig_spikelets_combined <-
  (col_spike | fig_coef_spike) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    legend.position = "none",
    axis.text         = element_text(size = 10),
    plot.tag          = element_text(size = 12, face = "plain"),
    plot.tag.position = "topleft"
  )

print(fig_spikelets_combined)

ggsave(
  file.path(fig_dir, "Fig_Spikelets_Combined.pdf"),
  plot   = fig_spikelets_combined,
  width  = 173,
  height = 82.8,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "Fig_Spikelets_Combined.tiff"),
  plot        = fig_spikelets_combined,
  width       = 173,
  height      = 82.8,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# PART 5 — STATISTICAL SUMMARY TABLES
# =============================================================================

posterior_table <- function(fit, sd_obj, model_label) {
  coef_names <- colnames(sd_obj$X)
  coef_names <- gsub("\\(Intercept\\)",      "Intercept",              coef_names)
  coef_names <- gsub("SymbiontS\\+",         "Sym[S+]",                coef_names)
  coef_names <- gsub("precip_std",           "Precip (std)",           coef_names)
  coef_names <- gsub("Sym\\[S\\+\\]:Precip", "Sym[S+] \u00d7 Precip", coef_names)
  coef_names <- gsub("age_std",              "Age (std)",              coef_names)

  draws_beta <- as.data.frame(as.matrix(fit, pars = "beta"))
  colnames(draws_beta) <- coef_names

  draws_beta %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "draw") %>%
    group_by(Parameter) %>%
    summarise(
      Median  = round(median(draw), 3),
      Mean    = round(mean(draw),   3),
      SD      = round(sd(draw),     3),
      `Q2.5`  = round(quantile(draw, 0.025), 3),
      `Q97.5` = round(quantile(draw, 0.975), 3),
      `P(>0)` = round(mean(draw > 0), 3),
      .groups = "drop"
    ) %>%
    mutate(Model = model_label, .before = 1)
}

tab_biomass <- posterior_table(fit_biomass, sd_biomass, "Biomass (Gaussian)")
tab_inflo   <- posterior_table(fit_inflo,   sd_inflo,   "Inflo count (NB)")
tab_totspi  <- posterior_table(fit_totspi,  sd_totspi,  "Total spikelets projected (NB)")
tab_avgspi  <- posterior_table(fit_avgspi,  sd_avgspi,  "Avg spikelets (NB)")

all_tabs <- bind_rows(tab_biomass, tab_inflo, tab_totspi, tab_avgspi)

cat("\n\n")
cat("=============================================================\n")
cat(" BAYESIAN POSTERIOR SUMMARY — All models\n")
cat(" Key parameter: Sym[S+] x Precip\n")
cat(" Null expectation for greenhouse: this coefficient ~ 0\n")
cat(" (soil-origin precip does not modify endophyte effect\n")
cat("  when plants grow in a common environment)\n")
cat("=============================================================\n\n")
for (mod_label in unique(all_tabs$Model)) {
  cat("--- ", mod_label, " ---\n", sep = "")
  all_tabs %>% filter(Model == mod_label) %>% dplyr::select(-Model) %>% print(n = Inf)
  cat("\n")
}

# write.csv(all_tabs,
#           file      = file.path(output_dir, "Table1_PosteriorSummary_precip.csv"),
#           row.names = FALSE)


# ── Full breakdown (for supplement) ───────────────────────────────────────────

ss_biomass <- aghysoils %>%
  filter(!is.na(calc_abg_mass_tot)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Biomass")

ss_inflo <- aghysoils %>%
  filter(!is.na(calc_total_inflo)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Inflorescence")

ss_spike <- aghysoils %>%
  filter(!is.na(calc_avg_spikelet)) %>%
  count(Site, Symbiont) %>%
  mutate(Response = "Avg spikelet")

sample_size_full <- bind_rows(ss_biomass, ss_inflo, ss_spike) %>%
  tidyr::complete(Response, Site, Symbiont, fill = list(n = 0)) %>%
  arrange(Response, Site)

# write.csv(sample_size_full,
#           file = file.path(output_dir, "Table_Sx_SampleSizes.csv"),
#           row.names = FALSE)

# =============================================================================
# PART 6 — SENSITIVITY ANALYSIS: DROP KER & SON
# =============================================================================
# Tom's question on the soil-origin precipitation result (plants from wetter-
# origin soils had more biomass/inflorescences even under common greenhouse
# conditions): is that signal being driven by KER and SON specifically (e.g.
# their limestone-derived soils), or does a precipitation/climate signal
# persist even without those two sites?
#
# We refit each response model (same formula, same priors, same compiled Stan
# programs as PART 2) on the subset of the data with KER and SON removed, and
# compare the posterior for the precip_std main effect (beta[3]) and the
# Symbiont x precip_std interaction (beta[5]) to the full-data models fit
# earlier. If the precip_std effect holds up (median away from 0, similar
# sign/magnitude) with KER/SON dropped, that supports a genuine climate
# signal rather than an artifact of those two sites' soils.
# =============================================================================

aghysoils_sens <- aghysoils %>%
  filter(!Site %in% c("KER", "SON")) %>%
  droplevels()

message(sprintf(
  "Sensitivity subset (KER & SON dropped): %d of %d rows retained",
  nrow(aghysoils_sens), nrow(aghysoils)
))

# ── 6.1  Stan data + refit (reuses the same compiled models from PART 2) ─────
sd_biomass_sens <- make_stan_data(aghysoils_sens, "calc_abg_mass_tot", "gaussian")
sd_inflo_sens   <- make_stan_data(aghysoils_sens, "calc_total_inflo",  "negbinom")
sd_totspi_sens  <- make_stan_data(aghysoils_sens, "calc_total_spik",   "negbinom")
sd_avgspi_sens  <- make_stan_data(
  aghysoils_sens %>% filter(!is.na(calc_avg_spikelet), calc_avg_spikelet > 0),
  "calc_avg_spikelet", "negbinom"
)

fit_biomass_sens <- run_stan(mod_gaussian, sd_biomass_sens, "gaussian")
fit_inflo_sens   <- run_stan(mod_negbinom, sd_inflo_sens,   "negbinom")
fit_totspi_sens  <- run_stan(mod_negbinom, sd_totspi_sens,  "negbinom")
fit_avgspi_sens  <- run_stan(mod_negbinom, sd_avgspi_sens,  "negbinom")

check_hmc(fit_biomass_sens, "M1s Biomass (no KER/SON)")
check_hmc(fit_inflo_sens,   "M2s Inflo (no KER/SON)")
check_hmc(fit_totspi_sens,  "M3s Total spikelets (no KER/SON)")
check_hmc(fit_avgspi_sens,  "M4s Avg spikelets (no KER/SON)")

# ── 6.2  Posterior summary table (same layout as Table 1) ────────────────────
tab_biomass_sens <- posterior_table(fit_biomass_sens, sd_biomass_sens,
                                    "Biomass, no KER/SON (Gaussian)")
tab_inflo_sens   <- posterior_table(fit_inflo_sens,   sd_inflo_sens,
                                    "Inflo count, no KER/SON (NB)")
tab_totspi_sens  <- posterior_table(fit_totspi_sens,  sd_totspi_sens,
                                    "Total spikelets, no KER/SON (NB)")
tab_avgspi_sens  <- posterior_table(fit_avgspi_sens,  sd_avgspi_sens,
                                    "Avg spikelets, no KER/SON (NB)")

all_tabs_sens <- bind_rows(tab_biomass_sens, tab_inflo_sens,
                           tab_totspi_sens, tab_avgspi_sens)

cat("\n\n")
cat("=============================================================\n")
cat(" SENSITIVITY ANALYSIS \u2014 KER & SON dropped\n")
cat(" Question: does the soil-origin precip signal survive without\n")
cat(" KER/SON (i.e. is it just their limestone-derived soils), or is\n")
cat(" there a climate signal even excluding those two sites?\n")
cat("=============================================================\n\n")
for (mod_label in unique(all_tabs_sens$Model)) {
  cat("--- ", mod_label, " ---\n", sep = "")
  all_tabs_sens %>% filter(Model == mod_label) %>% dplyr::select(-Model) %>% print(n = Inf)
  cat("\n")
}

write.csv(all_tabs_sens,
          file      = file.path(output_dir, "TableS_PosteriorSummary_precip_noKERSON.csv"),
          row.names = FALSE)

# ── 6.3  Full-data vs. sensitivity coefficient comparison figure ────────────
# Side-by-side caterpillar plot (biomass + inflorescence, the two responses
# Tom asked about) contrasting "All sites" vs "No KER/SON" posteriors for
# beta[2] (Symbiont), beta[3] (Precip), and beta[5] (Symbiont x Precip).
# Uses one flat color for all terms (see TERM_COLORS note above) so the
# comparison is read purely off dataset, not term type.
coef_full_vs_sens <- bind_rows(
  extract_focal_betas(fit_biomass,      "Aboveground\nbiomass (g)")  %>% mutate(dataset = "All sites"),
  extract_focal_betas(fit_inflo,        "Total\ninflorescences")    %>% mutate(dataset = "All sites"),
  extract_focal_betas(fit_biomass_sens, "Aboveground\nbiomass (g)")  %>% mutate(dataset = "No KER/SON"),
  extract_focal_betas(fit_inflo_sens,   "Total\ninflorescences")    %>% mutate(dataset = "No KER/SON")
)

coef_summary_sens <- coef_full_vs_sens %>%
  group_by(dataset, model, parameter) %>%
  summarise(
    median_est = median(estimate),
    lo90       = quantile(estimate, 0.05),
    hi90       = quantile(estimate, 0.95),
    lo95       = quantile(estimate, 0.025),
    hi95       = quantile(estimate, 0.975),
    .groups    = "drop"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = rev(paste0("beta[", FOCAL_IDX, "]"))),
    dataset   = factor(dataset, levels = c("All sites", "No KER/SON"))
  )

fig_sensitivity <- ggplot(coef_summary_sens,
                          aes(y = parameter, x = median_est, color = dataset)) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                position = position_dodge(width = 0.5),
                linewidth = 0.35, width = 0) +
  geom_errorbar(aes(xmin = lo90, xmax = hi90),
                position = position_dodge(width = 0.5),
                linewidth = 1.2, width = 0) +
  geom_point(position = position_dodge(width = 0.5), size = 2.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey25", linewidth = 0.45) +
  scale_color_manual(values = c("All sites" = "#4D4D4D", "No KER/SON" = "#0072B2"),
                     name = NULL) +
  facet_wrap(~ model, ncol = 1, scales = "free_x", strip.position = "right") +
  scale_y_discrete(labels = FOCAL_LABS, expand = expansion(add = 0.6)) +
  labs(x = "Posterior coefficient", y = NULL,
       title = "Sensitivity: all sites vs. KER & SON dropped") +
  vr_theme() +
  theme(
    legend.position     = "top",
    legend.key.size     = unit(8, "pt"),
    panel.grid.major.y  = element_blank(),
    strip.text.y.right  = element_text(size = 8, color = "black", angle = 90,
                                       lineheight = 0.85)
  )

print(fig_sensitivity)

ggsave(
  file.path(fig_dir, "FigS_Sensitivity_noKERSON.pdf"),
  plot   = fig_sensitivity,
  width  = 140,
  height = 120,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "FigS_Sensitivity_noKERSON.tiff"),
  plot        = fig_sensitivity,
  width       = 140,
  height      = 120,
  units       = "mm",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)

# =============================================================================
# End of script
# =============================================================================
