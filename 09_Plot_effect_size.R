# Project: Endo range limits
# Purpose: Caterpillar plots of ALL fixed-effect model coefficients (logit/log
#          scale) for survival, growth, inflorescence, and spikelet models.
#          One plot per vital rate, each exported as its own PDF, plus a
#          combined 2×2 supplement figure.
#
#   Filled circle ● = strong posterior support (> 90 % of mass on one side of 0)
#   Open   circle ○ = uncertain (CI straddles zero)
#   Colour           = term type (Intercept / Main / Two-way / Three-way)
#   Thick bar        = 90 % credible interval
#   Thin whisker     = 95 % credible interval
#
# Authors: Jacob Moutouama
# Date last modified (Y-M-D):
#
# Stan models (companion .stan files):
#   Survival      : survival_l.stan        – Bernoulli-logit, logit scale
#   Growth        : growth_l.stan          – Normal,          log   scale
#   Inflorescence : inflorescence_l.stan   – Zero-inflated NegBin, log scale
#   Spikelet      : spikelet_l.stan        – NegBin, log scale (ELVI & POAU only)
#
# Shared linear predictor (Eq. 1):
#   eta = b0 + bclim*clim + bendo*endo + bherb*herb
#           + bendoclim*(endo*clim)  + bherbclim*(herb*clim)
#           + bendoherb*(endo*herb)  + bendoherbclim*(endo*herb*clim)
# + random effects: site-year, plot, source population (non-centered).
#
# Fitted stanfit objects loaded from Dropbox — fully self-contained.

rm(list = ls())

# ── Packages ──────────────────────────────────────────────────────────────────
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(ggpubr)
library(Cairo)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD FITTED STAN OBJECTS & EXTRACT POSTERIORS
# ══════════════════════════════════════════════════════════════════════════════
fit_surv_ppt <- readRDS(url(
  "https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1"
))
posterior_samples_survival <- rstan::extract(fit_surv_ppt)

fit_grow_ppt <- readRDS(url(
  "https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1"
))
posterior_samples_grow <- rstan::extract(fit_grow_ppt)

fit_inf_ppt <- readRDS(url(
  "https://www.dropbox.com/scl/fi/5lfkgq6d5a2vzx5h2t9yr/fit_inf_abio_bio_endo_linear.rds?rlkey=pfqvm7sg8un5c14slypn1aqa9&dl=1"
))
posterior_samples_inf <- rstan::extract(fit_inf_ppt)

fit_spik_ppt <- readRDS(url(
  "https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1"
))
posterior_samples_spik <- rstan::extract(fit_spik_ppt)

# ══════════════════════════════════════════════════════════════════════════════
# 2. COEFFICIENT METADATA
# ══════════════════════════════════════════════════════════════════════════════
all_coef_names <- c(
  "b0", "bclim", "bendo", "bherb",
  "bendoclim", "bherbclim", "bendoherb",
  "bendoherbclim"
)

# Display labels — order here = top-to-bottom on y-axis
coef_labels_ordered <- c(
  "Intercept",
  "Climate",
  "Endophyte",
  "Herbivory",
  "Endophyte \u00d7 Climate",
  "Herbivory \u00d7 Climate",
  "Endophyte \u00d7 Herbivory",
  "Endophyte \u00d7 Herbivory \u00d7 Climate"
)

coef_meta <- tibble(
  parameter = all_coef_names,
  label     = coef_labels_ordered,
  term_type = c(
    "Intercept",
    "Main effect", "Main effect", "Main effect",
    "Two-way interaction", "Two-way interaction", "Two-way interaction",
    "Three-way interaction"
  )
) %>%
  mutate(
    # rev() puts Intercept at TOP of the discrete y-axis
    label     = factor(label, levels = rev(coef_labels_ordered)),
    term_type = factor(
      term_type,
      levels = c("Main effect",
                 "Two-way interaction", "Three-way interaction")
    )
  )

type_colors <- c(
  "Main effect"           = "#4575b4",
  "Two-way interaction"   = "#d95f02",
  "Three-way interaction" = "#7570b3"
)

# ══════════════════════════════════════════════════════════════════════════════
# 3. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

extract_coef_long <- function(posterior_samples, coef_names,
                              species_labels, trait_name) {
  bind_rows(lapply(coef_names, function(coef) {
    mat   <- posterior_samples[[coef]]
    n_spp <- ncol(mat)
    as.data.frame(mat) %>%
      setNames(as.character(seq_len(n_spp))) %>%
      pivot_longer(cols = everything(),
                   names_to  = "species_idx",
                   values_to = "estimate") %>%
      mutate(
        parameter   = coef,
        species_idx = as.integer(species_idx),
        species     = species_labels[species_idx],
        trait       = trait_name
      )
  }))
}

summarise_coef <- function(long_df, ci_width = 0.9) {
  lo <- (1 - ci_width) / 2
  hi <- 1 - lo
  long_df %>%
    group_by(trait, parameter, species) %>%
    summarise(
      median_est = median(estimate),
      lower_CI   = quantile(estimate, lo),
      upper_CI   = quantile(estimate, hi),
      lower_CI95 = quantile(estimate, 0.025),
      upper_CI95 = quantile(estimate, 0.975),
      prob_gt0   = mean(estimate > 0),
      prob_lt0   = mean(estimate < 0),
      .groups    = "drop"
    ) %>%
    mutate(strong_effect = (prob_gt0 > 0.9 | prob_lt0 > 0.9))
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. EXTRACT & SUMMARISE ALL POSTERIORS
# ══════════════════════════════════════════════════════════════════════════════
SPP3 <- c("A. hyemalis", "E. virginicus", "P. autumnalis")
SPP2 <- c("E. virginicus", "P. autumnalis")

all_long <- bind_rows(
  extract_coef_long(posterior_samples_survival, all_coef_names, SPP3, "Survival"),
  extract_coef_long(posterior_samples_grow,     all_coef_names, SPP3, "Growth"),
  extract_coef_long(posterior_samples_inf,      all_coef_names, SPP3, "Inflorescence"),
  extract_coef_long(posterior_samples_spik,     all_coef_names, SPP2, "Spikelet")
)

all_summary <- summarise_coef(all_long) %>%
  left_join(coef_meta, by = "parameter") %>%
  mutate(
    trait   = factor(trait,
                     levels = c("Survival", "Growth", "Inflorescence", "Spikelet")),
    species = factor(species, levels = unique(c(SPP3, SPP2)))
  )

# ══════════════════════════════════════════════════════════════════════════════
# 5. SHARED THEME
# ══════════════════════════════════════════════════════════════════════════════
caterpillar_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position    = "none",          # legend injected separately
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing      = unit(0.5, "lines"),
    axis.text.y        = element_text(size = 8.5, margin = margin(r = 4)),
    axis.text.x        = element_text(size = 7.5),
    axis.title.x       = element_text(size = 9),
    strip.text.x       = element_text(size = 9.5, face = "bold.italic"),
    strip.background   = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    panel.border       = element_rect(color = "black", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.margin        = margin(t = 6, r = 8, b = 4, l = 8, unit = "pt"),
    text               = element_text(family = "Arial")
  )

# ══════════════════════════════════════════════════════════════════════════════
# 6. SHARED LEGEND (built once from a dummy plot, extracted with get_legend)
#
#    Two rows of keys:
#      Row 1 – Term type  (4 filled circles, one per term_type colour)
#      Row 2 – Support    (1 filled + 1 open circle)
#    Approach: dummy ggplot with both aesthetics fully specified so
#    override.aes lengths are always exact — no facet-dependent mismatch.
# ══════════════════════════════════════════════════════════════════════════════
legend_data <- bind_rows(
  # Term-type keys (3 rows — Intercept has its own separate plots)
  tibble(
    x         = 0,
    y         = 1:3,
    col_key   = names(type_colors),
    fill_key  = names(type_colors),
    grp       = "Term type",
    lbl       = c("Main effect", "Two-way interaction", "Three-way interaction")
  ),
  # Posterior-support keys (2 rows)
  tibble(
    x         = 0,
    y         = 4:5,
    col_key   = c("grey40", "grey40"),
    fill_key  = c("filled_strong", "open_uncertain"),
    grp       = "Posterior support",
    lbl       = c("Strong  (P > 90 % on one side of 0)",
                  "Uncertain  (CI crosses 0)")
  )
) %>%
  mutate(
    col_key  = factor(col_key,  levels = c(names(type_colors), "grey40")),
    fill_key = factor(fill_key, levels = c(names(type_colors),
                                           "filled_strong", "open_uncertain"))
  )

dummy_legend_plot <- ggplot(legend_data,
                            aes(x = x, y = lbl,
                                color = col_key,
                                fill  = fill_key)) +
  geom_point(shape = 21, size = 3.5, stroke = 0.9) +
  scale_color_manual(
    name   = "Term type",
    values = c(type_colors, "grey40" = "grey40"),
    labels = c("Main effect", "Two-way interaction", "Three-way interaction",
               "grey40" = ""),
    guide  = "none"
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      type_colors,
      "filled_strong"  = "grey40",
      "open_uncertain" = "white"
    ),
    labels = c(
      "Main effect"           = "Main effect",
      "Two-way interaction"   = "Two-way interaction",
      "Three-way interaction" = "Three-way interaction",
      "filled_strong"         = "Strong  (P > 90 % on one side of 0)",
      "open_uncertain"        = "Uncertain  (CI crosses 0)"
    ),
    guide = guide_legend(
      title        = NULL,
      nrow         = 2,
      override.aes = list(
        color  = c(unname(type_colors), "grey40", "grey40"),
        fill   = c(unname(type_colors), "grey40", "white"),
        shape  = rep(21, 5),
        size   = rep(3.5, 5),
        stroke = rep(0.9, 5)
      )
    )
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.x = unit(0.3, "cm")
  )

shared_legend <- ggpubr::get_legend(dummy_legend_plot)

# ══════════════════════════════════════════════════════════════════════════════
# 7. PLOT-BUILDER FUNCTIONS
#    make_coef_plot()      – all coefficients except Intercept
#    make_intercept_plot() – Intercept (b0) only
# ══════════════════════════════════════════════════════════════════════════════
make_coef_plot <- function(trait_name, x_label, x_limits = NULL) {

  df <- all_summary %>%
    filter(trait == trait_name, parameter != "b0") %>%
    mutate(
      label   = droplevels(label),
      pt_fill = if_else(strong_effect, as.character(term_type), "white")
    )

  spp_in_trait <- levels(droplevels(df$species))
  spp_labeller <- setNames(
    paste0("italic('", spp_in_trait, "')"),
    spp_in_trait
  )

  coef_levels <- levels(df$label)
  n_coef      <- length(coef_levels)

  shade_rects <- lapply(seq_len(n_coef), function(i) {
    if (i %% 2 == 1)
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = i - 0.49, ymax = i + 0.49,
               fill = "grey93", color = NA, alpha = 0.55)
  })
  shade_rects <- shade_rects[!sapply(shade_rects, is.null)]

  present_types <- as.character(unique(df$term_type))
  fill_vals     <- c(type_colors[present_types], "white" = "white")

  p <- ggplot(df, aes(y = label, color = term_type)) +
    shade_rects +
    geom_errorbar(aes(xmin = lower_CI95, xmax = upper_CI95),
                  linewidth = 0.35, width = 0) +
    geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI),
                  linewidth = 1.2, width = 0) +
    geom_point(aes(x = median_est, fill = pt_fill),
               shape = 21, size = 3.0, stroke = 0.9) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_color_manual(values = type_colors, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey25", linewidth = 0.45) +
    facet_grid(. ~ species, scales = "free_x",
               labeller = labeller(
                 species = as_labeller(spp_labeller, label_parsed))) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(x = x_label, y = NULL, title = trait_name) +
    caterpillar_theme

  if (!is.null(x_limits)) p <- p + xlim(x_limits)
  p
}

make_intercept_plot <- function(trait_name, x_label) {

  df <- all_summary %>%
    filter(trait == trait_name, parameter == "b0") %>%
    mutate(
      label   = droplevels(label),
      pt_fill = if_else(strong_effect, as.character(term_type), "white")
    )

  spp_in_trait <- levels(droplevels(df$species))
  spp_labeller <- setNames(
    paste0("italic('", spp_in_trait, "')"),
    spp_in_trait
  )

  fill_vals <- c("Intercept" = "grey40", "white" = "white")

  ggplot(df, aes(y = label, color = term_type)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 0.51, ymax = 1.49,
             fill = "grey93", color = NA, alpha = 0.55) +
    geom_errorbar(aes(xmin = lower_CI95, xmax = upper_CI95),
                  linewidth = 0.35, width = 0) +
    geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI),
                  linewidth = 1.2, width = 0) +
    geom_point(aes(x = median_est, fill = pt_fill),
               shape = 21, size = 3.0, stroke = 0.9) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_color_manual(values = c("Intercept" = "grey40"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey25", linewidth = 0.45) +
    facet_grid(. ~ species, scales = "free_x",
               labeller = labeller(
                 species = as_labeller(spp_labeller, label_parsed))) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(x = x_label, y = NULL,
         title = paste0(trait_name, " — Intercept")) +
    caterpillar_theme
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. BUILD ALL PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Coefficient plots (no intercept) ─────────────────────────────────────────
Fig_surv <- make_coef_plot("Survival",      "Coefficient (logit scale)")
Fig_grow <- make_coef_plot("Growth",        "Coefficient (log scale)")
Fig_inf  <- make_coef_plot("Inflorescence", "Coefficient (log scale)")
Fig_spik <- make_coef_plot("Spikelet",      "Coefficient (log scale)")

# ── Intercept-only plots ──────────────────────────────────────────────────────
Fig_surv_int <- make_intercept_plot("Survival",      "Coefficient (logit scale)")
Fig_grow_int <- make_intercept_plot("Growth",        "Coefficient (log scale)")
Fig_inf_int  <- make_intercept_plot("Inflorescence", "Coefficient (log scale)")
Fig_spik_int <- make_intercept_plot("Spikelet",      "Coefficient (log scale)")

# ══════════════════════════════════════════════════════════════════════════════
# 9. EXPORT
# ══════════════════════════════════════════════════════════════════════════════
out_dir <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure"

# Helper: attach shared legend then write PDF
export_with_legend <- function(plot_obj, filename, width = 8, height = 6.5) {
  p_with_legend <- ggarrange(
    plot_obj,
    ggpubr::as_ggplot(shared_legend),
    ncol    = 1,
    heights = c(1, 0.12)
  )
  Cairo::CairoPDF(file.path(out_dir, filename), width = width, height = height)
  print(p_with_legend)
  dev.off()
  message("Saved: ", filename)
}

# ── Individual PDFs — coefficients ────────────────────────────────────────────
export_with_legend(Fig_surv, "FigS_surv_coef_caterpillar.pdf",  width = 8, height = 6.5)
export_with_legend(Fig_grow, "FigS_grow_coef_caterpillar.pdf",  width = 8, height = 6.5)
export_with_legend(Fig_inf,  "FigS_inf_coef_caterpillar.pdf",   width = 8, height = 6.5)
export_with_legend(Fig_spik, "FigS_spik_coef_caterpillar.pdf",  width = 6, height = 6.5)

# ── Individual PDFs — intercept only ─────────────────────────────────────────
export_with_legend(Fig_surv_int, "FigS_surv_int_caterpillar.pdf",  width = 8, height = 3)
export_with_legend(Fig_grow_int, "FigS_grow_int_caterpillar.pdf",  width = 8, height = 3)
export_with_legend(Fig_inf_int,  "FigS_inf_int_caterpillar.pdf",   width = 8, height = 3)
export_with_legend(Fig_spik_int, "FigS_spik_int_caterpillar.pdf",  width = 6, height = 3)

# ── Combined 2×2 — coefficients ───────────────────────────────────────────────
coef_body <- ggarrange(
  Fig_surv, Fig_grow,
  Fig_inf,  Fig_spik,
  ncol = 2, nrow = 2,
  labels     = c("(a)", "(b)", "(c)", "(d)"),
  font.label = list(size = 10, face = "plain")
)
FigS_coef_combined <- ggarrange(
  coef_body,
  ggpubr::as_ggplot(shared_legend),
  ncol = 1, heights = c(1, 0.07)
)
Cairo::CairoPDF(
  file.path(out_dir, "FigS_coef_caterpillar_combined.pdf"),
  width = 14, height = 11
)
print(FigS_coef_combined)
dev.off()
message("Saved: FigS_coef_caterpillar_combined.pdf")

# ── Combined 1×4 — intercept only ─────────────────────────────────────────────
int_body <- ggarrange(
  Fig_surv_int, Fig_grow_int,
  Fig_inf_int,  Fig_spik_int,
  ncol = 2, nrow = 2,
  labels     = c("(a)", "(b)", "(c)", "(d)"),
  font.label = list(size = 10, face = "plain")
)
FigS_int_combined <- ggarrange(
  int_body,
  ggpubr::as_ggplot(shared_legend),
  ncol = 1, heights = c(1, 0.15)
)
Cairo::CairoPDF(
  file.path(out_dir, "FigS_int_caterpillar_combined.pdf"),
  width = 14, height = 6
)
print(FigS_int_combined)
dev.off()
message("Saved: FigS_int_caterpillar_combined.pdf")
