# Purpose: Plot vital rate models (survival, growth, flowering and spikelet) as
#          function of climate or distance from niche center.
# Authors: Jacob Moutouama
# Date last modified (Y-M-D):

rm(list = ls())

# ── Packages ──────────────────────────────────────────────────────────────────
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(rmutil)
library(actuar)
library(LaplacesDemon)
library(ggpubr)
library(raster)
library(readxl)
library(ggsci)
library(BiocManager)
library(swfscMisc)
library(bayesplot)
library(extraDistr)
library(lubridate)
library(patchwork)
library(ggh4x)
library(xtable)
library(ggtext)

# ── Register fonts with Cairo ──────────────────────────────────────────────────
# Without this, Cairo::CairoPDF() falls back to a substitute font whose metrics
# don't match "Arial" — this breaks vertical offset for plotmath superscripts
# (e.g. italic(S)^{"+"}), causing them to collapse to baseline in the saved PDF
# even though they render correctly in RStudio's on-screen graphics device.
Cairo::CairoFonts(
  regular    = "Arial:style=Regular",
  bold       = "Arial:style=Bold",
  italic     = "Arial:style=Italic",
  bolditalic = "Arial:style=Bold Italic"
)

# ── Shared output path for manuscript figures ──────────────────────────────────
# Set this once here; every figure below is saved into it via file.path(FIG_DIR, ...)
# so you never have to retype the full path again.
FIG_DIR <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Manuscript/Ecology letters/Manuscript/EcoletterR1"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

set.seed(13)
# ── Shared color scheme (use consistently across ALL figures) ─────────────────
# S- = tomato, S+ = cornflowerblue
ENDO_COLORS <- c("0" = "#E2492D", "1" = "#3E6DC9")
ENDO_LABELS <- c("0" = "*S*<sup>\u2212</sup>", "1" = "*S*<sup>+</sup>")

# Display labels for the delta ("Δ (S+ - S-)") panel/facet strip text.
# The underlying data values are left untouched (still "Δ (S+ - S-)") so all
# existing filter()/subset()/case_when() comparisons keep working; only the
# on-plot label changes.
PANEL_LABELS <- c(
  "Pr (survival)"    = "paste(\"Pr (survival)\")",
  "Growth"           = "paste(\"Growth\")",
  "Inflorescences"   = "paste(\"Inflorescences\")",
  "Spikelets"        = "paste(\"Spikelets\")",
  "\u0394 (S+ - S-)" = "paste(Delta, \" (\", italic(S)^{\"+\"} - italic(S)^{\"\u2212\"}, \")\")"
)

# ── Shared theme ──────────────────────────────────────────────────────────────
vr_theme <- function() {
  theme_classic() +
    theme(
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line         = element_line(color = "black", linewidth = 0.1),
      axis.title.x      = element_text(size = 12),
      axis.title.y      = element_text(size = 12),
      axis.text         = element_text(size = 8),
      axis.ticks.x      = element_line(color = "black", linewidth = 0.2),
      axis.ticks.y      = element_line(color = "black", linewidth = 0.2),
      legend.title      = element_text(size = 8),
      legend.text       = element_markdown(size = 10),
      panel.spacing.y   = unit(0.2, "cm"),
      text              = element_text(family = "Arial"),
      strip.text.x      = element_text(size = 10, color = "black"),
      strip.text.y      = element_text(size = 10, color = "black"),
      strip.background  = element_rect(color = "black", fill = "grey90", linewidth = 0.2)
    )
}

# ── Data ──────────────────────────────────────────────────────────────────────
demography_climate <- readRDS(url("https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1"))
climate_max        <- readRDS(url("https://www.dropbox.com/scl/fi/h8m15t64sxkkghgsnsypp/prism_edge_yr_means.rds?rlkey=ljjf16kqvvoe15qmhrqh8qfah&dl=1"))
climate_scaled     <- readRDS(url("https://www.dropbox.com/scl/fi/irecsnoh3xrq6g8d5cysa/climate_site_scaled.rds?rlkey=63r7ugrtkuo5ncmqfoywbidps&dl=1"))

ppt_mean <- mean(climate_scaled$ppt_log)
ppt_sd   <- sd(climate_scaled$ppt_log)

# ══════════════════════════════════════════════════════════════════════════════
# SURVIVAL
# ══════════════════════════════════════════════════════════════════════════════
demography_climate_surv <- demography_climate %>%
  filter(tiller_t > 0) %>%
  dplyr::select(
    Species, Population, Site, site_species_plot, site_year,
    Endo, Herbivory, tiller_t, surv1, cum_ppt, ppt_scaled
  ) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_year         = as.integer(factor(site_year)),
    site_species_plot = as.integer(factor(site_species_plot)),
    log_size_t0       = log(tiller_t),
    surv_t1           = as.integer(surv1),
    ppt               = ppt_scaled
  )

demography_surv_ppt <- list(
  nSpp       = n_distinct(demography_climate_surv$Species),
  nSite      = n_distinct(demography_climate_surv$Site),
  nsite_year = n_distinct(demography_climate_surv$site_year),
  nPop       = n_distinct(demography_climate_surv$Population),
  nPlot      = n_distinct(demography_climate_surv$site_species_plot),
  Spp        = demography_climate_surv$Species,
  site       = demography_climate_surv$Site,
  site_year  = demography_climate_surv$site_year,
  pop        = demography_climate_surv$Population,
  plot       = demography_climate_surv$site_species_plot,
  clim       = as.vector(demography_climate_surv$ppt),
  endo       = demography_climate_surv$Endo,
  herb       = demography_climate_surv$Herbivory,
  size       = demography_climate_surv$log_size_t0,
  y          = demography_climate_surv$surv_t1,
  N          = nrow(demography_climate_surv)
)

fit_surv_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1"))

# Species-specific prediction ranges
climate_range_per_species_surv <- lapply(1:3, function(sp) {
  sp_clim <- demography_surv_ppt$clim[demography_surv_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})

predictions <- do.call(rbind, lapply(1:3, function(sp) {
  expand.grid(
    clim    = climate_range_per_species_surv[[sp]],
    endo    = c(0, 1),
    herb    = c(0, 1),
    species = sp
  )
}))

posterior_samples_survival <- rstan::extract(fit_surv_ppt)

get_predictions_survival <- function(clim, endo, herb, species_index, ps) {
  with(ps, {
    logit_preds <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index]    * endo * clim +
      bendoherb[, species_index]    * endo * herb +
      bherbclim[, species_index]    * herb * clim +
      bendoherbclim[, species_index] * endo * herb * clim
    plogis(logit_preds)
  })
}

n_post_survival <- nrow(posterior_samples_survival$b0)
pred_probs_matrix_survival <- matrix(NA, nrow = nrow(predictions), ncol = n_post_survival)
for (i in seq_len(nrow(predictions))) {
  pred_probs_matrix_survival[i, ] <- get_predictions_survival(
    predictions$clim[i], predictions$endo[i],
    predictions$herb[i], predictions$species[i],
    posterior_samples_survival
  )
}

pred_probs_df_survival  <- cbind(predictions, as.data.frame(pred_probs_matrix_survival))
pred_probs_long_df_survival <- pred_probs_df_survival %>%
  pivot_longer(cols = starts_with("V"), names_to = "Posterior_Sample", values_to = "Pred_Survival")

cred_intervals_survival <- pred_probs_long_df_survival %>%
  group_by(species, endo, herb, clim) %>%
  summarise(
    lower_90 = quantile(Pred_Survival, 0.05),
    upper_90 = quantile(Pred_Survival, 0.95),
    median   = quantile(Pred_Survival, 0.5),
    mean     = mean(Pred_Survival),
    .groups  = "drop"
  ) %>%
  mutate(panel = "Pr (survival)")

# ── Observed: plot means only (SD removed — cleaner with model ribbon present) ──
observed_data_survival <- demography_surv_ppt %>%
  data.frame(
    clim    = .$clim,
    endo    = .$endo,
    herb    = .$herb,
    species = .$Spp,
    plot    = .$plot,
    y       = .$y
  ) %>%
  group_by(plot, species, herb, clim, endo) %>%
  summarise(
    y_plot_mean = mean(y, na.rm = TRUE),
    n_obs       = sum(!is.na(y)),
    .groups     = "drop"
  ) %>%
  mutate(
    panel      = "Pr (survival)",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

diff_df_survival <- pred_probs_long_df_survival %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(
    diff = mean(Pred_Survival[endo == 1]) - mean(Pred_Survival[endo == 0]),
    .groups = "drop"
  )

diff_ci_survival <- diff_df_survival %>%
  group_by(species, herb, clim) %>%
  summarise(
    lower_90 = quantile(diff, 0.05),
    upper_90 = quantile(diff, 0.95),
    mean     = mean(diff),
    .groups  = "drop"
  ) %>%
  mutate(panel = "Δ (S+ - S-)")

plot_data_survival <- bind_rows(cred_intervals_survival, diff_ci_survival)

species_levels_3 <- c(
  "italic('Agrostis hyemalis')",
  "italic('Elymus virginicus')",
  "italic('Poa autumnalis')"
)

plot_data_survival$species <- factor(plot_data_survival$species, levels = 1:3, labels = species_levels_3)
observed_data_survival$species <- factor(observed_data_survival$species, levels = 1:3, labels = species_levels_3)
plot_data_survival <- plot_data_survival %>% mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))
plot_data_survival$panel <- trimws(as.character(plot_data_survival$panel))
observed_data_survival$panel <- trimws(as.character(observed_data_survival$panel))

# Pre-compute jitter offset so errorbar and point share the exact same x position
observed_data_survival <- observed_data_survival %>%
  mutate(jitter_x = climate_mm + runif(n(), -18, 18))

panel_labels_surv <- data.frame(
  species = rep(species_levels_3, each = 2),
  herb    = rep(c(0, 1), times = 3),
  label   = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  panel   = "Pr (survival)"
)
# Get n_obs range for caption
range(observed_data_survival$n_obs)

# ── Create visible dash segments for the zero line ───────────────────────────
make_dashes <- function(xmin, xmax, n = 16, dash_prop = 0.65) {
  step <- (xmax - xmin) / n
  
  tibble(
    x    = xmin + (0:(n - 1)) * step,
    xend = xmin + (0:(n - 1)) * step + step * dash_prop,
    y    = 0,
    yend = 0
  )
}

global_xmin <- min(plot_data_survival$climate_mm, na.rm = TRUE)
global_xmax <- max(plot_data_survival$climate_mm, na.rm = TRUE)

zero_dashes <- plot_data_survival %>%
  filter(panel == "Δ (S+ - S-)") %>%
  distinct(species, herb, panel) %>%
  mutate(
    dashes = purrr::map(
      seq_len(n()),
      ~ make_dashes(global_xmin, global_xmax)
    )
  ) %>%
  tidyr::unnest(dashes)

Cairo::CairoTIFF(
  file.path(FIG_DIR, "PrSurvival_diff.tiff"),
  width = 6,
  height = 7,
  units = "in",
  dpi = 600
)

ggplot(plot_data_survival) +
  # ── Predicted survival ──
  geom_line(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, y = mean, color = factor(endo), group = endo)
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Pr (survival)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90,
        fill = factor(endo), group = endo),
    alpha = 0.2, color = NA
  ) +
  # ── Observed: plot means only, size ∝ n_obs ──────────────────────────────
  geom_point(
    data = subset(observed_data_survival, panel == "Pr (survival)"),
    aes(
      x     = jitter_x,
      y     = y_plot_mean,
      color = factor(endo),
      size  = n_obs
    ),
    alpha       = 0.5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(0.5, 4)) +   # wide range: small n vs large n obvious
  # ── Δ panel ──
  geom_line(
    data = subset(plot_data_survival, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, y = mean),
    color = "black", linewidth = 0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_survival, panel == "Δ (S+ - S-)"),
    aes(x = climate_mm, ymin = lower_90, ymax = upper_90),
    fill = "#D9BFD6", alpha = 0.6
  ) +
  geom_segment(
    data = zero_dashes,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45,
    color = "black",
    lineend = "butt",
    inherit.aes = FALSE
  ) +

  # ── Facets ──
  ggh4x::facet_nested(
    species + panel ~ herb,
    scales = "free", space = "fixed",
    labeller = labeller(
      species = label_parsed,
      herb    = c("0" = "Herbivory access", "1" = "Herbivory exclusion"),
      panel   = ggplot2::as_labeller(PANEL_LABELS, default = label_parsed)
    )
  ) +
  ggh4x::facetted_pos_scales(
    x = list(
      panel == "Δ (S+ - S-)" ~ scale_x_continuous(
        expand = c(0, 0)
      )
    ),
    y = list(
      panel == "Δ (S+ - S-)" ~ scale_y_continuous(
        limits = c(-0.3, 0.35),
        expand = c(0, 0),
        breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
        labels = c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3"),
        minor_breaks = NULL
      ),
      panel == "Pr (survival)" ~ scale_y_continuous(
        limits = c(-0.05, 1.05),
        expand = c(0, 0),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      )
    )
  ) +
  labs(
    x = "Precipitation (mm)",
    y = expression(paste("Survival probability / ", Delta, " survival (",
                          italic(S)^{"+"} - italic(S)^{"\u2212"}, ")")),
    color = "Symbiont", fill = "Symbiont"
  ) +
  scale_color_manual(values = ENDO_COLORS, labels = ENDO_LABELS) +
  scale_fill_manual(values  = ENDO_COLORS, labels = ENDO_LABELS) +
  vr_theme() +
  theme(
    legend.position = c(0.18, 0.461),
    legend.direction = "horizontal",
    legend.justification = "center",
    axis.title.y = element_text(size = 8)
  )+
  geom_text(
    data = panel_labels_surv,
    aes(x = -Inf, y = Inf, label = label),
    fontface = "plain",
    size = 4,
    hjust = -0.3,
    vjust = 1.3,
    inherit.aes = FALSE
  )
dev.off()

# ── Delta survival summary (used downstream) ──────────────────────────────────
compute_delta_surv <- function(clim_val, ps, herb_values = c(0, 1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      pred_Splus <- plogis(
        ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
          ps$bendoclim[,sp]*1*clim_val + ps$bendoherb[,sp]*1*h +
          ps$bherbclim[,sp]*h*clim_val + ps$bendoherbclim[,sp]*1*h*clim_val
      )
      pred_Sminus <- plogis(
        ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
          ps$bendoclim[,sp]*0*clim_val + ps$bendoherb[,sp]*0*h +
          ps$bherbclim[,sp]*h*clim_val + ps$bendoherbclim[,sp]*0*h*clim_val
      )
      data.frame(species=sp, herb=h, clim=clim_val,
                 Posterior_Sample=1:n_post, delta=pred_Splus - pred_Sminus)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_surv_ppt$clim[demography_surv_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})
names(climate_range_per_species) <- 1:3

delta_surv_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_surv(cl, posterior_samples_survival) %>% filter(species == sp)
  }) %>% bind_rows()
}) %>% bind_rows()

delta_surv_summary <- delta_surv_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(
    median_delta   = median(delta),
    lower_90       = quantile(delta, 0.05),
    upper_90       = quantile(delta, 0.95),
    prob_delta_gt0 = mean(delta > 0),
    .groups = "drop"
  ) %>%
  mutate(
    species  = factor(species, levels=1:3,
                      labels=c("Agrostis hyemalis","Elymus virginicus","Poa autumnalis")),
    herb     = factor(herb, levels=c(0,1),
                      labels=c("Herbivory access","Herbivory exclusion")),
    clim_mm  = exp(clim * ppt_sd + ppt_mean)
  )

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()
for(sp in unique(delta_surv_summary$species)) {
  for(h in unique(delta_surv_summary$herb)) {
    df_sub <- delta_surv_summary %>%
      filter(species == sp, herb == h)
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}
delta_surv_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm        = round(clim_mm, 0),
    median_delta   = round(median_delta, 3),
    lower_90       = round(lower_90, 3),
    upper_90       = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

delta_long_surv <- delta_surv_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(cols = c(median_delta, prob_delta_gt0),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    "median_delta"   = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "*Agrostis hyemalis*",
      species == "Elymus virginicus" ~ "*Elymus virginicus*",
      species == "Poa autumnalis"    ~ "*Poa autumnalis*"
    )
  )

# ── Abbreviation map (shared across all vital rates) ─────────────────────────
species_abbrev <- c(
  "Agrostis hyemalis" = "A. hyemalis",
  "Elymus virginicus" = "E. virginicus",
  "Poa autumnalis"    = "P. autumnalis"
)

# ── LaTeX table: survival ─────────────────────────────────────────────────────
delta_surv_latex <- delta_surv_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm  = clim_mm,
    Median_delta      = median_delta,
    Lower_90          = lower_90,
    Upper_90          = upper_90,
    Prob_delta_gt0    = prob_delta_gt0
  )

delta_surv_xt <- xtable(
  delta_surv_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)
colnames(delta_surv_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)
print(
  delta_surv_xt,
  include.rownames       = FALSE,
  sanitize.text.function = identity,
  floating               = FALSE
)

# ══════════════════════════════════════════════════════════════════════════════
# GROWTH
# ══════════════════════════════════════════════════════════════════════════════
demography_climate %>%
  filter(tiller_t > 0 & tiller_t1 > 0) %>%
  dplyr::select(Species, Population, Site, Plot, site_year, site_species_plot,
                Endo, Herbivory, tiller_t, grow, cum_ppt, ppt_scaled) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    site_year         = as.integer(factor(site_year)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  ) -> demography_climate_grow

demography_grow_ppt <- list(
  nSpp       = n_distinct(demography_climate_grow$Species),
  nSite      = n_distinct(demography_climate_grow$Site),
  nsite_year = n_distinct(demography_climate_grow$site_year),
  nPop       = n_distinct(demography_climate_grow$Population),
  nPlot      = n_distinct(demography_climate_grow$site_species_plot),
  Spp        = demography_climate_grow$Species,
  site       = demography_climate_grow$Site,
  site_year  = demography_climate_grow$site_year,
  pop        = demography_climate_grow$Population,
  plot       = demography_climate_grow$site_species_plot,
  clim       = as.vector(demography_climate_grow$ppt),
  endo       = demography_climate_grow$Endo,
  herb       = demography_climate_grow$Herbivory,
  size       = demography_climate_grow$log_size_t0,
  y          = demography_climate_grow$grow,
  N          = nrow(demography_climate_grow)
)

fit_grow_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1"))

climate_range_per_species_grow <- lapply(1:3, function(sp) {
  sp_clim <- demography_grow_ppt$clim[demography_grow_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out = 30)
})

predictions <- do.call(rbind, lapply(1:3, function(sp) {
  expand.grid(clim=climate_range_per_species_grow[[sp]],
              endo=c(0,1), herb=c(0,1), species=sp)
}))

posterior_samples_grow <- rstan::extract(fit_grow_ppt)

get_predictions_grow <- function(clim, endo, herb, species_index, ps) {
  with(ps, {
    b0[,species_index] + bendo[,species_index]*endo + bherb[,species_index]*herb +
      bclim[,species_index]*clim + bendoclim[,species_index]*clim*endo +
      bendoherb[,species_index]*endo*herb + bherbclim[,species_index]*herb*clim +
      bendoherbclim[,species_index]*endo*herb*clim
  })
}

n_post_grow <- nrow(posterior_samples_grow$b0)
pred_matrix_grow <- matrix(NA, nrow=nrow(predictions), ncol=n_post_grow)
for (i in seq_len(nrow(predictions))) {
  pred_matrix_grow[i,] <- get_predictions_grow(
    predictions$clim[i], predictions$endo[i],
    predictions$herb[i], predictions$species[i], posterior_samples_grow
  )
}

pred_grow_long <- cbind(predictions, as.data.frame(pred_matrix_grow)) %>%
  pivot_longer(cols=starts_with("V"), names_to="Posterior_Sample", values_to="Pred_Growth")

cred_intervals_grow <- pred_grow_long %>%
  group_by(species, endo, herb, clim) %>%
  summarise(lower_90=quantile(Pred_Growth,0.05), upper_90=quantile(Pred_Growth,0.95),
            median=quantile(Pred_Growth,0.5), mean=mean(Pred_Growth), .groups="drop") %>%
  mutate(panel = "Growth")

diff_ci_grow <- pred_grow_long %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(diff = mean(Pred_Growth[endo==1]) - mean(Pred_Growth[endo==0]), .groups="drop") %>%
  group_by(species, herb, clim) %>%
  summarise(lower_90=quantile(diff,0.05), upper_90=quantile(diff,0.95),
            mean=mean(diff), .groups="drop") %>%
  mutate(panel = "Δ (S+ - S-)")

plot_data_grow <- bind_rows(cred_intervals_grow, diff_ci_grow)
plot_data_grow$species <- factor(plot_data_grow$species, levels=1:3, labels=species_levels_3)

# ── Observed growth: plot means only ─────────────────────────────────────────
observed_data_grow <- demography_grow_ppt %>%
  data.frame(clim=.$clim, endo=.$endo, herb=.$herb,
             species=.$Spp, plot=.$plot, y=.$y) %>%
  group_by(plot, species, herb, clim, endo) %>%
  summarise(
    y_plot_mean = mean(y, na.rm=TRUE),
    n_obs       = sum(!is.na(y)),
    .groups     = "drop"
  ) %>%
  mutate(
    panel      = "Growth",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

observed_data_grow$species <- factor(observed_data_grow$species, levels=1:3, labels=species_levels_3)
plot_data_grow <- plot_data_grow %>% mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# Pre-compute jitter offset
observed_data_grow <- observed_data_grow %>%
  mutate(jitter_x = climate_mm + runif(n(), -20, 20))

y_limits <- plot_data_grow %>%
  filter(panel == "Δ (S+ - S-)") %>%
  group_by(species) %>%
  summarise(ymin=min(lower_90, na.rm=TRUE), ymax=max(upper_90, na.rm=TRUE))

panel_labels_grow <- data.frame(
  species = rep(species_levels_3, each=2),
  herb    = rep(c(0,1), times=3),
  label   = c("(a)","(b)","(c)","(d)","(e)","(f)"),
  panel   = "Growth",
  ymax    = c(1, 1, 1.25, 1.25, 2.2, 2.2)
)

# Get n_obs range for caption
range(observed_data_grow$n_obs)

# ── Zero-line dashes (same style as survival, for harmony across figures) ─────
global_xmin_grow <- min(plot_data_grow$climate_mm, na.rm = TRUE)
global_xmax_grow <- max(plot_data_grow$climate_mm, na.rm = TRUE)

zero_dashes_grow <- plot_data_grow %>%
  filter(panel == "Δ (S+ - S-)") %>%
  distinct(species, herb, panel) %>%
  mutate(
    dashes = purrr::map(
      seq_len(n()),
      ~ make_dashes(global_xmin_grow, global_xmax_grow)
    )
  ) %>%
  tidyr::unnest(dashes)

Cairo::CairoTIFF(
  file.path(FIG_DIR, "Growth_diff.tiff"),
  width = 7, height = 8, units = "in", dpi = 600
)
ggplot(plot_data_grow) +
  geom_line(
    data = subset(plot_data_grow, panel=="Growth"),
    aes(x=climate_mm, y=mean, color=factor(endo), group=endo), linewidth=0.5
  ) +
  geom_ribbon(
    data = subset(plot_data_grow, panel=="Growth"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90, fill=factor(endo), group=endo),
    alpha=0.2, color=NA
  ) +
  # ── Observed: plot means only, size ∝ n_obs ──────────────────────────────
  geom_point(
    data = subset(observed_data_grow, panel=="Growth"),
    aes(x=jitter_x, y=y_plot_mean, color=factor(endo), size=n_obs),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(0.7, 5.3))+   # was c(0.5, 4)
  geom_line(
    data=subset(plot_data_grow, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, y=mean), color="black", linewidth=0.5
  ) +
  geom_ribbon(
    data=subset(plot_data_grow, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90),
    fill="#D9BFD6", alpha=0.6
  ) +
  geom_segment(
    data = zero_dashes_grow,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45,
    color = "black",
    lineend = "butt",
    inherit.aes = FALSE
  ) +
  ggh4x::facet_nested(
    species + panel ~ herb, scales="free",
    labeller=labeller(species=label_parsed,
                      herb=c("0"="Herbivory access","1"="Herbivory exclusion"),
                      panel=ggplot2::as_labeller(PANEL_LABELS, default = label_parsed))
  ) +
  ggh4x::facetted_pos_scales(y=list(
    panel=="Δ (S+ - S-)" & species=="italic('Agrostis hyemalis')" ~
      scale_y_continuous(minor_breaks=NULL,
        limits=c(y_limits$ymin[y_limits$species=="italic('Agrostis hyemalis')"],
                 y_limits$ymax[y_limits$species=="italic('Agrostis hyemalis')"]),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Δ (S+ - S-)" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(minor_breaks=NULL,
        limits=c(y_limits$ymin[y_limits$species=="italic('Elymus virginicus')"],
                 y_limits$ymax[y_limits$species=="italic('Elymus virginicus')"]),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Δ (S+ - S-)" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(minor_breaks=NULL,
        limits=c(y_limits$ymin[y_limits$species=="italic('Poa autumnalis')"],
                 y_limits$ymax[y_limits$species=="italic('Poa autumnalis')"]),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Growth" & species=="italic('Agrostis hyemalis')" ~
      scale_y_continuous(limits=c(-2.5, 1),  expand=c(0,0)),
    panel=="Growth" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(limits=c(-1.1, 1.25), expand=c(0,0)),
    panel=="Growth" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(limits=c(-2.8, 2.2), expand=c(0,0))
  )) +
  labs(x="Precipitation (mm)",
       y = expression(paste("Log size ratio / ", Delta, " growth (",
                             italic(S)^{"+"} - italic(S)^{"\u2212"}, ")")),
       color="Symbiont", fill="Symbiont") +
  scale_color_manual(values=ENDO_COLORS, labels=ENDO_LABELS) +
  scale_fill_manual(values=ENDO_COLORS,  labels=ENDO_LABELS) +
  vr_theme() +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.27),
    axis.line         = element_line(color = "black", linewidth = 0.13),
    axis.title        = element_markdown(size = 11),
    axis.title.y      = element_text(size = 11),
    axis.text         = element_text(size = 8),
    axis.ticks.x      = element_line(color = "black", linewidth = 0.27),
    axis.ticks.y      = element_line(color = "black", linewidth = 0.27),
    legend.title      = element_text(size = 10),
    legend.text       = element_markdown(size = 10),
    panel.spacing.y   = unit(0.2, "cm"),
    text              = element_text(family = "Arial"),
    strip.text.x      = element_text(size = 13, color = "black"),
    strip.text.y      = element_text(size = 11, color = "black"),
    strip.background  = element_rect(color = "black", fill = "grey90", linewidth = 0.27),
    legend.position   = c(0.82, 0.37),
    legend.direction  = "horizontal",
    legend.justification = "center"
  ) +
  geom_text(data = panel_labels_grow,
            aes(x = 490, y = ymax * 0.70, label = label),
            hjust = 0, size = 4.7, inherit.aes = FALSE)   # was size = 3.5
dev.off()

# ── Delta growth summary ───────────────────────────────────────────────────────
compute_delta_grow <- function(clim_val, ps, herb_values=c(0,1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      pred_Eplus  <- ps$b0[,sp] + ps$bendo[,sp]*1 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
        ps$bendoclim[,sp]*clim_val*1 + ps$bendoherb[,sp]*1*h + ps$bendoherbclim[,sp]*1*h*clim_val
      pred_Eminus <- ps$b0[,sp] + ps$bendo[,sp]*0 + ps$bherb[,sp]*h + ps$bclim[,sp]*clim_val +
        ps$bendoclim[,sp]*clim_val*0 + ps$bendoherb[,sp]*0*h + ps$bendoherbclim[,sp]*0*h*clim_val
      data.frame(species=sp, herb=h, clim=clim_val,
                 Posterior_Sample=1:n_post, delta=pred_Eplus - pred_Eminus)
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
}

climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_grow_ppt$clim[demography_grow_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out=30)
})
names(climate_range_per_species) <- 1:3

delta_grow_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_grow(cl, posterior_samples_grow) %>% filter(species==sp)
  }) %>% bind_rows()
}) %>% bind_rows()

delta_grow_summary <- delta_grow_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(median_delta=median(delta), lower_90=quantile(delta,0.05),
            upper_90=quantile(delta,0.95), prob_delta_gt0=mean(delta>0), .groups="drop") %>%
  mutate(
    species = factor(species, levels=1:3,
                     labels=c("Agrostis hyemalis","Elymus virginicus","Poa autumnalis")),
    herb    = factor(herb, levels=c(0,1),
                     labels=c("Herbivory access","Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  )

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()
for(sp in unique(delta_grow_summary$species)) {
  for(h in unique(delta_grow_summary$herb)) {
    df_sub <- delta_grow_summary %>%
      filter(species == sp, herb == h)
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}
delta_grow_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm        = round(clim_mm, 0),
    median_delta   = round(median_delta, 3),
    lower_90       = round(lower_90, 3),
    upper_90       = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

delta_long_grow <- delta_grow_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(cols=c(median_delta, prob_delta_gt0), names_to="metric", values_to="value") %>%
  mutate(
    metric = recode(metric,
                    "median_delta"   = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "*Agrostis hyemalis*",
      species == "Elymus virginicus" ~ "*Elymus virginicus*",
      species == "Poa autumnalis"    ~ "*Poa autumnalis*"
    )
  )

# ── LaTeX table: growth ───────────────────────────────────────────────────────
delta_grow_latex <- delta_grow_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm  = clim_mm,
    Median_delta      = median_delta,
    Lower_90          = lower_90,
    Upper_90          = upper_90,
    Prob_delta_gt0    = prob_delta_gt0
  )

delta_grow_xt <- xtable(
  delta_grow_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)
colnames(delta_grow_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)
print(
  delta_grow_xt,
  include.rownames       = FALSE,
  sanitize.text.function = identity,
  floating               = FALSE
)

# ══════════════════════════════════════════════════════════════════════════════
# INFLORESCENCE
# ══════════════════════════════════════════════════════════════════════════════
demography_climate %>%
  filter(tiller_t1 > 0) %>%
  dplyr::select(Species, Population, Site, site_year, Plot, site_species_plot,
                Endo, Herbivory, tiller_t, inf_t1, cum_ppt, ppt_scaled) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    site_year         = as.integer(factor(site_year)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  ) -> demography_climate_inf

demography_inf_ppt <- list(
  nSpp       = n_distinct(demography_climate_inf$Species),
  nSite      = n_distinct(demography_climate_inf$Site),
  nsite_year = n_distinct(demography_climate_inf$site_year),
  nPop       = n_distinct(demography_climate_inf$Population),
  nPlot      = n_distinct(demography_climate_inf$site_species_plot),
  Spp        = demography_climate_inf$Species,
  site       = demography_climate_inf$Site,
  site_year  = demography_climate_inf$site_year,
  pop        = demography_climate_inf$Population,
  plot       = demography_climate_inf$site_species_plot,
  clim       = as.vector(demography_climate_inf$ppt),
  endo       = demography_climate_inf$Endo,
  herb       = demography_climate_inf$Herbivory,
  size       = demography_climate_inf$log_size_t0,
  y          = demography_climate_inf$inf_t1,
  N          = nrow(demography_climate_inf)
)

fit_inf_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/v5kcirie79uae1dex6l1a/fit_inf_abio_bio_endo_hurdle_linear.rds?rlkey=pxp0pdptn95j5ubkjs819ou3x&dl=1"))

climate_range_per_species_inf <- lapply(1:3, function(sp) {
  sp_clim <- demography_inf_ppt$clim[demography_inf_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out=30)
})

predictions <- do.call(rbind, lapply(1:3, function(sp) {
  expand.grid(clim=climate_range_per_species_inf[[sp]],
              endo=c(0,1), herb=c(0,1), species=sp)
}))

posterior_samples_inf <- rstan::extract(fit_inf_ppt)

# get_predictions_inf <- function(clim, endo, herb, species_index, ps) {
#   with(ps, {
#     eta <- b0[,species_index] + bendo[,species_index]*endo + bherb[,species_index]*herb +
#       bclim[,species_index]*clim + bendoclim[,species_index]*clim*endo +
#       bendoherb[,species_index]*endo*herb + bherbclim[,species_index]*herb*clim +
#       bendoherbclim[,species_index]*endo*herb*clim
#     (1 - zi) * exp(eta)
#   })
# }

get_predictions_inf <- function(clim, endo, herb, species_index, ps) {
  with(ps, {
    eta <- b0[, species_index] +
      bendo[, species_index] * endo +
      bherb[, species_index] * herb +
      bclim[, species_index] * clim +
      bendoclim[, species_index] * clim * endo +
      bendoherb[, species_index] * endo * herb +
      bherbclim[, species_index] * herb * clim +
      bendoherbclim[, species_index] * endo * herb * clim
    
    exp(eta)
  })
}

n_post_inf <- nrow(posterior_samples_inf$b0)
pred_matrix_inf <- matrix(NA, nrow=nrow(predictions), ncol=n_post_inf)
for (i in seq_len(nrow(predictions))) {
  pred_matrix_inf[i,] <- get_predictions_inf(
    predictions$clim[i], predictions$endo[i],
    predictions$herb[i], predictions$species[i], posterior_samples_inf
  )
}

pred_inf_long <- cbind(predictions, as.data.frame(pred_matrix_inf)) %>%
  pivot_longer(cols=starts_with("V"), names_to="Posterior_Sample", values_to="Pred_Inf")

cred_intervals_inf <- pred_inf_long %>%
  group_by(species, endo, herb, clim) %>%
  summarise(lower_90=quantile(Pred_Inf,0.05), upper_90=quantile(Pred_Inf,0.95),
            median=quantile(Pred_Inf,0.5), mean=mean(Pred_Inf), .groups="drop") %>%
  mutate(panel="Inflorescences")

diff_ci_inf <- pred_inf_long %>%
  group_by(species, herb, clim, Posterior_Sample) %>%
  summarise(diff=mean(Pred_Inf[endo==1])-mean(Pred_Inf[endo==0]), .groups="drop") %>%
  group_by(species, herb, clim) %>%
  summarise(lower_90=quantile(diff,0.05), upper_90=quantile(diff,0.95),
            mean=mean(diff), .groups="drop") %>%
  mutate(panel="Δ (S+ - S-)")

plot_data_inf <- bind_rows(cred_intervals_inf, diff_ci_inf)
plot_data_inf$species <- factor(plot_data_inf$species, levels=1:3, labels=species_levels_3)

# ── Observed inflorescence: plot means only ───────────────────────────────────
observed_data_inf <- demography_inf_ppt %>%
  data.frame(clim=.$clim, endo=.$endo, herb=.$herb,
             species=.$Spp, plot=.$plot, y=.$y) %>%
  group_by(plot, species, herb, clim, endo) %>%
  summarise(
    y_plot_mean = mean(y, na.rm=TRUE),
    n_obs       = sum(!is.na(y)),
    .groups     = "drop"
  ) %>%
  mutate(
    panel      = "Inflorescences",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

observed_data_inf$species <- factor(observed_data_inf$species, levels=1:3, labels=species_levels_3)
plot_data_inf <- plot_data_inf %>% mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# Pre-compute jitter offset

observed_data_inf <- observed_data_inf %>%
  mutate(jitter_x = climate_mm + runif(n(), -18, 18))

y_limits_inf <- plot_data_inf %>%
  filter(panel=="Δ (S+ - S-)") %>%
  group_by(species) %>%
  summarise(ymin=min(lower_90,na.rm=TRUE), ymax=max(upper_90,na.rm=TRUE), .groups="drop")

panel_labels_inf <- data.frame(
  species = rep(species_levels_3, each=2),
  herb    = rep(c(0,1), times=3),
  label   = c("(a)","(b)","(c)","(d)","(e)","(f)"),
  panel   = "Inflorescences",
  ymax    = rep(c(20, 7, 65), each=2)
)

# ── Zero-line dashes (same style as survival, for harmony across figures) ─────
global_xmin_inf <- min(plot_data_inf$climate_mm, na.rm = TRUE)
global_xmax_inf <- max(plot_data_inf$climate_mm, na.rm = TRUE)

zero_dashes_inf <- plot_data_inf %>%
  filter(panel == "Δ (S+ - S-)") %>%
  distinct(species, herb, panel) %>%
  mutate(
    dashes = purrr::map(
      seq_len(n()),
      ~ make_dashes(global_xmin_inf, global_xmax_inf)
    )
  ) %>%
  tidyr::unnest(dashes)

Cairo::CairoTIFF(
  file.path(FIG_DIR, "Inflorescence_diff_v.tiff"),
  width = 6, height = 7, units = "in", dpi = 600
)
ggplot(plot_data_inf) +
  geom_line(
    data=subset(plot_data_inf, panel=="Inflorescences"),
    aes(x=climate_mm, y=mean, color=factor(endo), group=endo), linewidth=0.5
  ) +
  geom_ribbon(
    data=subset(plot_data_inf, panel=="Inflorescences"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90, fill=factor(endo), group=endo),
    alpha=0.3, color=NA
  ) +
  # ── Observed: plot means only, size ∝ n_obs ──────────────────────────────
  geom_point(
    data=subset(observed_data_inf, panel=="Inflorescences"),
    aes(x=jitter_x, y=y_plot_mean, color=factor(endo), size=n_obs),
    alpha=0.5, show.legend=FALSE
  ) +
  scale_size_continuous(range=c(0.5, 4)) +
  geom_line(
    data=subset(plot_data_inf, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, y=mean), color="black", linewidth=0.5
  ) +
  geom_ribbon(
    data=subset(plot_data_inf, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90),
    fill="#D9BFD6", alpha=0.6
  ) +
  geom_segment(
    data = zero_dashes_inf,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45,
    color = "black",
    lineend = "butt",
    inherit.aes = FALSE
  ) +
  ggh4x::facet_nested(
    species + panel ~ herb, scales="free_y",
    labeller=labeller(species=label_parsed,
                      herb=c("0"="Herbivory access","1"="Herbivory exclusion"),
                      panel=ggplot2::as_labeller(PANEL_LABELS, default = label_parsed))
  ) +
  ggh4x::facetted_pos_scales(y=list(
    panel=="Δ (S+ - S-)" & species=="italic('Agrostis hyemalis')" ~
      scale_y_continuous(minor_breaks=NULL,
        limits=c(y_limits_inf$ymin[y_limits_inf$species=="italic('Agrostis hyemalis')"],
                 y_limits_inf$ymax[y_limits_inf$species=="italic('Agrostis hyemalis')"]),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Δ (S+ - S-)" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(minor_breaks=NULL, limits=c(-1.5,3.2),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Δ (S+ - S-)" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(minor_breaks=NULL,
        limits=c(y_limits_inf$ymin[y_limits_inf$species=="italic('Poa autumnalis')"],
                 y_limits_inf$ymax[y_limits_inf$species=="italic('Poa autumnalis')"]),
        breaks = scales::pretty_breaks(n = 4),
        expand=c(0,0)),
    panel=="Inflorescences" & species=="italic('Agrostis hyemalis')" ~
      scale_y_continuous(limits=c(0,20)),
    panel=="Inflorescences" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(limits=c(0,7)),
    panel=="Inflorescences" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(limits=c(0,65))
  )) +
  labs(x="Precipitation (mm)",
       y = expression(paste("Number of inflorescences / ", Delta, " inflorescences (",
                             italic(S)^{"+"} - italic(S)^{"\u2212"}, ")")),
       color="Symbiont", fill="Symbiont") +
  scale_color_manual(values=ENDO_COLORS, labels=ENDO_LABELS) +
  scale_fill_manual(values=ENDO_COLORS,  labels=ENDO_LABELS) +
  vr_theme() +
  theme(legend.position=c(0.8, 0.46),
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.title.y    = element_text(size = 10)) +
  geom_text(data=panel_labels_inf, aes(x=490, y=ymax*0.8, label=label),
            hjust=0, size=3.5, inherit.aes=FALSE)
dev.off()

# ── Delta inflorescence summary ───────────────────────────────────────────────
compute_delta_inf <- function(clim, ps, herb_values=c(0,1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      eta_p <- ps$b0[,sp]+ps$bendo[,sp]*1+ps$bherb[,sp]*h+ps$bclim[,sp]*clim+
        ps$bendoclim[,sp]*1*clim+ps$bendoherb[,sp]*1*h+
        ps$bherbclim[,sp]*h*clim+ps$bendoherbclim[,sp]*1*h*clim
      eta_m <- ps$b0[,sp]+ps$bendo[,sp]*0+ps$bherb[,sp]*h+ps$bclim[,sp]*clim+
        ps$bendoclim[,sp]*0*clim+ps$bendoherb[,sp]*0*h+
        ps$bherbclim[,sp]*h*clim+ps$bendoherbclim[,sp]*0*h*clim
      data.frame(species=sp, herb=h, clim=clim,
                 Posterior_Sample=1:n_post, delta=exp(eta_p)-exp(eta_m))
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
}

climate_range_per_species <- lapply(1:3, function(sp) {
  sp_clim <- demography_inf_ppt$clim[demography_inf_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out=30)
})
names(climate_range_per_species) <- 1:3

delta_inf_species_range <- lapply(1:3, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_inf(cl, posterior_samples_inf) %>% filter(species==sp)
  }) %>% bind_rows()
}) %>% bind_rows()

delta_inf_summary <- delta_inf_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(median_delta=median(delta), lower_90=quantile(delta,0.05),
            upper_90=quantile(delta,0.95), prob_delta_gt0=mean(delta>0), .groups="drop") %>%
  mutate(
    species = factor(species, levels=1:3,
                     labels=c("Agrostis hyemalis","Elymus virginicus","Poa autumnalis")),
    herb    = factor(herb, levels=c(0,1),
                     labels=c("Herbivory access","Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  )

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()
for(sp in unique(delta_inf_summary$species)) {
  for(h in unique(delta_inf_summary$herb)) {
    df_sub <- delta_inf_summary %>%
      filter(species == sp, herb == h)
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}
delta_inf_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm        = round(clim_mm, 0),
    median_delta   = round(median_delta, 3),
    lower_90       = round(lower_90, 3),
    upper_90       = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

delta_long_inf <- delta_inf_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(cols=c(median_delta, prob_delta_gt0), names_to="metric", values_to="value") %>%
  mutate(
    metric = recode(metric,
                    "median_delta"   = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "*Agrostis hyemalis*",
      species == "Elymus virginicus" ~ "*Elymus virginicus*",
      species == "Poa autumnalis"    ~ "*Poa autumnalis*"
    )
  )

# ── LaTeX table: inflorescence ────────────────────────────────────────────────
delta_inf_latex <- delta_inf_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm  = clim_mm,
    Median_delta      = median_delta,
    Lower_90          = lower_90,
    Upper_90          = upper_90,
    Prob_delta_gt0    = prob_delta_gt0
  )

delta_inf_xt <- xtable(
  delta_inf_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)
colnames(delta_inf_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)
print(
  delta_inf_xt,
  include.rownames       = FALSE,
  sanitize.text.function = identity,
  floating               = FALSE
)

# ══════════════════════════════════════════════════════════════════════════════
# SPIKELET
# ══════════════════════════════════════════════════════════════════════════════
demography_climate %>%
  filter(Species %in% c("ELVI","POAU"), tiller_t1 > 0, inf_t1 > 0) %>%
  dplyr::select(Species, Population, Site, site_year, Plot, site_species_plot,
                Endo, Herbivory, tiller_t, spikelet_t1, cum_ppt, ppt_scaled) %>%
  na.omit() %>%
  mutate(
    Site              = as.integer(factor(Site)),
    site_year         = as.integer(factor(site_year)),
    Species           = as.integer(factor(Species)),
    Population        = as.integer(factor(Population)),
    site_species_plot = as.integer(factor(site_species_plot)),
    Endo              = as.integer(Endo),
    Herbivory         = as.integer(Herbivory),
    log_size_t0       = log(tiller_t),
    ppt               = ppt_scaled
  ) -> demography_climate_spik

demography_spik_ppt <- list(
  nSpp       = n_distinct(demography_climate_spik$Species),
  nSite      = n_distinct(demography_climate_spik$Site),
  nsite_year = n_distinct(demography_climate_spik$site_year),
  nPop       = n_distinct(demography_climate_spik$Population),
  nPlot      = n_distinct(demography_climate_spik$site_species_plot),
  Spp        = demography_climate_spik$Species,
  site       = demography_climate_spik$Site,
  site_year  = demography_climate_spik$site_year,
  pop        = demography_climate_spik$Population,
  plot       = demography_climate_spik$site_species_plot,
  clim       = as.vector(demography_climate_spik$ppt),
  endo       = demography_climate_spik$Endo,
  herb       = demography_climate_spik$Herbivory,
  size       = demography_climate_spik$log_size_t0,
  y          = demography_climate_spik$spikelet_t1,
  N          = nrow(demography_climate_spik)
)

fit_spik_ppt <- readRDS(url("https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1"))
posterior_samples_spik <- rstan::extract(fit_spik_ppt)

climate_range_per_species_spik <- lapply(1:2, function(sp) {
  sp_clim <- demography_spik_ppt$clim[demography_spik_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out=30)
})

predictions <- do.call(rbind, lapply(1:2, function(sp) {
  expand.grid(clim=climate_range_per_species_spik[[sp]],
              endo=c(0,1), herb=c(0,1), species=sp)
}))

get_predictions_spik <- function(clim, endo, herb, species_index, ps) {
  with(ps, {
    eta <- b0[,species_index]+bendo[,species_index]*endo+bherb[,species_index]*herb+
      bclim[,species_index]*clim+bendoclim[,species_index]*clim*endo+
      bendoherb[,species_index]*endo*herb+bherbclim[,species_index]*herb*clim+
      bendoherbclim[,species_index]*endo*herb*clim
    exp(eta)
  })
}

n_post_spik <- nrow(posterior_samples_spik$b0)
pred_matrix_spik <- matrix(NA, nrow=nrow(predictions), ncol=n_post_spik)
for (i in seq_len(nrow(predictions))) {
  pred_matrix_spik[i,] <- get_predictions_spik(
    predictions$clim[i], predictions$endo[i],
    predictions$herb[i], predictions$species[i], posterior_samples_spik
  )
}

pred_spik_long <- cbind(predictions, as.data.frame(pred_matrix_spik)) %>%
  pivot_longer(cols=starts_with("V"), names_to="Posterior_Sample", values_to="Prediction")

plot_data_spik <- pred_spik_long %>%
  group_by(species, clim, endo, herb) %>%
  summarise(mean=mean(Prediction), lower_90=quantile(Prediction,0.05),
            upper_90=quantile(Prediction,0.95), .groups="drop")

delta_spik <- plot_data_spik %>%
  pivot_wider(names_from=endo, values_from=c(mean, lower_90, upper_90)) %>%
  mutate(mean=mean_1-mean_0, lower_90=lower_90_1-upper_90_0,
         upper_90=upper_90_1-lower_90_0, panel="Δ (S+ - S-)") %>%
  dplyr::select(species, clim, herb, mean, lower_90, upper_90, panel)

plot_data_spik <- plot_data_spik %>% mutate(panel="Spikelets") %>% bind_rows(delta_spik)

species_levels_2 <- c("italic('Elymus virginicus')", "italic('Poa autumnalis')")
plot_data_spik <- plot_data_spik %>%
  mutate(species = factor(species, labels=species_levels_2))

# ── Observed spikelets: plot means only ───────────────────────────────────────
observed_data_spik <- demography_spik_ppt %>%
  data.frame(clim=.$clim, endo=.$endo, herb=.$herb,
             species=.$Spp, plot=.$plot, y=.$y) %>%
  group_by(plot, species, herb, clim, endo) %>%
  summarise(
    y_plot_mean = mean(y, na.rm=TRUE),
    n_obs       = sum(!is.na(y)),
    .groups     = "drop"
  ) %>%
  mutate(
    panel      = "Spikelets",
    climate_mm = exp(clim * ppt_sd + ppt_mean)
  )

observed_data_spik$species <- factor(observed_data_spik$species,
                                     levels=1:2, labels=species_levels_2)
plot_data_spik <- plot_data_spik %>% mutate(climate_mm = exp(clim * ppt_sd + ppt_mean))

# Pre-compute jitter offset
set.seed(13)
observed_data_spik <- observed_data_spik %>%
  mutate(jitter_x = climate_mm + runif(n(), -20, 20))

y_limits_spik <- plot_data_spik %>%
  filter(panel == "Δ (S+ - S-)") %>%
  group_by(species) %>%
  summarise(
    ymin = min(0, min(lower_90, na.rm = TRUE)),
    ymax = max(0, max(upper_90, na.rm = TRUE)),
    .groups = "drop"
  )

panel_labels_spik <- data.frame(
  species = rep(species_levels_2, each=2),
  herb    = rep(c(0,1), times=2),
  label   = c("(a)","(b)","(c)","(d)"),
  panel   = "Spikelets"
)

# ── Zero-line segments ────────────────────────────────────────────────────────
make_dashes <- function(xmin, xmax, n = 8, dash_prop = 0.45) {
  step <- (xmax - xmin) / n
  tibble(
    x = xmin + (0:(n - 1)) * step,
    xend = xmin + (0:(n - 1)) * step + step * dash_prop,
    y = 0,
    yend = 0
  )
}

global_xmin_spik <- min(plot_data_spik$climate_mm, na.rm = TRUE)
global_xmax_spik <- max(plot_data_spik$climate_mm, na.rm = TRUE)

zero_dashes_spik <- expand.grid(
  species = levels(plot_data_spik$species),
  herb = c(0, 1)
) %>%
  mutate(
    dashes = purrr::map(
      seq_len(n()),
      ~ make_dashes(global_xmin_spik, global_xmax_spik)
    )
  ) %>%
  tidyr::unnest(dashes)


Cairo::CairoTIFF(
  file.path(FIG_DIR, "Spikelet_diff.tiff"),
  width = 7, height = 5.5, units = "in", dpi = 600
)
ggplot(plot_data_spik) +
  geom_line(
    data=subset(plot_data_spik, panel=="Spikelets"),
    aes(x=climate_mm, y=mean, color=factor(endo), group=endo), linewidth=0.5
  ) +
  geom_ribbon(
    data=subset(plot_data_spik, panel=="Spikelets"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90, fill=factor(endo), group=endo),
    alpha=0.3, color=NA
  ) +
  # ── Observed: plot means only, size ∝ n_obs ──────────────────────────────
  geom_point(
    data=subset(observed_data_spik, panel=="Spikelets"),
    aes(x=jitter_x, y=y_plot_mean, color=factor(endo), size=n_obs),
    alpha=0.5, show.legend=FALSE
  ) +
  scale_size_continuous(range = c(0.8, 6.7)) +  # was c(0.5, 4)
  geom_line(
    data=subset(plot_data_spik, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, y=mean), color="black", linewidth=0.5
  ) +
  geom_ribbon(
    data=subset(plot_data_spik, panel=="Δ (S+ - S-)"),
    aes(x=climate_mm, ymin=lower_90, ymax=upper_90),
    fill="#D9BFD6", alpha=0.5
  ) +
  geom_segment(
    data = zero_dashes_spik,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.55,
    color = "black",
    lineend = "butt",
    inherit.aes = FALSE
  ) +
  ggh4x::facet_nested(
    panel ~ species + herb, scales="free_y", space="fixed",
    labeller=labeller(species=label_parsed,
                      herb=c("0"="Herbivory access","1"="Herbivory exclusion"),
                      panel=ggplot2::as_labeller(PANEL_LABELS, default = label_parsed))
  ) +
  ggh4x::facetted_pos_scales(y=list(
    panel=="Δ (S+ - S-)" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(minor_breaks=NULL,
                         limits=c(y_limits_spik$ymin[y_limits_spik$species=="italic('Elymus virginicus')"],
                                  y_limits_spik$ymax[y_limits_spik$species=="italic('Elymus virginicus')"]),
                         breaks = scales::pretty_breaks(n = 4),
                         expand=c(0,0)),
    panel=="Δ (S+ - S-)" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(minor_breaks=NULL,
                         limits=c(y_limits_spik$ymin[y_limits_spik$species=="italic('Poa autumnalis')"],
                                  y_limits_spik$ymax[y_limits_spik$species=="italic('Poa autumnalis')"]),
                         breaks = scales::pretty_breaks(n = 4),
                         expand=c(0,0)),
    panel=="Spikelets" & species=="italic('Elymus virginicus')" ~
      scale_y_continuous(limits=c(0,50), expand=c(0,0)),
    panel=="Spikelets" & species=="italic('Poa autumnalis')" ~
      scale_y_continuous(limits=c(0,60), expand=c(0,0))
  )) +
  labs(x="Precipitation (mm)",
       y = expression(paste("Number of spikelets per inflorescence / ", Delta, " spikelets (",
                            italic(S)^{"+"} - italic(S)^{"\u2212"}, ")")),
       color="Symbiont", fill="Symbiont") +
  scale_color_manual(values=ENDO_COLORS, labels=ENDO_LABELS) +
  scale_fill_manual(values=ENDO_COLORS,  labels=ENDO_LABELS) +
  vr_theme() +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.27),
    axis.line         = element_line(color = "black", linewidth = 0.13),
    axis.title        = element_markdown(size = 11),
    axis.title.y      = element_text(size = 11),
    axis.text         = element_text(size = 6),
    axis.ticks.x      = element_line(color = "black", linewidth = 0.27),
    axis.ticks.y      = element_line(color = "black", linewidth = 0.27),
    legend.title      = element_text(size = 10),
    legend.text       = element_markdown(size = 10),
    panel.spacing.y   = unit(0.2, "cm"),
    text              = element_text(family = "Arial"),
    strip.text.x      = element_text(size = 12, color = "black"),
    strip.text.y      = element_text(size = 11, color = "black"),
    strip.background  = element_rect(color = "black", fill = "grey90", linewidth = 0.27),
    legend.position   = c(0.12, 0.88),
    #legend.direction  = "horizontal",
    legend.justification = "center"
  )+
  geom_text(data = panel_labels_spik,
            aes(x = 490, y = 47, label = label),
            fontface = "plain", size = 4, hjust = 0, inherit.aes = FALSE)
dev.off()

# ── Delta spikelet summary ────────────────────────────────────────────────────
compute_delta_spik <- function(clim, ps, herb_values=c(0,1)) {
  n_species <- dim(ps$b0)[2]; n_post <- dim(ps$b0)[1]
  lapply(1:n_species, function(sp) {
    lapply(herb_values, function(h) {
      eta_p <- ps$b0[,sp]+ps$bendo[,sp]*1+ps$bherb[,sp]*h+ps$bclim[,sp]*clim+
        ps$bendoclim[,sp]*clim*1+ps$bendoherb[,sp]*1*h+
        ps$bherbclim[,sp]*h*clim+ps$bendoherbclim[,sp]*1*h*clim
      eta_m <- ps$b0[,sp]+ps$bendo[,sp]*0+ps$bherb[,sp]*h+ps$bclim[,sp]*clim+
        ps$bendoclim[,sp]*clim*0+ps$bendoherb[,sp]*0*h+
        ps$bherbclim[,sp]*h*clim+ps$bendoherbclim[,sp]*0*h*clim
      data.frame(species=sp, herb=h, clim=clim,
                 Posterior_Sample=1:n_post, delta=exp(eta_p)-exp(eta_m))
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
}

climate_range_per_species <- lapply(1:2, function(sp) {
  sp_clim <- demography_spik_ppt$clim[demography_spik_ppt$Spp == sp]
  seq(min(sp_clim), max(sp_clim), length.out=30)
})
names(climate_range_per_species) <- 1:2

delta_spik_species_range <- lapply(1:2, function(sp) {
  lapply(climate_range_per_species[[sp]], function(cl) {
    compute_delta_spik(cl, posterior_samples_spik) %>% filter(species==sp)
  }) %>% bind_rows()
}) %>% bind_rows()

delta_spik_summary <- delta_spik_species_range %>%
  group_by(species, herb, clim) %>%
  summarise(median_delta=median(delta), lower_90=quantile(delta,0.05),
            upper_90=quantile(delta,0.95), prob_delta_gt0=mean(delta>0), .groups="drop") %>%
  mutate(
    species = factor(species, levels=1:3,
                     labels=c("Agrostis hyemalis","Elymus virginicus","Poa autumnalis")),
    herb    = factor(herb, levels=c(0,1),
                     labels=c("Herbivory access","Herbivory exclusion")),
    clim_mm = exp(clim * ppt_sd + ppt_mean)
  )

# Prepare table with 5 representative climate points per species × herb
filtered_rows <- list()
for(sp in unique(delta_spik_summary$species)) {
  for(h in unique(delta_spik_summary$herb)) {
    df_sub <- delta_spik_summary %>%
      filter(species == sp, herb == h)
    clim_pts <- quantile(df_sub$clim_mm, probs = c(0, 0.25, 0.5, 0.75, 1))
    for(pt in clim_pts) {
      closest_row <- df_sub[which.min(abs(df_sub$clim_mm - pt)), ]
      filtered_rows <- append(filtered_rows, list(closest_row))
    }
  }
}
delta_spik_filtered <- bind_rows(filtered_rows) %>%
  mutate(
    clim_mm        = round(clim_mm, 0),
    median_delta   = round(median_delta, 3),
    lower_90       = round(lower_90, 3),
    upper_90       = round(upper_90, 3),
    prob_delta_gt0 = round(prob_delta_gt0, 3)
  )

delta_long_spik <- delta_spik_summary %>%
  dplyr::select(species, herb, clim_mm, median_delta, prob_delta_gt0) %>%
  pivot_longer(cols=c(median_delta, prob_delta_gt0), names_to="metric", values_to="value") %>%
  mutate(
    metric = recode(metric,
                    "median_delta"   = "Median Δ (S+ − S−)",
                    "prob_delta_gt0" = "Pr (Δ > 0)"),
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "*Agrostis hyemalis*",
      species == "Elymus virginicus" ~ "*Elymus virginicus*",
      species == "Poa autumnalis"    ~ "*Poa autumnalis*"
    )
  )

# ── LaTeX table: spikelet ─────────────────────────────────────────────────────
delta_spik_latex <- delta_spik_filtered %>%
  mutate(
    Species = species_abbrev[as.character(species)],
    Species = paste0("\\textit{", Species, "}"),
    Herbivore_treatment = case_when(
      herb == "Herbivory access"    ~ "Access",
      herb == "Herbivory exclusion" ~ "Exclusion",
      TRUE ~ as.character(herb)
    )
  ) %>%
  arrange(Species, Herbivore_treatment, clim_mm) %>%
  dplyr::select(
    Species,
    Herbivore_treatment,
    Precipitation_mm  = clim_mm,
    Median_delta      = median_delta,
    Lower_90          = lower_90,
    Upper_90          = upper_90,
    Prob_delta_gt0    = prob_delta_gt0
  )

delta_spik_xt <- xtable(
  delta_spik_latex,
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)
colnames(delta_spik_xt) <- c(
  "Species",
  "\\makecell{Herbivore \\\\ treatment}",
  "Precipitation",
  "Median $\\Delta$",
  "Lower 90\\%",
  "Upper 90\\%",
  "P($\\Delta>0$)"
)
print(
  delta_spik_xt,
  include.rownames       = FALSE,
  sanitize.text.function = identity,
  floating               = FALSE
)

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED Pr(Δ > 0) PANEL FIGURE
# ══════════════════════════════════════════════════════════════════════════════
delta_long_surv$trait <- "Survival"
delta_long_grow$trait <- "Growth"
delta_long_inf$trait  <- "Inflorescence"

delta_long_all <- dplyr::bind_rows(delta_long_surv, delta_long_grow, delta_long_inf) %>%
  mutate(
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('A. hyemalis')",
      species == "Elymus virginicus" ~ "italic('E. virginicus')",
      species == "Poa autumnalis"    ~ "italic('P. autumnalis')",
      TRUE ~ as.character(species)
    )
  )

panel_labels <- delta_long_all %>%
  filter(metric == "Pr (Δ > 0)") %>%
  distinct(trait, species_label) %>%
  arrange(trait, species_label) %>%
  mutate(label = paste0("(", letters[1:n()], ")"))

make_dashes <- function(xmin, xmax, y = 0.5, n = 14, dash_prop = 0.70) {
  step <- (xmax - xmin) / n
  tibble(
    x    = xmin + (0:(n - 1)) * step,
    xend = xmin + (0:(n - 1)) * step + step * dash_prop,
    y    = y,
    yend = y
  )
}

global_xmin <- min(delta_long_all$clim_mm, na.rm = TRUE)
global_xmax <- max(delta_long_all$clim_mm, na.rm = TRUE)

reference_dashes <- make_dashes(
  global_xmin,
  global_xmax,
  y = 0.5
)


p_lower <- delta_long_all %>%
  filter(metric == "Pr (Δ > 0)") %>%
  ggplot(aes(x=clim_mm, y=value, color=herb, group=herb)) +
  geom_line(linewidth=0.6) +
  geom_segment(
    data = reference_dashes,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35,
    color = "black",
    lineend = "butt",
    inherit.aes = FALSE
  ) +
  facet_grid(
    trait ~ species_label, scales="free_x",
    labeller=labeller(species_label=label_parsed, trait=label_value)
  ) +
  scale_color_manual(values=c("Herbivory access"="#E69F00",
                               "Herbivory exclusion"="#009E73")) +
  labs(x="Precipitation (mm)", y="P(Δ > 0)", color="Herbivore treatment") +
  theme_classic(base_size=10) +
  geom_text(
    data = panel_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.2,
    vjust = 1.2,
    size = 4,
    fontface = "plain"
  )+
  theme(
    panel.border     = element_rect(color="black", fill=NA, linewidth=0.2),
    axis.line        = element_line(color="black", linewidth=0.1),
    legend.position  = c(0.18, 0.76),
    legend.title     = element_text(size=8),
    legend.text      = element_text(size=8),
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.3, "cm"),
    panel.spacing.y  = unit(0.2, "cm"),
    axis.title       = element_text(size=10),
    axis.text        = element_text(size=6),
    axis.ticks.x     = element_line(color="black", linewidth=0.2),
    axis.ticks.y     = element_line(color="black", linewidth=0.2),
    text             = element_text(family="Arial"),
    strip.text.x     = element_text(size=12, color="black"),
    strip.text.y     = element_text(size=10, color="black"),
    strip.background = element_rect(color="black", fill="grey90", linewidth=0.2)
  )

Cairo::CairoTIFF(
  file.path(FIG_DIR, "All_traits_diff_stat_lower.tiff"),
  width = 7, height = 6, units = "in", dpi = 600
)
print(p_lower)
dev.off()

# Make sure only the two species are retained
delta_long_spik <- delta_long_spik %>%
  filter(species %in% c("Agrostis hyemalis", "Elymus virginicus")) %>%
  mutate(
    species_label = case_when(
      species == "Agrostis hyemalis" ~ "italic('Agrostis hyemalis')",
      species == "Elymus virginicus" ~ "italic('Elymus virginicus')"
    )
  )

# Create dashed reference line
make_dashes <- function(xmin, xmax, y, n = 14, dash_prop = 0.70) {
  step <- (xmax - xmin) / n
  
  tibble(
    x    = xmin + (0:(n - 1)) * step,
    xend = xmin + (0:(n - 1)) * step + step * dash_prop,
    y    = y,
    yend = y
  )
}

global_xmin <- min(delta_long_spik$clim_mm, na.rm = TRUE)
global_xmax <- max(delta_long_spik$clim_mm, na.rm = TRUE)

# Dashes at 0 for Median Δ panel
reference_dashes_zero <- make_dashes(
  global_xmin,
  global_xmax,
  y = 0
) %>%
  mutate(metric = "Median Δ (S+ − S−)")

# Dashes at 0.5 for Pr(Δ > 0) panel
reference_dashes_half <- make_dashes(
  global_xmin,
  global_xmax,
  y = 0.5
) %>%
  mutate(metric = "Pr (Δ > 0)")


Cairo::CairoTIFF(
  file.path(FIG_DIR, "Spike_diff_stat.tiff"),
  width = 6, height = 5, units = "in", dpi = 600
)

ggplot(
  delta_long_spik %>%
    filter(!is.na(species_label)),
  aes(x = clim_mm, y = value, color = herb, group = herb)
) +
  
  geom_line(linewidth = 0.5) +
  
  # Reference line: Median Δ = 0
  geom_segment(
    data = reference_dashes_zero,
    aes(
      x = x, xend = xend,
      y = y, yend = y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  
  # Reference line: Pr(Δ > 0) = 0.5
  geom_segment(
    data = reference_dashes_half,
    aes(
      x = x, xend = xend,
      y = y, yend = y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  
  facet_grid(
    metric ~ species_label,
    scales = "free_y",
    labeller = labeller(
      species_label = label_parsed,
      metric = label_value
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Herbivory access" = "#E69F00",
      "Herbivory exclusion" = "#009E73"
    )
  ) +
  
  labs(
    x = "Precipitation (mm)",
    y = NULL,
    color = "Herbivore exclusion"
  ) +
  
  vr_theme() +
  
  theme(
    legend.position = "bottom",
    panel.spacing.y = unit(0.02, "cm")
  )

dev.off()