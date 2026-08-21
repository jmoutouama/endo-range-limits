library(ggplot2)
library(patchwork)

# ── Palette (colorblind-safe, Ecology Letters compatible) ──────────────────────
COL_HERB <- "#E69F00"
COL_EXCL <- "#009E73"

# ── Shared theme ───────────────────────────────────────────────────────────────
theme_el <- function(show_legend = FALSE) {
  theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(
      # Panel
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(colour = "#3C3C3C", fill = NA,
                                       linewidth = 0.55),
      panel.grid        = element_blank(),
      axis.line         = element_blank(),

      # Title – italic, not bold, lighter weight
      plot.title = element_text(
        face   = "italic",
        size   = 10.5,
        hjust  = 0.5,
        margin = margin(b = 5),
        colour = "#1A1A1A"
      ),

      # Axes
      axis.title       = element_text(size = 12, colour = "#1A1A1A"),
      axis.text        = element_text(size = 9,    colour = "#3C3C3C"),
      axis.ticks       = element_line(colour = "#3C3C3C", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),

      # Legend – bottom of full figure (controlled at plot_annotation level)
      legend.position   = if (show_legend) c(0.70, 0.25) else "none",
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.width  = unit(1.4, "cm"),
      legend.text       = element_text(size = 9.5, family = "Helvetica"),
      legend.title      = element_blank(),

      # Margins
      plot.margin = margin(6, 8, 4, 8)
    )
}

# ── Panel function ─────────────────────────────────────────────────────────────
make_panel <- function(title,
                       herb_i, herb_s,
                       excl_i, excl_s,
                       show_x      = FALSE,
                       show_y      = FALSE,
                       show_legend = FALSE) {

  x <- seq(0, 1, length.out = 300)

  df <- data.frame(
    x = rep(x, 2),
    y = c(herb_i + herb_s * x,
          excl_i + excl_s * x),
    treatment = rep(c("Herbivory access", "Herbivory exclusion"), each = 300)
  )

  df$treatment <- factor(df$treatment,
                         levels = c("Herbivory access", "Herbivory exclusion"))

  ggplot(df, aes(x = x, y = y,
                 colour   = treatment,
                 linetype = treatment)) +

    # Zero reference line – more visible than before
    geom_hline(yintercept = 0,
               colour    = "#888888",
               linewidth = 0.55,
               linetype  = "dashed") +

    # Slightly thicker lines for print legibility
    geom_line(linewidth = 1.3) +

    scale_colour_manual(
      values = c("Herbivory access"    = COL_HERB,
                 "Herbivory exclusion" = COL_EXCL)
    ) +

    scale_linetype_manual(
      values = c("Herbivory access"    = "solid",
                 "Herbivory exclusion" = "solid")
    ) +

    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Dry", "Wet"),
      expand = c(0.05, 0.05)
    ) +

    scale_y_continuous(
      breaks = 0,
      labels = "0",
      limits = c(-0.65, 0.82),
      expand = c(0.03, 0.03)
    ) +

    labs(
      title = title,
      x     = if (show_x) "Precipitation" else NULL,
      y     = if (show_y) "Symbiont effect (S+ \u2212 S-)" else NULL
    ) +

    theme_el(show_legend = show_legend)
}

# ── Panels ─────────────────────────────────────────────────────────────────────
p1 <- make_panel("Context-independent mutualism",
                 0.58,  0.00,
                 0.50,  0.00,
                 show_y = TRUE)

p2 <- make_panel("Herbivory-dependent",
                 0.58,  0.00,
                 0.10,  0.00,
                 show_legend = TRUE)   # legend pulled to figure bottom below

p3 <- make_panel("Climate-dependent",
                 0.55, -0.80,
                 0.45, -0.80,
                 show_x = TRUE, show_y = TRUE)

p4 <- make_panel("Herbivory \u00d7 Climate-dependent",
                 0.60, -1.00,
                -0.03, -0.10,
                 show_x = TRUE)

# ── Combine ────────────────────────────────────────────────────────────────────
# Extract legend from p2 and place it at the bottom of the combined figure
final <- (p1 | p2) / (p3 | p4) +

  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        face   = "bold",
        size   = 11,
        family = "Helvetica",
        colour = "#1A1A1A"
      ),
      plot.tag.position = c(0.01, 0.99)
    )
  ) &

  # Apply consistent spacing across all panels
  theme(plot.margin = margin(6, 10, 4, 10))

# ── Save ───────────────────────────────────────────────────────────────────────
# Update path to your local directory
out_path <- "/Users/jacobmoutouama/Dropbox/Miller Lab/github/endo-range-limits/Figure/conceptual_figure.pdf"

ggsave(out_path,
       plot   = final,
       width  = 8,
       height = 6.0,        # slightly tighter than before
       device = cairo_pdf)

# Also save a high-res PNG for quick preview
ggsave(sub("\\.pdf$", ".png", out_path),
       plot   = final,
       width  = 7.5,
       height = 6.0,
       dpi    = 300)

message("Saved: conceptual_figure.pdf and conceptual_figure.png")
