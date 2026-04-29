library(ggplot2)
library(patchwork)
library(cowplot)

# ── Palette ────────────────────────────────────────────────────────────────────
COL_HERB <- "#E69F00"
COL_EXCL <- "#009E73"

# ── Panel function (NO manual labels anymore) ───────────────────────────────────
make_panel <- function(title,
                       herb_i, herb_s,
                       excl_i, excl_s,
                       show_x = FALSE,
                       show_y = FALSE,
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
                 colour = treatment,
                 linetype = treatment)) +
    
    geom_hline(yintercept = 0,
               colour = "#BDBDBD",
               linewidth = 0.6,
               linetype = "dashed") +
    
    geom_line(linewidth = 1.0, alpha = 0.95) +
    
    scale_colour_manual(values = c("Herbivory access" = COL_HERB,
                                   "Herbivory exclusion" = COL_EXCL)) +
    
    scale_linetype_manual(values = c("Herbivory access" = "solid",
                                     "Herbivory exclusion" = "dashed")) +
    
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
      x = if (show_x) "Precipitation" else NULL,
      y = if (show_y) "Symbiont effect (\u0394)" else NULL
    ) +
    
    theme_classic(base_size = 11, base_family = "serif") +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5,
        margin = margin(b = 6),
        colour = "#1A1A1A"
      ),
      
      axis.title = element_text(size = 9.5),
      axis.text  = element_text(size = 9, colour = "#2c2c2c"),
      
      axis.line = element_blank(),
      
      panel.background = element_rect(fill = "#F7F7F7", colour = NA),
      panel.border     = element_rect(colour = "#D0D0D0", fill = NA, linewidth = 0.6),
      
      axis.ticks = element_line(colour = "#2c2c2c", linewidth = 0.4),
      axis.ticks.length = unit(3.5, "pt"),
      
      legend.position = if (show_legend) c(0.72, 0.25) else "none",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 9.5, family = "serif"),
      legend.title = element_blank()
    )
}

# ── Panels ────────────────────────────────────────────────────────────────────
p1 <- make_panel("Context-independent\nmutualism",
                 0.58, 0.00,
                 0.50, 0.00,
                 show_y = TRUE)

p2 <- make_panel("Herbivory-dependent",
                 0.58, 0.00,
                 0.10, 0.00,
                 show_legend = TRUE)

p3 <- make_panel("Climate-dependent",
                 0.55, -0.80,
                 0.45, -0.80,
                 show_x = TRUE, show_y = TRUE)

p4 <- make_panel("Herbivory × Climate-\ndependent",
                 0.60, -1.00,
                 -0.03, -0.10,
                 show_x = TRUE)

# ── COMBINE with OUTSIDE PANEL LABELS ─────────────────────────────────────────
final <- (p1 | p2) / (p3 | p4) +
  
  plot_annotation(
    tag_levels = "a",   # THIS gives (a), (b), (c), (d)
    tag_prefix = "(",
    tag_suffix = ")",
    
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 12,
        family = "serif"
      ),
      plot.tag.position = c(0, 1)  # top-left OUTSIDE panel
    )
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("conceptual_figure.pdf",
       plot = final,
       width = 7.5, height = 6.5,
       device = cairo_pdf)

ggsave("conceptual_figure.png",
       plot = final,
       width = 7.5, height = 6.5,
       dpi = 300,
       bg = "white")

message("Saved: conceptual_figure.pdf and conceptual_figure.png")