# Plot 2D symmetric jet profile vs. axial position z for the frame prior to breakup.

import rpy2.robjects as ro

# load the data (columns: node, z, h, v)
ro.r('source("plotTemplates.R")')
ro.r('lastFrame <- read.csv("microjet_last_frame.csv")')

# convert to micrometers
ro.r('lastFrame$z_um <- lastFrame$z * 1e6')
ro.r('lastFrame$h_um <- lastFrame$h * 1e6')

ro.r('''
# 2D
p_profile <- ggplot(lastFrame, aes(x = z_um)) +
    geom_ribbon(aes(ymin = -h_um, ymax = h_um), fill = "#3366cc", color = "#1a4080", linewidth = 0.5) +
    geom_line(aes(y = h_um), color = "#1a4080", linewidth = 0.6) +
    geom_line(aes(y = -h_um), color = "#1a4080", linewidth = 0.6) +
    coord_fixed(ratio = 1) +
    labs(
        title = "Microjet Profile at First Breakup",
        subtitle = "Van Hoeve et al. (18.5 µm glycerol-saline jet, v = 5.4 m/s)",
        x = "Axial Distance z (µm)",
        y = "Radius r (µm)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        text = element_text(family = plot_config$font_family)
    )

ggsave("microjet_last_frame.png", p_profile, width = 12, height = 3.5, dpi = 300)

# 1D
p_radius <- create_line_plot(lastFrame, "z_um", "h_um",
    title = "Jet Radius Profile at First Breakup",
    x_label = "z (µm)",
    y_label = "Radius h (µm)") +
    expand_limits(y = 0)
ggsave("microjet_last_frame_radius.png", p_radius, width = 10, height = 5, dpi = 300)
''')

print("Saved.")