# Plot jet radius vs. axial position z for the last frame of the simulation.
# Reads microjet_last_frame.csv (written by main.cpp) and uses the
# create_line_plot template from plotTemplates.R.
# Requires: pip install rpy2

import rpy2.robjects as ro

# load the templates and the data (columns: node, z, h, v)
ro.r('source("plotTemplates.R")')
ro.r('lastFrame <- read.csv("microjet_last_frame.csv")')

# convert to micrometers so the axes are readable
ro.r('lastFrame$z_um <- lastFrame$z * 1e6')
ro.r('lastFrame$h_um <- lastFrame$h * 1e6')

# make the plot with the template and save it
ro.r('''
p <- create_line_plot(lastFrame, "z_um", "h_um",
                      title = "Jet Radius at Final Time",
                      x_label = "z (µm)",
                      y_label = "Radius h (µm)")
ggsave("microjet_last_frame.png", p, width = 10, height = 6, dpi = 300)
''')

print("Saved microjet_last_frame.png")
