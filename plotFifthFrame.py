# Plot jet radius vs. axial position z for the fifth frame of the simulation.
# Reads microjet_fifth_frame.csv (written by main.cpp) and uses the
# create_line_plot template from plotTemplates.R.

import rpy2.robjects as ro

# load the templates and the data (columns: node, z, h, v)
ro.r('source("plotTemplates.R")')
ro.r('fifthFrame <- read.csv("microjet_5th_frame.csv")')

# convert to micrometers so the axes are readable
ro.r('fifthFrame$z_um <- fifthFrame$z * 1e6')
ro.r('fifthFrame$h_um <- fifthFrame$h * 1e6')

ro.r('''
p <- create_line_plot(fifthFrame, "z_um", "h_um",
    title = "Jet Radius at Final Time",
    x_label = "z (µm)",
    y_label = "Radius h (µm)") +
    expand_limits(y = 0)
    ggsave("microjet_fifth_frame.png", p, width = 10, height = 6, dpi = 300)
''')

print("Saved microjet_fifth_frame.png")
