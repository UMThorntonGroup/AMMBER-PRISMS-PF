#!/usr/bin/env python3

# run with:
# visit -cli -nowin -s plot_and_save.py
import sys
import os
from visit import *
from visit_utils import encoding

#DeleteAllPlots()
dir_name = "frames"
phases = {"FCC_A1_0_0": "Reds", "AL2CU_1_0": "Blues"} # "LIQUID_0": "Purples"
concentration = "c_CU"
conc_color = "gray" # "turbo"
conc_minmax = (0, 0.333333) # min, max
hide_legend = True

# Step 1: Open a database (the whole .vtu time series)
dbname="solution-*.vtu database"
OpenDatabase(dbname)

#Step 2B: Pseudocolor plots
for phase, color in phases.items():
    AddPlot("Pseudocolor", phase, 1, 1)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 1
    PseudocolorAtts.colorTableName = color
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.Ramp
    SetPlotOptions(PseudocolorAtts)

AddPlot("Pseudocolor", concentration, 1, 1)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = conc_minmax[0]
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = conc_minmax[1]
PseudocolorAtts.colorTableName = conc_color
PseudocolorAtts.invertColorTable = 0
SetPlotOptions(PseudocolorAtts)

#Step 2C: Adjust annotation and view
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.axes2D.xAxis.title.visible = 0
AnnotationAtts.axes2D.yAxis.title.visible = 0
AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
AnnotationAtts.legendInfoFlag = not hide_legend
SetAnnotationAttributes(AnnotationAtts)

View2DAtts = View2DAttributes()
View2DAtts.viewportCoords = (0.05, 0.95, 0.05, 0.90) # (xmin, xmax, ymin, ymax)
View2DAtts.fullFrameActivationMode = View2DAtts.On  # On, Off, Auto
SetView2D(View2DAtts)
ResetView()

#ViewCurveAtts = ViewCurveAttributes()
#View2DAtts.viewportCoords = (0.2, 0.95, 0.10, 0.95)
#SetView2D(View2DAtts)

# Step 3: Draw the plots
DrawPlots()

sys.exit()
