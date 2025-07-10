#!/usr/bin/env python3

# run with:
# visit -cli -nowin -s plot_and_save.py
import sys
import os
from visit import *
from visit_utils import encoding

#DeleteAllPlots()
dir_name = "frames"
phases = {"Fe2O3_0": "Reds", "MoO3_0": "Blues", "Fe2MoO43_0": "Oranges", "liquid_0": "Purples"}

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
    PseudocolorAtts.max = 2
    PseudocolorAtts.colorTableName = color
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.Ramp
    SetPlotOptions(PseudocolorAtts)

#Step 2C: Adjust annotation and view
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.axes2D.xAxis.title.visible = 0
AnnotationAtts.axes2D.yAxis.title.visible = 0
AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
AnnotationAtts.legendInfoFlag = 0
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

# Step 4: Animate through time and save results

# Clear any existing frames in the directory
if os.path.exists(dir_name):
    for f in os.listdir(dir_name):
        file_path = os.path.join(dir_name, f)
        if os.path.isfile(file_path):
            os.remove(file_path)
# Make the directory 
os.makedirs(dir_name, exist_ok=True)
for states in range(TimeSliderGetNStates()):
    #Set slider to state
    SetTimeSliderState(states)
    # Get the time corresponding to the state
    Query("Time")
    # Assign this time to the variable "t"
    t = GetQueryOutputValue()
    # Print the state number, time and phase fraction to
    # screen and to files
    print("Saving frame % d, time %.1f" %(states, t))
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.fileName = f'{dir_name}/frame_'
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
DeleteAllPlots()
CloseDatabase(dbname)

# Make movie
#for ext in ["mpg","wmv","mov"]:
#   encodinprintg.encode(f'{dir_name}/frame_*',"movie." + ext)

sys.exit()
