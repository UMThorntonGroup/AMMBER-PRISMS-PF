#!/usr/bin/env python3

# run with:
# visit -cli -nowin -s calculate_fractions.py
import sys
import os
from visit import *
from visit_utils import encoding

#DeleteAllPlots()
dir_name = "frames"
phases = ["Fe2O3_0", "MoO3_0", "Fe2MoO43_0", "liquid_0"]

# Step 1: Open a database (the whole .vtu time series)
dbname="solution-*.vtu database"
OpenDatabase(dbname)

# Step 2: Add plots (using variable "n")
AddPlot("Pseudocolor", phases[0], 1, 1)
DrawPlots()

stride = 20  # Only sample every 10th state
# Get total number of states
n_states = TimeSliderGetNStates()
states = list(range(0, n_states, stride))

print(f"Number of time states: {n_states}\n")
time = [0.0]*len(states)
for idx, state in enumerate(states):
    print(state)
    SetTimeSliderState(state)
    Query("Time")
    time[idx] = GetQueryOutputValue()

fractions = {}
for phase in phases:
    ChangeActivePlotsVar(phase)
    print(f"Processing phase: {phase}...\n")
    fractions[phase] = [0.0]*len(states)
    for idx, state in enumerate(states):
        print(state)
        SetTimeSliderState(state)
        Query("Average Value")
        fractions[phase][idx] = GetQueryOutputValue()

# Step 3: Write results to file
output_file = "phase_fractions.csv"#os.path.join(dir_name, "phase_fractions.csv")
with open(output_file, "w") as f:
    f.write("Time," + ",".join(phases) + "\n")
    for idx, t in enumerate(time):
        f.write(f"{t}," + ",".join(str(fractions[phase][idx]) for phase in phases) + "\n")

DeleteAllPlots()
CloseDatabase(dbname)

sys.exit()
