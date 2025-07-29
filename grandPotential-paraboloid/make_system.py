import ammber.BinarySystems as bi
from ammber.utils import add_to_dict
import json
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database

tdb_file = "NIST-solder.tdb"
db = Database(tdb_file)
elements = ["PB", "BI"]
component = "BI"
solution_component = "PB"
eutectic_temperature = 400.0
eutectic_comp = 0.56

undercooling = 20.0
c0 = 0.56

temperature = eutectic_temperature - undercooling
PbBi_Sys = bi.BinaryIsothermalDiscreteSystem()
PbBi_Sys.fromTDB(db, component, solution_component, temperature)
PbBi_Sys_neareq = PbBi_Sys.resample_near_equilibrium(eutectic_comp)

PbBi_Sys_neareq.phases['LIQUID'] = PbBi_Sys.phases['LIQUID'].resample_near_xpoint(eutectic_comp)
PbBi_Fit = bi.BinaryIsothermal2ndOrderSystem()
PbBi_Fit.from_discrete(PbBi_Sys_neareq)

D0 = 1e-9  # m^2/s
E0 = 1e5  # J/mol
k = 8.314  # J/(mol*K)
T = temperature  # K
diffusion_coeff_0 = D0 * np.exp(E0/(k*T))  # m^2/s
# Load a system file with desired kinetics
#with open("system.json", "r") as f:
#    system = json.load(f)
system = {}

add_to_dict(PbBi_Fit, system, add_templates=True, c0={"HCP_A3": 0.2, "RHOMBO_A7_1": 1.0, "LIQUID": c0}, Vm=3.1e-5)

with open(f"system.json", "w") as f:
    json.dump(system, f, indent=4)


import matplotlib.pyplot as plt
from pycalphad import binplot
import pycalphad.variables as v

all_phases = list(db.phases.keys())

# Create a matplotlib Figure object and get the active Axes
fig = plt.figure(figsize=(9,6))
axes = fig.gca()

# Compute the phase diagram and plot it on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
binplot(db, ['PB', 'BI', 'VA'] , all_phases, {v.X('BI'):(0,1,0.02), v.T: (200, 600, 10), v.P:101325, v.N: 1}, plot_kwargs={'ax': axes})

#plt.show()