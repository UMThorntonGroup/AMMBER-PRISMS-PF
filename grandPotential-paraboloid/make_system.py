import ammber.BinarySystems as bi
from ammber.utils import add_to_dict
import json
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database

tdb_file = "NIST-solder.tdb"
#tdb_file = "mmc4.tdb"
db = Database(tdb_file)
elements = ["PB", "BI"]
component = "BI"
solution_component = "PB"
eutectic_temperature = 399.06
eutectic_comp = 0.553

undercooling = 2.0
c0 = 0.553

temperature = eutectic_temperature - undercooling
PbBi_Sys = bi.BinaryIsothermalDiscreteSystem()
PbBi_Sys.fromTDB(db, component, solution_component, temperature)
PbBi_Sys_neareq = PbBi_Sys.resample_near_equilibrium(eutectic_comp)
eq = PbBi_Sys.get_equilibrium(eutectic_comp)
print(eq)

PbBi_Sys_neareq.phases['LIQUID'] = PbBi_Sys.phases['LIQUID'].resample_near_xpoint(eutectic_comp)
PbBi_Fit = bi.BinaryIsothermal2ndOrderSystem()
PbBi_Fit.from_discrete(PbBi_Sys_neareq)

T = temperature  # K
# source https://doi.org/10.1063/1.1730280
D0_PB = 0.432e-4  # m^2/s
b_PB = 20.7
Tm_PB = 600 # K
diffusion_coeff_PB = D0_PB * np.exp(-b_PB*Tm_PB/T)  # m^2/s

# source https://doi.org/10.1007/s11669-021-00868-y
# estimate
diffusion_coeff_LIQ = 1.1e-9#1.1e-9  # m^2/s

# assuming the same for solid phases
diffusion_coeff_HCP_A3 = diffusion_coeff_LIQ #diffusion_coeff_PB  # m^2/s
diffusion_coeff_RHOMBO_A7_1 = diffusion_coeff_LIQ #diffusion_coeff_PB  # m^2/s

# assume 
sigma = 0.1  # J/m^2
mu_int = 2.0e-10 #0.1 *diffusion_coeff_LIQ / sigma
l_int = 2e-7  # m, interfacial length scale

# Load a system file with desired kinetics
#with open("system.json", "r") as f:
#    system = json.load(f)
system = {}

add_to_dict(PbBi_Fit, system, add_templates=True, c0={"HCP_A3": eq["HCP_A3"][0], "RHOMBO_A7_1": eq["RHOMBO_A7_1"][0], "LIQUID": c0}, Vm=2.0e-5)
system["phases"]["HCP_A3"]["D"] = diffusion_coeff_HCP_A3
system["phases"]["RHOMBO_A7_1"]["D"] = diffusion_coeff_RHOMBO_A7_1
system["phases"]["LIQUID"]["D"] = diffusion_coeff_LIQ
system["phases"]["HCP_A3"]["sigma"] = sigma
system["phases"]["RHOMBO_A7_1"]["sigma"] = sigma
system["phases"]["LIQUID"]["sigma"] = sigma
system["phases"]["HCP_A3"]["mu_int"] = mu_int
system["phases"]["RHOMBO_A7_1"]["mu_int"] = mu_int
system["phases"]["LIQUID"]["mu_int"] = mu_int

system["l_int"] = l_int

with open(f"system.json", "w") as f:
    json.dump(system, f, indent=4)

import matplotlib.pyplot as plt
def quad_phase(x, quad_phase):
    return 0.5* quad_phase.kwell * (x - quad_phase.cmin)**2 + quad_phase.fmin
x = np.linspace(0, 1, 100)
plt.plot(x, quad_phase(x, PbBi_Fit.phases['HCP_A3']), label=f"HCP_A3 at {temperature} K", linestyle='--', color='orange')
plt.plot(x, quad_phase(x, PbBi_Fit.phases['RHOMBO_A7_1']), label=f"RHOMBO_A7_1 at {temperature} K", linestyle='--', color='blue')
plt.plot(x, quad_phase(x, PbBi_Fit.phases['LIQUID']), label=f"LIQUID at {temperature} K", linestyle='--', color='green')
print(list(PbBi_Sys.phases.keys()))
print(list(PbBi_Sys_neareq.phases.keys()))
print(PbBi_Sys_neareq.phases['HCP_A3'].xdata.shape, PbBi_Sys_neareq.phases['HCP_A3'].Gdata.shape)
plt.plot(PbBi_Sys.phases['HCP_A3'].xdata, PbBi_Sys.phases['HCP_A3'].Gdata, label=f"HCP_A3 at {temperature} K")
plt.plot(PbBi_Sys.phases['RHOMBO_A7'].xdata, PbBi_Sys.phases['RHOMBO_A7'].Gdata, label=f"RHOMBO_A7 at {temperature} K")
plt.plot(PbBi_Sys.phases['LIQUID'].xdata, PbBi_Sys.phases['LIQUID'].Gdata, label=f"LIQUID at {temperature} K")
plt.plot(PbBi_Sys_neareq.phases['HCP_A3'].xdata, PbBi_Sys_neareq.phases['HCP_A3'].Gdata, label=f"HCP_A3 at {temperature} K", color='orange')
plt.plot(PbBi_Sys_neareq.phases['RHOMBO_A7_1'].xdata, PbBi_Sys_neareq.phases['RHOMBO_A7_1'].Gdata, label=f"RHOMBO_A7 at {temperature} K", color='blue')
plt.plot(PbBi_Sys_neareq.phases['LIQUID'].xdata, PbBi_Sys_neareq.phases['LIQUID'].Gdata, label=f"LIQUID at {temperature} K", color='green')
plt.ylim(top=-10000, bottom = -30000)
plt.show()

#import matplotlib.pyplot as plt
#from pycalphad import binplot
#import pycalphad.variables as v
#
#all_phases = list(db.phases.keys())
#
## Create a matplotlib Figure object and get the active Axes
#fig = plt.figure(figsize=(9,6))
#axes = fig.gca()
#
## Compute the phase diagram and plot it on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
#binplot(db, ['PB', 'BI', 'VA'] , all_phases, {v.X('BI'):(0,1,0.02), v.T: (200, 600, 10), v.P:101325, v.N: 1}, plot_kwargs={'ax': axes})
#
##plt.show()