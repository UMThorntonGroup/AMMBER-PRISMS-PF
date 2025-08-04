import ammber.BinarySystems as bi
from ammber.utils import add_to_dict
import json
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database

tdb_file = "AlCu.TDB"
db = Database(tdb_file)
elements = ["AL", "CU"]
component = "CU"
solution_component = "AL"
eutectic_temperature = 819.0  # K
eutectic_comp = 0.177

undercooling = 8.0
c0 = 0.177

temperature = eutectic_temperature - undercooling
AlCu_Sys = bi.BinaryIsothermalDiscreteSystem()
AlCu_Sys.fromTDB(db, component, solution_component, temperature)
AlCu_Sys_neareq = AlCu_Sys.resample_near_equilibrium(eutectic_comp)
eq = AlCu_Sys.get_equilibrium(eutectic_comp)
print(eq)

AlCu_Sys_neareq.phases['LIQUID'] = AlCu_Sys.phases['LIQUID'].resample_near_xpoint(eutectic_comp)
AlCu_Fit = bi.BinaryIsothermal2ndOrderSystem()
AlCu_Fit.from_discrete(AlCu_Sys_neareq)

#T = temperature  # K
## source https://doi.org/10.1063/1.1730280
#D0_PB = 0.432e-4  # m^2/s
#b_PB = 20.7
#Tm_PB = 600 # K
#diffusion_coeff_PB = D0_PB * np.exp(-b_PB*Tm_PB/T)  # m^2/s

# source https://doi.org/10.1007/s11669-021-00868-y
# estimate
diffusion_coeff_LIQ = 1.1e-9#1.1e-9  # m^2/s

# assuming the same for solid phases
diffusion_coeff_FCC_A1_0 = diffusion_coeff_LIQ/100 # m^2/s
diffusion_coeff_AL2CU_1 = diffusion_coeff_LIQ /100 # m^2/s

# assume 
sigma = 0.2  # J/m^2
mu_int = 2.0e-9 # 0.1 * diffusion_coeff_LIQ / sigma
l_int = 1e-8  # m, interfacial length scale

# Load a system file with desired kinetics
#with open("system.json", "r") as f:
#    system = json.load(f)
system = {}

add_to_dict(AlCu_Fit, system, add_templates=True, c0={"FCC_A1_0": eq["FCC_A1_0"][0], "AL2CU_1": eq["AL2CU_1"][0], "LIQUID": c0}, Vm=1.0e-5)
system["phases"]["FCC_A1_0"]["D"] = diffusion_coeff_FCC_A1_0
system["phases"]["AL2CU_1"]["D"] = diffusion_coeff_AL2CU_1
system["phases"]["LIQUID"]["D"] = diffusion_coeff_LIQ
system["phases"]["FCC_A1_0"]["sigma"] = sigma
system["phases"]["AL2CU_1"]["sigma"] = sigma
system["phases"]["LIQUID"]["sigma"] = sigma
system["phases"]["FCC_A1_0"]["mu_int"] = mu_int
system["phases"]["AL2CU_1"]["mu_int"] = mu_int
system["phases"]["LIQUID"]["mu_int"] = mu_int

system["l_int"] = l_int

with open(f"system.json", "w") as f:
    json.dump(system, f, indent=4)

import matplotlib.pyplot as plt
def quad_phase(x, quad_phase):
    return 0.5* quad_phase.kwell * (x - quad_phase.cmin)**2 + quad_phase.fmin
x = np.linspace(0, 1, 100)
plt.plot(x, quad_phase(x, AlCu_Fit.phases['FCC_A1_0']), label=f"FCC_A1_0 at {temperature} K", linestyle='--', color='orange')
plt.plot(x, quad_phase(x, AlCu_Fit.phases['AL2CU_1']), label=f"AL2CU_1 at {temperature} K", linestyle='--', color='blue')
plt.plot(x, quad_phase(x, AlCu_Fit.phases['LIQUID']), label=f"LIQUID at {temperature} K", linestyle='--', color='green')
print(list(AlCu_Sys.phases.keys()))
print(list(AlCu_Sys_neareq.phases.keys()))
plt.plot(AlCu_Sys.phases['FCC_A1'].xdata, AlCu_Sys.phases['FCC_A1'].Gdata, label=f"FCC_A1 at {temperature} K", color='orange')
plt.plot(AlCu_Sys.phases['AL2CU'].xdata, AlCu_Sys.phases['AL2CU'].Gdata, label=f"AL2CU at {temperature} K", color='blue')
plt.plot(AlCu_Sys.phases['LIQUID'].xdata, AlCu_Sys.phases['LIQUID'].Gdata, label=f"LIQUID at {temperature} K", color='green')
plt.plot(AlCu_Sys_neareq.phases['FCC_A1_0'].xdata, AlCu_Sys_neareq.phases['FCC_A1_0'].Gdata, label=f"FCC_A1_0 at {temperature} K", color='orange', marker='o', markersize=2)
plt.plot(AlCu_Sys_neareq.phases['AL2CU_1'].xdata, AlCu_Sys_neareq.phases['AL2CU_1'].Gdata, label=f"AL2CU_1 at {temperature} K", color='blue', marker='o', markersize=2)
plt.plot(AlCu_Sys_neareq.phases['LIQUID'].xdata, AlCu_Sys_neareq.phases['LIQUID'].Gdata, label=f"LIQUID at {temperature} K", color='green', marker='o', markersize=2)
plt.ylim(top=-25000, bottom = -50000)
plt.xlim(0,0.34)
plt.savefig("AlCu_energy.png", dpi=300, bbox_inches='tight')

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