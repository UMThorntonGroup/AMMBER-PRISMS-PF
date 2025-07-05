import ammber.BinarySystems as bi
from ammber.utils import add_to_dict
import json
import numpy as np
import matplotlib.pyplot as plt

# Load the original free energy data
with open("equations_1.json", "r") as f:
    data = json.load(f)

def make_func(expr):
    # Replace '^' with '**' for Python syntax
    expr = expr.replace('^', '**')
    # Allow log(x) to work
    expr = expr.replace('log(', 'np.log(')
    # Return a lambda that evaluates the expression for a given x (scalar or array)
    return lambda x: eval(expr, {"x": x, "np": np})

objects = []
for entry in data:
    obj = {"Temperature": entry["Temperature"]}
    for key in entry:
        if key != "Temperature":
            obj[key] = make_func(entry[key])
    objects.append(obj)

x = np.linspace(0.01, 0.99, 100)
system_set = []
for obj in objects:
    # Create a BinarySystem object for each entry
    GFE2O3_disc = bi.BinaryIsothermalDiscretePhase("GFE2O3", x, obj["GFE2O3"](x))
    GMoO3_disc = bi.BinaryIsothermalDiscretePhase("GMoO3", x, obj["GMoO3"](x))
    GFE2MoO43_disc = bi.BinaryIsothermalDiscretePhase("GFE2MoO43", x, obj["GFE2MoO43"](x))
    Gliquid_disc = bi.BinaryIsothermalDiscretePhase("Gliquid", x, obj["Gliquid"](x))

    # plt.plot(x, obj["GFE2O3"](x), label=f"GFE2O3 at {obj['Temperature']} K")
    # plt.savefig(f"{obj['Temperature']}.png")

    Gliquid_resample = Gliquid_disc.resample(np.linspace(0.2, 0.3, 100))

    # discrete_system = bi.BinaryIsothermalDiscreteSystem(component="Fe2O3", solution_component="MoO3")
    # discrete_system.phases["GFE2O3"] = GFE2O3_disc
    # discrete_system.phases["GMoO3"] = GMoO3_disc
    # discrete_system.phases["GFE2MoO43"] = GFE2MoO43_disc
    # discrete_system.phases["Gliquid"] = Gliquid_resample

    # fit_system = bi.BinaryIsothermal2ndOrderSystem()
    # fit_system.from_discrete(discrete_system, kwellmax=1.0e11)

    # alternatively
    liq_fit = bi.BinaryIsothermal2ndOrderPhase("Gliquid")
    liq_fit.fit_phase(Gliquid_resample.xdata, Gliquid_resample.Gdata, kwellmax=1.0e11)
    big_k = 5000000.0
    phases = {
        "GFE2O3": bi.BinaryIsothermal2ndOrderPhase("GFE2O3", obj["GFE2O3"](1.0), big_k, 1.0 ),
        "GMoO3": bi.BinaryIsothermal2ndOrderPhase("GMoO3", obj["GMoO3"](0.0), big_k, 0.0 ),
        "GFE2MoO43": bi.BinaryIsothermal2ndOrderPhase("GFE2MoO43", obj["GFE2MoO43"](0.25), big_k, 0.25 ),
        "Gliquid": liq_fit }
    
    fit_system = bi.BinaryIsothermal2ndOrderSystem(phases=phases)
    fit_system.component="Fe2O3"
    fit_system.solution_component="MoO3"

    # For some reson, this breaks if put outside the loop
    # Load a system file with desired kinetics
    with open("system.json", "r") as f:
        base_system = json.load(f)

    output = {}
    output = base_system.copy()

    add_to_dict(fit_system, output, add_templates=False, c0={"GFE2O3": 1.0, "GMoO3": 0.0, "GFE2MoO43": 0.25, "Gliquid": 0.25}, Vm=3.1e-5)
    system_set.append({obj['Temperature']: output})

for system in system_set:
    for temp, output in system.items():
        with open(f"system_{temp}.json", "w") as f:
            json.dump(output, f, indent=4)