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
    Fe2O3_disc = bi.BinaryIsothermalDiscretePhase("Fe2O3", x, obj["GFE2O3"](x))
    MoO3_disc = bi.BinaryIsothermalDiscretePhase("MoO3", x, obj["GMoO3"](x))
    Fe2MoO43_disc = bi.BinaryIsothermalDiscretePhase("Fe2MoO43", x, obj["GFE2MoO43"](x))
    liquid_disc = bi.BinaryIsothermalDiscretePhase("liquid", x, obj["Gliquid"](x))

    # plt.plot(x, obj["Fe2O3"](x), label=f"Fe2O3 at {obj['Temperature']} K")
    # plt.savefig(f"{obj['Temperature']}.png")

    liquid_resample = liquid_disc.resample(np.linspace(0.2, 0.3, 100))

    # discrete_system = bi.BinaryIsothermalDiscreteSystem(component="Fe2O3", solution_component="MoO3")
    # discrete_system.phases["Fe2O3"] = Fe2O3_disc
    # discrete_system.phases["MoO3"] = MoO3_disc
    # discrete_system.phases["Fe2MoO43"] = Fe2MoO43_disc
    # discrete_system.phases["liquid"] = liquid_resample

    # fit_system = bi.BinaryIsothermal2ndOrderSystem()
    # fit_system.from_discrete(discrete_system, kwellmax=1.0e11)

    # alternatively
    liq_fit = bi.BinaryIsothermal2ndOrderPhase("liquid")
    liq_fit.fit_phase(liquid_resample.xdata, liquid_resample.Gdata, kwellmax=1.0e11)
    big_k = 5000000.0
    phases = {
        "Fe2O3": bi.BinaryIsothermal2ndOrderPhase("Fe2O3", obj["GFE2O3"](1.0), big_k, 1.0 ),
        "MoO3": bi.BinaryIsothermal2ndOrderPhase("MoO3", obj["GMoO3"](0.0), big_k, 0.0 ),
        "Fe2MoO43": bi.BinaryIsothermal2ndOrderPhase("Fe2MoO43", obj["GFE2MoO43"](0.25), big_k, 0.25 ),
        "liquid": liq_fit }
    
    fit_system = bi.BinaryIsothermal2ndOrderSystem(phases=phases)
    fit_system.component="Fe2O3"
    fit_system.solution_component="MoO3"

    # For some reson, this breaks if put outside the loop
    # Load a system file with desired kinetics
    with open("system.json", "r") as f:
        base_system = json.load(f)

    output = {}
    output = base_system.copy()

    add_to_dict(fit_system, output, add_templates=False, c0={"Fe2O3": 1.0, "MoO3": 0.0, "Fe2MoO43": 0.25, "liquid": 0.25}, Vm=3.1e-5)
    system_set.append({obj['Temperature']: output})

for system in system_set:
    for temp, output in system.items():
        with open(f"system_{temp}.json", "w") as f:
            json.dump(output, f, indent=4)