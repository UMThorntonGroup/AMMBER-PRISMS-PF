# What is AMMBER?
The AI-assisted Microstructure Model BuildER (AMMBER) is an ongoing project at the University of Michigan.

AMMBER-PRISMS-PF is a phase-field simulation code.

Phase-field models, which incorporate thermodynamic and kinetic data from atomistic calculations and experiments, have become a key computational tool for understanding microstructural evolution and providing a path to control and optimize morphologies and topologies of structures from nanoscale to microscales. However, due to the complexity of interactions between multiple species, these models are difficult to parameterize. In this project, we developed algorithms and software that automate and optimize the selection of thermodynamic and kinetic parameters for phase-field simulations of microstructure evolution in multicomponent systems.

Presently, the framework consists of two modules: [AMMBER_python](https://github.com/UMThorntonGroup/AMMBER_python), which is used to extract phase-field usable free energies from general data sources, and [AMMBER-PRISMS-PF](https://github.com/UMThorntonGroup/AMMBER-PRISMS-PF), which provides an open-source suite of multi-component, multi-phase-field model implementations with a simple, flexible interface for defining a system of thermodynamic and kinetic parameters.

AMMBER-PRISMS-PF Features:
- Open-source
- Multi-component
- Multi-phase
- Simple, flexible input file
- Automatic parameter selection
- High-performance Code
- Adaptive Mesh Refinement (AMR)

# Quick Start Guide

### Install:
AMMBER-PRISMS-PF can be installed on Linux and MacOS. <br>
1. Install [PRISMS-PF](https://github.com/prisms-center/phaseField). For more information on installing PRISMS-PF and its dependencies, see [PRISMS-PF installation instructions](https://prisms-center.github.io/phaseField/doxygen/install.html).

2. Install AMMBER-PRISMS-PF 
First, set an installation path
```bash
# In your .rc or .profile file
export AMMBER_DIR='/path/to/where/to/install'
```
In the terminal, clone this repository, navigate inside, and install with cmake.
```bash
git clone https://github.com/UMThorntonGroup/AMMBER-PRISMS-PF.git
cd AMMBER-PRISMS-PF
cmake -B build
cmake --install build --prefix=$AMMBER_DIR
```

#### Recommended: Install [AMMBER_python](https://github.com/UMThorntonGroup/AMMBER_python)
```bash
pip install ammber
```

### Running an application
Each application in this suite has a more detailed README, explaining how to use each model. PRISMS-PF also has [documentation](https://prisms-center.github.io/phaseField/doxygen/app_structure.html) explaining the requirements of an app.
To run an application without any modifications, you just need to navigate to the application directory and compile first. For example:
```bash
cd examples/grand-potential/paraboloid
cmake -B release -DCMAKE_BUILD_TYPE=Release
cmake --build release
```
Next, you can run the simulation in parallel using
```bash
mpirun -n <nprocs> release/main
```
(nprocs=8 for most desktops) or just `./main` for serial.

### Visualization:

Output of the fields is in standard vtk
format (parallel:*.pvtu, serial:*.vtu files) which can be visualized with the
following open source applications:

1. VisIt (https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
2. Paraview (http://www.paraview.org/download/)

# License:
GNU Lesser General Public License (LGPL). Please see [LICENSE](LICENSE) for details.

# Acknowledgement 
This project is made possible by funding from the National Science Foundation (NSF) Award No. OAC-2209423

# Links
[AMMBER-PRISMS-PF Repository](https://github.com/UMThorntonGroup/AMMBER-PRISMS-PF) <br>
[AMMBER_python Repository](https://github.com/UMThorntonGroup/AMMBER_python) <br>
[PRISMS-PF Homepage](https://prisms-center.github.io/phaseField/) <br>
[PRISMS-PF Repository](https://github.com/prisms-center/phaseField) <br>
[PRISMS-PF User Forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>