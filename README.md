# What is AMMBER?
A cool software made by X.

# Quick Start Guide
Blah
### Install:
```bash
git clone --recurse-submodules https://github.com/UMThorntonGroup/AMMBER-PRISMS-PF.git
cd AMMBER-PRISMS-PF
```
If you have already cloned the repository without the submodule, you can run
```bash
git submodule update --init --recursive
```
to initialize them.

Next, you will need to install the PRISMS-PF library included as a submodule.
For more information on installing PRISMS-pF and its dependencies, see https://prisms-center.github.io/phaseField/doxygen/install.html.

If you already have the dependencies installed, you can run
```bash
cd phaseField
cmake .
make -j <nprocs>
```
or simply, `make`.

# License:
GNU Lesser General Public License (LGPL). Please see the file [LICENSE](LICENSE) for details.

# Links
[PRISMS Center Homepage](http://www.prisms-center.org/#/home) <br>
[PRISMS-PF Homepage](https://prisms-center.github.io/phaseField/) <br>
[Code Repository](https://github.com/prisms-center/phaseField) <br>
[User Registration Form](http://goo.gl/forms/GXo7Im8p2Y) <br>
[PRISMS-PF User Forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>