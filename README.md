# ewok
a simple 2D Particle-In-Cell code

This software was my first essay in ptyhon-fortran-mpi development.

In 2011 the code was succefully used to model amplification of the radiation scattered by an electron beam in the optical lattice --
so-called [optical lattice Raman XFEL](https://doi.org/10.1103/PhysRevLett.109.244802), which was simulated using the Lorentz boosted frame.

Sweet feature of EWOK is a GUI written with Tkinter package, which allows to setup and launch the simulation.

<p align="center"><img src="https://github.com/hightower8083/ewok/blob/master/docs/gui-snap.png" width="500"/></p>

#### To use:
- start ewokGUI
- select the input file
- choose the output folder (the path should exist)
- run the code

#### To analyse the results:
- go to the output folder
- use the default plotter. The possible choices are:
  - dns (electron density)
  - idns (ion density)
  - fld (electromagnetic potential)
  - estat (electrostatic potential)
  - phs (analysis of particles phase distributions)
  
EWOK is not to be used for any scientific research, and it does numerous simplifications of the PIC method and may contain bugs.

