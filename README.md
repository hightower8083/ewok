# ewok
a simple 2D Particle-In-Cell code

This software was my first essay in python-fortran-mpi development.

In 2011 the code was succefully used to model amplification of the radiation scattered by an electron beam in the optical lattice -- so-called optical lattice Raman XFEL [<cite>[1]</cite>], which was simulated using the Lorentz boosted frame.

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
  - phs (analysis of particles phase distributions: choices are x,y,vx,vy,gamma)
  - RETURN key to exit
  - `input time step` and `input limits` fields may be left blank
  
EWOK is not to be used for any scientific research, and its Maxwell is very specific for the case considered in [<cite>[1]</cite>], and it does numerous simplifications of the PIC method (it alo may contain new or old bugs). Some documetation can be found in `./docs/`

\[[1]\] I. A. Andriyash, E. d’Humières, V. T. Tikhonchuk, and Ph. Balcou Phys. Rev. Lett. 109, 244802

[1]:https://doi.org/10.1103/PhysRevLett.109.244802
