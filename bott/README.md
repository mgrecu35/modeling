### Collision-coalescence emulator

This is a python wrapper around the collision-coalescence solver of Bott (1998.) To compile it use:
f2py -c -m bottCC coad1dpy.f 
This will produce a python module that can be imported in python using.


```python
# make the fortran module using f2py if you haven't already
# f2py -c -m bottCC coad1dpy.f 

import bottCC as bcc # import the collision coalescence module

# set parameters 
# rq0 is the radius mode in microns
# xmw is the total water content in g/m3
# nbins is the number of bins (this is hard coded in the fortran code to be 400)
# dt is the time step in seconds 
rq0_in=10.
xmw_in=1.0
nbins=400
dt=1.0
g_out,r_out,dlnr_out = bcc.coad1d_init(rq0_in,xmw_in,nbins,dt)  # this initializes the arrays and sets the initial mass distribution
# r_out is the range of radii in microns
# g_out is the initial mass distribution in g/m3/ln(r) 
# g_out.sum()*dlnr_out is the total water content
# while subroutine coad1d_init analytically sets the initial mass distribution, 
# one can set the initial mass distribution to any distribution they want and 
# reset the initial mass distribution using subroutine set_g_initial
t=0.0
g_out=2*g_out
bcc.set_g_initial(g_out) # this resets the initial mass distribution to 2 times the initial mass distribution

while t<3600:
    g,t = bcc.integrate(nbins,dt,t)

```

Model Name: coasd1d
Model Version: 1.0.0
Licence: GNU GPLv3   

Author: Bott, Andreas

Model Description: 

- coad1d is a model that calculates the evolution of a droplet spectra due to the collision-coalescence process.

- The model follows the bin microphysics approach.

- A full description of the model formulation can be found at  (Bott, A., 1998).



******* DEPENDENCIES ******************************************************************

The model is written in FORTRAN77. A compatible FORTRAN compiler is needed in order to compile and run the code.


******* HOW TO RUN **********************************************************************
1- Compile coad1d.f with the FORTRAN compiler of your preference.

2- Run the *.o executable file generated










 