gfortran -c -fPIC src_f90/band_dble.f90
gcc -c -I/Users/mgrecu/multiscatter-1.2.10/include/ -fPIC src_f90/multiscatter2_ascii.c
f2py -c -m libScatter src_f90/multiscatter.f90 src_f90/radtran_tau_dble.f src_f90/rosen.f band_dble.o multiscatter2_ascii.o /Users/mgrecu/multiscatter-1.2.10/lib/libmultiscatter.a
