export FC=/Users/mgrecu/homebrew/bin/gfortran
export F90=/Users/mgrecu/homebrew/bin/gfortran
export FFLAGS="-fPIC -g -ffpe-trap=invalid,zero,overflow"
gfortran -c -fPIC -g src/typeKind.f90 
cpp -DDEF_NZ=120 -DDEF_NX=210 src/parameters.f90 > src/parameters_cpp.f90
gfortran -c -fPIC -g -ffpe-trap=invalid src/parameters_cpp.f90
gfortran -c -fPIC -g -ffpe-trap=invalid src/class_species.f90
gfortran -c -fPIC -g -fcheck=all -ffpe-trap=invalid src/column_variables.f90
gfortran -c -fPIC -g -fcheck=all -ffpe-trap=invalid src/header_data.f90
gfortran -c -fPIC -g -fcheck=all -ffpe-trap=invalid src/switches.f90
gfortran -c -fPIC -g -fcheck=all -ffpe-trap=invalid src/switches_bin.f90
gfortran -c -fPIC -g -fcheck=all -ffpe-trap=invalid src/physconst.f90
gfortran -c -fPIC -g -ffpe-trap=invalid,zero,overflow -O2 -fcheck=all -finit-real=nan -fopenmp src/module_mp_morr_two_moment_nodg.f90
#gfortran -c -fPIC -g -ffpe-trap=invalid,zero,overflow -O2 -finit-real=nan src/mphys_morr_two_moment_nodg.f90
#gfortran -shared -o libmorisson.so  src/morrison_py.f90  *.o
gfortran -shared -o libmorisson.so   *.o
