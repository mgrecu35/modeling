/Users/mgrecu/homebrew/bin/gfortran -I build/include -c -fPIC -fopenmp keras_def.f90
/Users/mgrecu/homebrew/bin/gfortran -I build/include -I kid -c -fPIC -fopenmp readFile.f90
/Users/mgrecu/homebrew/bin/gfortran -fopenmp -L build/lib -I kid -c -fPIC emulator_interfaceC.f90
/Users/mgrecu/homebrew/bin/gfortran -shared -o keras_interface.so emulator_interfaceC.o readFile.o keras_def.o build/lib/libneural.a kid/libmorisson.so
