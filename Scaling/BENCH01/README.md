In order to run any of these cases, you must first
generate the particle_input.dat file for each case.

To do this, first build the stand-alone executable of par-gen.f90,
e.g. if using the GNU compiler simply type

gfortran par-gen.f90

Then, starting in BENCH1, for each Size*** case where ***

cd Size0001
../a.out --size 1
cd ..

cp Size0008
../a.out --size 8
cd ..

cp Size0027
../a.out --size 27
cd ..

cp Size0064
../a.out --size 64
cd ..

cp Size0216
../a.out --size 216
cd ..

cp Size1000
../a.out --size 1000
cd ..
