#rm ./t.o
ifort element_table.f90 me_calculate_dipole.f90 test.f90 -c
ifort element_table.o me_calculate_dipole.o test.o -o t.o
./t.o
#python ./draw.py

