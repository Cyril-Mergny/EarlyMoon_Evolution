#!/bin/sh

mod_path="modules/"
out_path="out"

gfortran -J$mod_path -I$mod_path -o $out_path/main.out main.f90 write.f90 structure.f90 tides.f90 misc.f90 conduction.f90 orbitevol.f90 


./$out_path/main.out

#cd ../Outputs
#python draw_structure.py
