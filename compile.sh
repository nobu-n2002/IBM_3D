#!/bin/bash
SRC=main_3d.f90

FC=gfortran
# FC=nvfortran

case "$FC" in
#[MEMO] IF "2055 Segmentation fault (core dumped)" APPEARS, ADD "-fno-automatic".
	# gfortran) FC_FLAG='-O3 -fopenmp -fno-automatic -foffload=nvptx-none';;
	gfortran) FC_FLAG='-O3 -fopenmp -fno-automatic';;
	
#[MEMO] IF "Segmentation fault" APPEARS, ADD "-Msave".	
	nvfortran) FC_FLAG='-O3 -mp=multicore -Msave';;


esac

echo 'rm -f *.o *.mod *.exe *.out'
rm -f *.o *.mod *.exe *.out

# In UNIX systems, excutable will be 'a.out' in defaults.
# echo $FC $FC_FLAG $PATH_SRC/$SRC
echo $FC $FC_FLAG $SRC
# $FC $FC_FLAG $PATH_SRC/$SRC
$FC $FC_FLAG $SRC
