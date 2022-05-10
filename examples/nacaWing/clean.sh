#!/bin/bash
rm -f run.pbs.o*
rm -f CAPS_struct/Scratch/capsLock
rm -f CAPS_fluid/Scratch/capsLock
rm -f *.out
rm -f *.txt
rm -f funtofem/*.out
rm -f funtofem/*.txt
rm -f funtofem/*.f5

#remove previous ugrid file
rm -f steady/Flow/*.ugrid

#fun3d outputs
rm -f steady/Flow/core*
rm -f steady/Flow/*.vtk
rm -f steady/Flow/*tec_boundary*
rm -f steady/Flow/*.grid_info
rm -f steady/Flow/*.flow
rm -f steady/Flow/*.forces
rm -f steady/Flow/baseline.flow
rm -f steady/Flow/baseline_volume*
rm -f steady/Flow/baseline_hist*
rm -f steady/Flow/baseline.forces
rm -f steady/Flow/baseline.forces_imaginary
rm -f steady/Flow/*.vtk
rm -f steady/Flow/*.f5
rm -f steady/Flow/*.dat
rm -f steady/Adjoint/*
