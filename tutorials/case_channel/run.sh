#! /bin/bash -l
#$ -cwd
#$ -V 
#$ -pe smp.pe 16
#$ -N membr 


   module purge
   module load apps/gcc/openfoam/8


source $foamDotFile


./Allclean
./Allrun
