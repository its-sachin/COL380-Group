cd $PBS_O_WORKDIR
module load compiler/gcc/9.1/openmpi/4.0.2

sh compile.sh
sh DataSetup.sh  /scratch/cse/phd/anz198717/TA/COL380/A3/to_students outFolder 4 4
 
