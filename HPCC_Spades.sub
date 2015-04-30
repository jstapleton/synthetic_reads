#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=2
#PBS -j oe
#PBS -N Spades_parallel

# Memory needed by the job
#PBS -l mem=10gb

#PBS -t 1-900

### load necessary modules
module load Biopython
module swap GNU GNU/4.7.1
module load SPAdes/3.1.1

cd $PBS_O_WORKDIR
export WORKDIR="/mnt/scratch/${USER}/$PBS_O_WORKDIR/${PBS_JOBNAME}"
mkdir -p $WORKDIR
cp $PBS_O_WORKDIR/splitDicts/splitDict${PBS_ARRAYID}.txt $WORKDIR
cd $WORKDIR

##do the work here
/usr/bin/time -v python $HOME/source/dictReader.py splitDict${PBS_ARRAYID}.txt --runSpades --HPCC

##copy data back
cp -R $WORKDIR ${PBS_O_WORKDIR}/${PBS_JOBID}
rm -r $WORKDIR
