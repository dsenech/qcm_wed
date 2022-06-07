#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=44
#SBATCH --mem-per-cpu=128M
#SBATCH --constraint=cascade

source ENV/bin/activate
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
OUT_FILE=$(date  +"%H%M%S_%d%m%y").txt
echo "Core,Time (s)" > $OUT_FILE
for i in 44 32 22 16 12 8 6 4 2 1
do
  export OMP_NUM_THREADS=$i
  #export OPENBLAS_NUM_THREADS=$i
  #export CUBACORES=$i
  echo "Threads: $i"
  #time python qcm_tests/debug.py
  \echo -n "$i," >> $OUT_FILE
  \time -f "%e" python debug_10min.py |& tail -n 1 >> $OUT_FILE
done
