#! /bin/bash

source testfunctions.sh

export OMP_NUM_THREADS=1

FASTPM=`dirname $0`/../src/fastpm

set -x

#mpirun -n 4 $FASTPM fastPM20.lua ic || fail

#mpirun -n 4 $FASTPM fastPM20.lua fastpm lineark || fail
#mpirun -n 4 $FASTPM fastPM20.lua fastpm whitenoisek || fail
#mpirun -n 4 $FASTPM fastPM20.lua fastpm || fail


mpirun -n 4 $FASTPM fastPM20.lua fastpm || fail 


