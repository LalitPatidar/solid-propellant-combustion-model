#!/bin/bash

cd 1atm/
dos2unix 1atm.pbs
qsub -A open 1atm.pbs
cd ..
cd 5atm/
dos2unix 5atm.pbs
qsub -A open 5atm.pbs
cd ..
cd 20atm/
dos2unix 20atm.pbs
qsub -A open 20atm.pbs
cd ..
cd 50atm/
dos2unix 50atm.pbs
qsub -A open 50atm.pbs
cd ..
cd 90atm/
dos2unix 90atm.pbs
qsub -A open 90atm.pbs
cd ..

exit 0

