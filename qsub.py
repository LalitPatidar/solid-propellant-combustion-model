# Program for creating .pbs file for each gaussian input file i.e. .com file
				  
import os
#Path1 = "./M062X-minima" 
#Path2 = "./M062X-minima/pbs"

dirs = ['./1atm','./5atm','./20atm','./50atm','./90atm']

bash_file_name = 'bash.sh'
version = 'for-Chapter-6/All-final-simulations/Tinit-298K-test'

for item in dirs:
    with open(os.path.join(item,item[2:]+'.pbs'), "w") as File:	
        File.write('#PBS -l nodes=1:ppn=1\n')
        File.write('#PBS -l walltime=48:00:00\n')
        File.write('#PBS -j oe\n')
        File.write('\n')
        File.write('# change the current working directory to the directory where\n')
        File.write('# the input deck foo.com can be found\n')
        File.write('cd /gpfs/group/umt/default/HMX/combustion-model/HMX-final-combustion/'+version+item[1:]+' " "\n')
        File.write('echo "Starting job on `hostname` at `date`"\n')
        File.write('echo " "\n')
        File.write('\n')
        File.write('# load the g09 module\n')
        File.write('module load python\n')
        File.write('\n')
        File.write('source activate /storage/home/lkp5147/work/sw/cantera-cython\n')
        File.write('\n')
        File.write('# start g09 with input deck foo.com\n')
        File.write('python run.py')
        File.write('\n')       
                
        File.write('\n')
        File.write('echo " "\n')
        File.write('echo "Completing job on `hostname` at `date`"\n')
        File.write('echo " "\n')

with open(bash_file_name, "w") as File:
    File.write('#!/bin/bash\n')
    File.write('\n')
    
    for item in dirs:
        dos2unix = 'dos2unix '+ item[2:]+'.pbs'
        #qsub = 'qsub -A umt_a_g_sc_default '+ item[2:]+'.pbs'
        qsub = 'qsub -A umt_b_g_sc_default '+ item[2:]+'.pbs'
        qsub = 'qsub -A open '+ item[2:]+'.pbs'
        File.write('cd ' + item[2:] + '/\n')
        File.write(dos2unix)
        File.write('\n')
        File.write(qsub)
        File.write('\n')
        File.write('cd ..\n')
         
    File.write('\n')
    File.write('exit 0\n')
    File.write('\n')      
   
