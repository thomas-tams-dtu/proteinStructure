#!/bin/sh
### General options
### –- specify queue --
#BSUB -q gpuv100
### -- set the job Name --
#BSUB -J <complex>
### -- ask for number of cores (default: 1) --
#BSUB -n 4
### -- Select the resources: 1 gpu in exclusive process mode --
#BSUB -gpu "num=1:mode=exclusive_process"
### -- set walltime limit: hh:mm --  maximum 24 hours for GPU-queues right now
#BSUB -W 01:00
# request 5GB of system-memory
#BSUB -R "rusage[mem=5GB]"
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o logs/fold/gpu_<complex>.out
#BSUB -e logs/fold/gpu_<complex>.err
# -- end of LSF options --

# Get localcolabfold env
source /dtu/projects/RFdiffusion/setup.sh
module load colabfold

# Get openmm env
source /zhome/99/d/155947/scratch/miniconda3/bin/activate openmm 

colabfold_batch --templates --num-recycle 1 --num-models 1  data/complexes/<complex>/<complex>_complex.fasta data/complexes/<complex>/
#colabfold_batch --templates --num-recycle 5 --num-models 5  data/complexes/<complex>/<complex>_complex.fasta data/complexes/<complex>/
