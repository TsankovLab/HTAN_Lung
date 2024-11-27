#!/bin/bash
#BSUB -J visium_cnv
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 1
#BSUB -W 12:00
#BSUB -gpu num=1
#BSUB -R v100
#BSUB -R rusage[mem=64000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10
ml cuda/11.7.0;
ml cudnn/8.9.5-11;
ml proxies;
ml java/11.0.2;

source activate cnv

cd /sc/arion/projects/Tsankov_Normal_Lung/users/tp53_resubmission/final_codes/

python run_infercnv.py
