#!/bin/bash
#SBATCH --job-name=ag-1
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=dpetrov,hns
#################

nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace -with-dag preprocessing.png
