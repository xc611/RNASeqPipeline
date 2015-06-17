#!/bin/bash
#
#PBS -N rnaseqwf
#PBS -m be

cd $PBS_O_WORKDIR
module load python/3.4.0
module load R/3.1.0

snakemake --jobname 's.{jobid}.{rulename}' --js jobscript.sh -s path2pipeline/rnaseq_wf_v2.py  -k --stats snakemake.stats -T --rerun-incomplete -j 300 --cluster 'qsub {params.batch}'  
