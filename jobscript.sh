#!/bin/bash
#rule: {job}
#input: {job.input}
#output: {job.output}
cd $PBS_O_WORKDIR
module load python/3.4.0
module load R/3.1.0

snakemake --snakefile /data/CCBR/dev/RNA-Seq-pipeline/ver2/rnaseq_wf_v2.py \
--force -j 300 \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
> /dev/null && touch "{jobfinished}" || touch "{jobfailed}"
exit 0
