## workflow

https://www.lucidchart.com/documents/view/6a805a1f-7d82-43cf-b0a4-4e71faed9ce1

##  Instructions
# Step 0: get the pipeline
   git clone https://github.com/CCBR/RNASeqPipeline

# Step 1: login mode
1/ add a BASH_ENV to  your .bashrc file

  1-1 Create a file called ~/.bash_script_exports.
      The contents of the file should be:
      source /usr/local/Modules/default/init/bash

  1-2 Add to your ~/.bashrc file this line:
      export BASH_ENV=~/.bash_script_exports

2/ allow X11 display
- login with option -X
  example: ssh -X user@biowulf.nih.gov 

# Step 2: data, config and scripts 

3/ cp *.sh path2workingdirectory
4/ cp *.json path2workingdirectory
  
  ## edit the config file & rename it config.json


# Step 3: run pipeline from your working directory

## test and check your configuration (fake run)
5/ bash dryrun.sh

## main call (real run)
6/ qsub -V -l nodes=1:gpfs submit.sh
