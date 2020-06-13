#!/bin/sh

#SBATCH --partition=normal
#SBATCH --mem 15000

module load BEDTools
coverageBed -sorted -b $(echo $1) -a $(echo $2) > $(echo $3)
