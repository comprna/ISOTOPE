#!/bin/sh

#SBATCH --mem 15000
#SBATCH --partition=normal

module load BEDTools
coverageBed -sorted -b $(echo $1) -a $(echo $2) > $(echo $3)
