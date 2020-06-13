"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_coverageBed_adapter: this is an adapter, for running a job per sample in the slurm cluster

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
from statsmodels.distributions.empirical_distribution import ECDF
import subprocess



# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


def extract_number(id):
    '''
    Extract the number id. If there is not, return a 0:
    '''
    try:
        return int(id.split("_")[3])
    except:
        return 0

# def get_job_ids(list):
#     for x in list:



def get_coverageBed_adapter(input_path, gtf_path, coverage_path, output_path, name_user):

    try:
        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # coverage_path = sys.argv[3]
        # output_path = sys.argv[4]
        # temp_path = sys.argv[5]

        # input_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered2.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/random_introns.bed"
        # coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/intron_retention/IR_v5"
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered3.tab"
        # temp_path = "/users/genomics/juanluis/SCLC_cohorts/George/coverageBed/scripts"

        # Load the gtf file
        # logger.info("Loading gtf with the randomizations...")
        # gtf = pd.read_table(gtf_path, delimiter="\t")
        # gtf.columns = ['chr','start','end','id','strand','score']

        # Load the input file
        logger.info("Loading exonizations...")
        exonizations = pd.read_table(input_path, delimiter="\t")
        unique_exonizations = exonizations.drop_duplicates(['Sample_id'])
        unique_sample_ids = unique_exonizations.loc[:,"Sample_id"].tolist()
        dir_path = os.path.dirname(os.path.realpath(__file__))

        # Split the input_path and the gtf_path by sample
        dict_jobs = {}
        for sample in unique_sample_ids:
            # Format the sample
            sample_formatted = sample.rstrip()
            logger.info("Processing "+sample_formatted+"...")
            # Code for SCLC analysis: Remove T's and X's and replace _ by .
            # Only the samples from George or Peifer
            # Remove T's and X's and replace _ by .
            # sample_formatted = sample.replace("T","").replace("X","").replace(".","-")
            # Create a separated input file with only the exonizations associated to sample_formatted
            command1 = "head -1 " + input_path + " > " + output_path + "/input.aux." + sample_formatted + ".tab" + ";" \
                        "awk '{if ($1==\"" + sample_formatted + "\") print }' " + input_path + " >> " + output_path + "/input.aux." + sample_formatted + ".tab"
            os.system(command1)
            # Create an auxiliary script
            command3 = "module load Python/3.5.2; python "+dir_path+"/get_coverageBed.py " \
                       + output_path+"/input.aux."+sample_formatted+".tab " + gtf_path + " " + coverage_path + " " + \
                       output_path + "/get_coverageBed_results." + sample_formatted + ".tab True"
            open_peptides_file = open(output_path + "/aux.sh", "w")
            open_peptides_file.write("#!/bin/sh\n")
            open_peptides_file.write("#SBATCH --partition=normal\n")
            open_peptides_file.write("#SBATCH --mem 1000\n")
            open_peptides_file.write("#SBATCH -e " + output_path + "/" + "get_coverageBed" + "_" + sample_formatted + ".err" + "\n")
            open_peptides_file.write("#SBATCH -o " + output_path + "/" + "get_coverageBed" + "_" + sample_formatted + ".out" + "\n")
            open_peptides_file.write(command3 + ";\n")
            open_peptides_file.close()
            command4 = "sbatch -J "+sample_formatted+"_coverageBed " + output_path + "/aux.sh; sleep 0.5;"
            # os.system(command4)
            job_message = subprocess.check_output(command4, shell=True)
            #Get the job id and store it
            job_id = (str(job_message).rstrip().split(" ")[-1])[:-3]
            dict_jobs[job_id] = 1

        logger.info("Waiting for all the jobs to finished...")
        flag_exit = False
        while(not flag_exit):
            # Initialize the dictionary with the pending jobs in the cluster
            pending_jobs = {}
            os.system("sleep 10")
            p = subprocess.Popen(["squeue","-u", name_user], stdout=subprocess.PIPE)
            # Skip the first line (the header)
            line = p.stdout.readline()
            for line in p.stdout:
                flag_exit = True
                #Get the id of the job
                job_id_aux = str(line).rstrip().split()[1]
                #Save the id of the jobs
                pending_jobs[job_id_aux] = 1
                #If there is any job on the cluster on dict_jobs, break the loop and wait for another 10 seconds
                # to check the status of the jobs in the cluster
                if(job_id_aux in dict_jobs):
                    flag_exit = False
                    break

        logger.info("All jobs finished.")


    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
