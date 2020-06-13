"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

run_netMHC-4.0: run netMHC-4.0 on each sample
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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


def run_netMHC_classI_slurm_part1(input_list_path, HLAclass_path, HLAtypes_path, input_sequence_pieces_path, output_netMHC_path,
                                  output_peptides_path,output_peptides_all_path,output_peptides_path2, output_peptides_all_path2,
                                  output_list_path,netMHC_path):
    try:
        logger.info("Starting execution")

        # Load the list of accepted HLA types
        logger.info("Load the list of accepted HLA types")
        HLA_accepted_types = set()
        with open(HLAtypes_path) as f:
            for line in f:
                tokens = line.rstrip()
                HLA_accepted_types.add(tokens)

        # Assign to each sample their corresponding HLA types according to the results with seq2HLA
        logger.info("Assigning the prediction for the HLA types to each sample")
        HLA_samples = {}
        with open(HLAclass_path) as f:
            next(f)
            cont = 0
            for line in f:
                cont += 1
                tokens = line.rstrip().split("\t")
                # Check if the HLA_types are significant and if that type exists
                aux = "HLA-" + tokens[1].replace("'", "").replace("*", "").replace(":", "")
                # A1 class
                if (float(tokens[2]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)
                aux = "HLA-" + tokens[3].replace("'", "").replace("*", "").replace(":", "")
                # A2 class
                if (float(tokens[4]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)
                aux = "HLA-" + tokens[5].replace("'", "").replace("*", "").replace(":", "")
                # B1 class
                if (float(tokens[6]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)
                aux = "HLA-" + tokens[7].replace("'", "").replace("*", "").replace(":", "")
                # B2 class
                if (float(tokens[8]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)
                aux = "HLA-" + tokens[9].replace("'", "").replace("*", "").replace(":", "")
                # C1 class
                if (float(tokens[10]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)
                aux = "HLA-" + tokens[11].replace("'", "").replace("*", "").replace(":", "")
                # C2 class
                if (float(tokens[12]) <= 0.05 and aux in HLA_accepted_types):
                    if (tokens[0] not in HLA_samples):
                        HLA_samples[tokens[0]] = [aux]
                    else:
                        if (aux not in HLA_samples[tokens[0]]):
                            HLA_samples[tokens[0]].append(aux)

        # Go over the input file, running netMHC
        logger.info("Processing samples for running netMHC")
        # open_peptides_file = open(output_peptides_path, "w")
        # open_peptides_file.write("Exonization\tSample_ID\tHLA_type\tPredicted_neoantigen\tBinding_affinity\n")
        # open_peptides_file2 = open(output_peptides_path2, "w")
        # open_peptides_file2.write("Exonization\tSample_ID\tHLA_type\tPredicted_neoantigen\tBinding_affinity\n")
        # open_peptides_all_file = open(output_peptides_all_path, "w")
        # open_peptides_all_file.write("Exonization\tSample_ID\tHLA_type\tPredicted_neoantigen\tBinding_affinity\n")
        # open_peptides_all_file2 = open(output_peptides_all_path2, "w")
        # open_peptides_all_file2.write("Exonization\tSample_ID\tHLA_type\tPredicted_neoantigen\tBinding_affinity\n")
        status_neoantigen = []
        path1 = "/".join(output_peptides_path.split("/")[:-1])
        cont = 0
        with open(input_list_path) as f:
            next(f)
            for line in f:
                cont += 1
                # if(cont==3):
                #     break
                logger.info("Index: " + str(cont))
                tokens1 = line.rstrip().split("\t")
                index = tokens1[0]
                exonization = tokens1[3]
                sample = tokens1[2].rstrip()
                results_by_exon = []
                # Get the HLA types associated. Run netMHC for each HLA type
                if (sample in HLA_samples):
                    HLA_types = HLA_samples[sample]
                    cont2 = 0
                    for x in HLA_types:
                        logger.info(
                            "Running job HLA-type: " + output_netMHC_path + "/" + index + "_" + x + ".out")
                        cont2 += 1
                        command1 = netMHC_path + " -a " + x + " -f " + input_sequence_pieces_path + "/" + index + \
                                   ".fa > " + output_netMHC_path + "/" + index + "_" + x + ".out"
                        # Output this to an auxiliary script
                        open_peptides_file = open(path1 + "/aux.sh", "w")
                        open_peptides_file.write("#!/bin/sh\n")
                        open_peptides_file.write("#SBATCH --partition=lowmem\n")
                        open_peptides_file.write("#SBATCH --mem 2000\n")
                        open_peptides_file.write(
                            "#SBATCH -e " + path1 + "/" + index + "_" + x + ".err" + "\n")
                        open_peptides_file.write(
                            "#SBATCH -o " + path1 + "/" + index + "_" + x + ".out" + "\n")
                        open_peptides_file.write(command1 + ";\n")
                        open_peptides_file.close()
                        command2 = "module load Perl; sbatch -J " + index + "_" + x + " " + path1 + "/aux.sh; sleep 0.5"
                        os.system(command2)


                else:
                    pass

        logger.info("When all the jobs have finished, run part2.")
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
