"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

extract_significant_IR: identify the IR events with a minimum of expression

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import numpy as np

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


def extract_significant_IR(input_path, threshold, output_path):

    try:
        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # threshold = int(sys.argv[2])
        # output_path = sys.argv[3]

        # input_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_George_Peifer_Rudin_Yokota.txt"
        # threshold = 1
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant.tab"

        # input_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_Rudin_Normal.txt"
        # threshold = 1
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_Rudin_Normal.tab"

        # 1. Read the input file, extracting all the IR over the threshold
        outfile = open(output_path,"w")
        with open(input_path) as f:
            header = next(f).rstrip().split("\t")
            sample_ids = header[5:]
            Index_pos = header.index("Index")
            outfile.write("Event_id\tSample_id\tTPM\n")
            for line in f:
                #By position, check if the tpm_values are over the threshold
                tokens = line.rstrip().split("\t")
                Index = tokens[Index_pos]
                tpm_values = line.rstrip().split("\t")[5:]
                for i in range(0,len(tpm_values)):
                    if(float(tpm_values[i])>threshold):
                        #Output the line
                        sample_id = sample_ids[i]
                        outfile.write(Index+"\t"+sample_id+"\t"+tpm_values[i]+"\n")

        outfile.close()
        logger.info("Created "+output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
