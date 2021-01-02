"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_significant_exonizations: given the table of the neoskipping events with the reads counts,
get those that are over a threshold

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re, csv

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

def get_significant_exonizations(exonization_file, threshold, fold, output_file):

    try:
        logger.info("Starting execution")

        # Load the reads with the exonizations
        exonizations_reads = pd.read_table(exonization_file, delimiter="\t")
        fold_flag = False

        if(fold!=0):
            fold_flag = True

        # Iterate the file and extract all the samples above a certain threshold
        logger.info("Filtering reads...")
        if(not fold_flag):
            df_filtered = exonizations_reads[exonizations_reads['mean_skipped_ReadCounts'] >= threshold]
        else:
            df_filtered = exonizations_reads[exonizations_reads['mean_skipped_ReadCounts'] >= threshold and
                                             exonizations_reads['Fold'] >= fold]
        df_filtered.to_csv(output_file,index=False,header=False,sep ='\t',quoting=csv.QUOTE_NONE)
        logger.info("Saved "+ output_file)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
