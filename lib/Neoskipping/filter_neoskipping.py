"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_neoskipping: from our list of neoskipping, we are gonna remove those ones that appears in the Rudin
and Intropolis data. We will exclude a junction if we see that there is at least a sample in any data that has a
minimum number of reads
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re

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


def filter_neoskipping(neoskipping_path, rudin_path, intropolis_path, output_path, flag_Rudin):

    try:
        logger.info("Starting execution")

        if(flag_Rudin):

            #Load the neoskipping
            neoskipping = pd.read_table(neoskipping_path, delimiter="\t")
            #Load the Rudin junction
            rudin = pd.read_table(rudin_path, delimiter="\t")
            #Load the Intropolis junctions
            intropolis = pd.read_table(intropolis_path, delimiter="\t")

            #Merge the files and extract the exons that are not in the other files
            merge1 = neoskipping.merge(rudin, on='Neoskipping_junction', how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Sample_id_y.notnull()]
            merge2 = merge1_f.merge(intropolis, on='Neoskipping_junction', how='left', suffixes=('', '_z'))
            merge2_f = merge2[~merge2.Sample_id_z.notnull()]

            #Remove the extra columns
            df = merge2_f.iloc[:, 0:7].copy()
            df.to_csv(output_path, sep="\t", index=False, header=True)
            logger.info("Saved "+output_path)
            logger.info("Done. Exiting program.")

        else:

            # Load the neoskipping
            neoskipping = pd.read_table(neoskipping_path, delimiter="\t")
            # Load the Intropolis junctions
            intropolis = pd.read_table(intropolis_path, delimiter="\t")

            # Merge the files and extract the exons that are not in the other files
            merge1 = neoskipping.merge(intropolis, on='Neoskipping_junction', how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Sample_id_y.notnull()]

            # Remove the extra columns
            df = merge1_f.iloc[:, 0:7].copy()
            # df = merge1_f.copy()
            df.to_csv(output_path, sep="\t", index=False, header=True)
            logger.info("Saved " + output_path)
            logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)