"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_IR: from our list opf IR, we are gonna remove those ones that appears in the Rudin data.
We will exclude a junction if we see that there is at least a sample in any data that has a
minimum expression
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


def filter_IR(exonizations_path, rudin_path, output_path):

    try:
        logger.info("Starting execution")

        # exonizations_path = sys.argv[1]
        # rudin_path = sys.argv[2]
        # output_path = sys.argv[3]
        #
        # exonizations_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes.tab"
        # rudin_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_Rudin_Normal.tab"
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered.tab"

        #Load the exonizations
        exonizations = pd.read_table(exonizations_path, delimiter="\t")
        #Load the Rudin exons
        rudin = pd.read_table(rudin_path, delimiter="\t")

        #Merge the files and extract the exons that are not in the other files
        merge1 = exonizations.merge(rudin, on='Event_id', how='left', suffixes=('', '_y'))
        merge1_f = merge1[~merge1.Sample_id_y.notnull()]

        #Remove the extra columns
        n_col = len(exonizations.columns)
        df = merge1_f.iloc[:, :n_col].copy()
        df.to_csv(output_path, sep="\t", index=False, header=True)
        logger.info("Saved "+output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
