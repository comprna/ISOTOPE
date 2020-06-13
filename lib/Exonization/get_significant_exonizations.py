"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_significant_exonizations: given the table of the exonizations with the reads counts,
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

def get_significant_exonizations(exonization_file, threshold, output_file):

    try:
        logger.info("Starting execution")

        # exonization_file = sys.argv[1]
        # threshold = int(sys.argv[2])
        # output_file = sys.argv[3]

        # exonization_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_exonized_junctions_reads_repeatitions.tab"
        # threshold = 5
        # output_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample.tab"

        # exonization_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/new_exonized_junctions_Rudin_normal_reads_repeatitions.tab"
        # threshold = 10
        # output_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/exonizations_by_sample_Rudin_normal.tab"

        # exonization_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_exonized_Intropolis_repeatitions.tab"
        # threshold = 10
        # output_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_Intropolis_by_sample.tab"

        # Load the reads with the exonizations
        exonizations_reads = pd.read_table(exonization_file, delimiter="\t")

        # Iterate the file and extract all the samples above a certain threshold
        logger.info("Obtaining significant reads...")
        columns = list(exonizations_reads.columns.values)
        f = open(output_file, 'w')
        f.write("Sample_id\tNew_exon\tExon_length\tGene\tJunction_id3\tJunction_id4\tsplice_site5\tsplice_site3\tRepeats\tReadCounts\n")
        logger.info(str(len(exonizations_reads.index)))
        for index, row in exonizations_reads.iterrows():
            # logger.info(str(index))
            for i in range(8,len(columns)):
                if(row.iloc[i]>=threshold):
                    #Output this value
                    f.write(row.index[i]+"\t"+row.iloc[3]+"\t"+str(row.iloc[4])+"\t"+row.iloc[0]+"\t"+
                            row.iloc[1]+"\t"+row.iloc[2]+"\t"+row.iloc[5]+"\t"+row.iloc[6]+"\t"+row.iloc[7]+"\t"+
                            str(row.iloc[i])+"\n")

        f.close()
        logger.info("Saved "+ output_file)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
