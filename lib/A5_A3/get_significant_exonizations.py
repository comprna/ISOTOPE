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

        # Load the reads with the exonizations
        exonizations_reads = pd.read_table(exonization_file, delimiter="\t")

        # Iterate the file and extract all the samples above a certain threshold
        logger.info("Obtaining significant reads...")
        columns = list(exonizations_reads.columns.values)
        f = open(output_file, 'w')
        f.write("Sample_id\tGene\tCanonical_Junction_id\tAlt_Junction_id\tCanonical_Exon\tAlt_Exon_id\tstrand\tOffset\tRepeats\tNew_Exon_length\tmotif\tsplice_site_type\tReadCounts\n")
        logger.info(str(len(exonizations_reads.index)))
        for index, row in exonizations_reads.iterrows():
            # logger.info(str(index))
            for i in range(11,len(columns)):
                if(row.iloc[i]>=threshold):
                    #Output this value
                    f.write(row.index[i]+"\t"+row.iloc[0]+"\t"+row.iloc[1]+"\t"+row.iloc[2]+"\t"+
                            row.iloc[3]+"\t"+row.iloc[4]+"\t"+row.iloc[5]+"\t"+str(row.iloc[6])+"\t"+
                            str(row.iloc[7])+"\t"+str(row.iloc[8])+"\t"+str(row.iloc[9])+"\t"+str(row.iloc[10])+"\t"+
                            str(row.iloc[i])+"\n")

        f.close()
        logger.info("Saved "+ output_file)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
