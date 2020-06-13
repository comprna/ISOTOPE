"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_reads_exonizations: given the list with the possible exonizations, get the reads associate to each of them

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re, csv
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

def get_reads_exonizations(exonization_file, readCounts_file, output_file):

    try:
        logger.info("Starting execution")

        # exonization_file = sys.argv[1]
        # readCounts_file = sys.argv[2]
        # output_file = sys.argv[3]

        # exonization_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/new_exonized_junctions.tab"
        # readCounts_file = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
        # output_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/new_exonized_junctions_Rudin_normal_reads.tab"

        # exonization_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_exonized_junctions.tab"
        # readCounts_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/readCounts_George_Peifer_Rudin_Yokota.tab"
        # output_file = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_exonized_junctions_reads.tab"

        # Load the readCounts
        junction_reads, junction_reads1, junction_reads2 = {},{},{}
        cont = 0
        with open(readCounts_file) as f:
            header = next(f).rstrip()
            logger.info("Loading readCounts...")
            for line in f:
                # cont += 1
                # print(str(cont))
                # if (cont == 1000):
                #     break
                tokens = line.rstrip().split("\t")
                junction_id = tokens[0]
                reads = list(map(int,tokens[8:]))
                # Save this reads in the dictionary
                if(junction_id not in junction_reads):
                    junction_reads[junction_id] = reads
                else:
                    raise Exception("Junction id "+junction_id+" repeated")

        #Load the exonizations
        logger.info("Loading exonizations...")
        exonizations = pd.read_table(exonization_file, delimiter="\t")
        #Get the lists of junctions associated to each exonization
        Junction_id3 = exonizations["Junction_id3"].tolist()
        Junction_id4 = exonizations["Junction_id4"].tolist()
        final_counts_list = []
        logger.info("Processing files...")
        cont = 0
        #By line, check all the possible reads associated to each junction.
        for junction_list1,junction_list2 in zip(Junction_id3,Junction_id4):
            # cont += 1
            # print(str(cont))
            # if (cont == 3):
            #     break
            # By sample, sum all the reads by junctions
            junction_list1_aux = junction_list1.split(",")
            totalCounts1 = None
            for junction in junction_list1_aux:
                #Get the associated reads
                # If the junctions is not in the junction reads file, assign a 0
                if(junction in junction_reads):
                    # print("Encontrada junction!!!")
                    readCounts1 = junction_reads[junction]
                else:
                    readCounts1 = list(np.repeat(0,len(tokens[8:])))
                if(totalCounts1 is None):
                    totalCounts1 = readCounts1
                else:
                    totalCounts1 = [x + y for x, y in zip(totalCounts1, readCounts1)]
            # totalCounts1_list.append(totalCounts1)

            junction_list2_aux = junction_list2.split(",")
            totalCounts2 = None
            for junction in junction_list2_aux:
                #Get the associated reads
                # If the junctions is not in the junction reads file, assign a 0
                if(junction in junction_reads):
                    # print("Encontrada junction!!!")
                    readCounts2 = junction_reads[junction]
                else:
                    readCounts2 = list(np.repeat(0,len(tokens[8:])))
                if (totalCounts2 is None):
                    totalCounts2 = readCounts2
                else:
                    totalCounts2 = [x + y for x, y in zip(totalCounts2, readCounts2)]
            # totalCounts2_list.append(totalCounts2)

            #From both lists, take the minimum value
            final_counts = [str(min(x,y)) for x, y in zip(totalCounts1, totalCounts2)]
            aux = "\t".join(final_counts)
            final_counts_list.append(aux)

        #Paste the readCounts
        logger.info("Creating output file...")
        exonizations["readCounts"] = final_counts_list
        columns = list(exonizations.columns.values)
        columns[-1] = "\t".join(header.split("\t")[8:])
        exonizations.columns = columns
        #Save the file
        exonizations.to_csv(output_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar=' ')
        # #Save the file without using to_csv (it gave problems regarding the removal of quotations)
        # f = open(output_file, 'w')
        # f.write("\t".join(list(exonizations.columns.values)))
        # for index, row in exonizations.iterrows():
        #     f.write(row)

        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
