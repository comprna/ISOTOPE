"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

generate_random_intronic_positions: generate a number of random position by exonization
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import matplotlib as mtplot
mtplot.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)
import random
import copy

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

def overlap_any_region(coordinates,list_of_exons):
    '''
    Go over the list_of_exons. If coordinates overlap on any exon, return True.
    If dont, return False
    '''
    for x in list_of_exons:
        if(int(x[0])<=coordinates[0]<=int(x[1]) or int(x[0])<=coordinates[1]<=int(x[1])):
            return True
    return False

def get_hg_chromosome_id(x):
    '''
    Get the corresponding numeric chromosome
    '''
    if(x[3:]=="X"):
        return 23
    elif(x[3:]=="Y"):
        return 24
    else:
        return int(x[3:])


def generate_random_intronic_positions(input_path, gtf_path, n, output_path, output_path2):

    try:
        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # n = int(sys.argv[3])
        # output_path = sys.argv[4]
        # output_path2 = sys.argv[5]

        # input_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # n = 100
        # output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/random_exonizations.gtf"
        # output_path2 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/random_exonizations.bed"

        # all_sample_ids = ["S00035T","S00213","S00825","S00829","S00831","S00832","S00838","S01242","S01297T","S01366",
        #                   "S01512","S01524","S01542","S01578","S01861T","S01864","S01873T","S02065","S02093","S02120",
        #                   "S02139","S02163","S02194","S02209","S02234","S02241","S02242","S02243","S02244","S02246",
        #                   "S02248","S02249","S02256","S02284","S02285","S02286","S02287","S02288","S02289","S02290",
        #                   "S02291","S02293","S02294","S02295","S02296","S02297","S02298","S02299","S02322","S02328",
        #                   "S02360","S02375T","S02376T","S02378T","S02382T","S02397","S00022","S00050","S00472","S00501",
        #                   "S00827","S00830","S00836","S00837","S01453","S01494","S01556","S01563","S01698","S01728",
        #                   "S01792","13368","31052","31056","31064","31076","31092","34398","34413","34417","34421",
        #                   "34426","34427","34430","85203","85205","85208","85210","85223","85258","85260","85267","85268",
        #                   "85270","85272","85276","85281","8687","8711","8735","8742","88421","SM09-001T","SM09-002T",
        #                   "SM09-003T","SM09-004T","SM09-005T","SM09-006T","SM09-007T","SM09-008T","SM09-010T",
        #                   "SM09-011T1","SM09-011T2","SM09-012T","SM09-013T","SM09-014T","SM09-015T","SM09-016T",
        #                   "SM09-017T","SM09-018T","SM09-019T","SM09-020T"]

        #Read the GTF data. Save the info of the exons by gene
        logger.info("Loading GTF...")
        gene_exons = {}
        with open(gtf_path) as f:
            for line in f:
                tokens = line.rstrip().split("\t")
                start = int(tokens[3])
                end = int(tokens[4])
                strand = tokens[6]
                info = tokens[8]
                gene = info.split(";")[0].split("\"")[1]
                if(gene not in gene_exons):
                    gene_exons[gene] = [(start,end)]
                else:
                    gene_exons[gene].append((start,end))

        # For each gene, obtain the genomic boundaries. First, sort all genes by start coordinates
        logger.info("Sorting coordinates...")
        gene_start_end = {}
        gene_exons_copy = copy.deepcopy(gene_exons)
        for key,values in gene_exons_copy.items():
            #If the gene only has 1 exon, delete it
            if(len(values)==1):
                del gene_exons[key]
            else:
                #Sort the values
                gene_exons[key] = sorted(values, key=lambda k: [k[0], k[1]])
                #Get the end of the first exon and the start of the last
                end = gene_exons[key][0][1]
                start = gene_exons[key][len(values)-1][0]
                if (key not in gene_start_end):
                    gene_start_end[key] = (end, start)
                else:
                    logger.info("Repeated gene "+key+"!!!")

        #Read the input data
        exonizations_dict = {}
        cont = 0
        output_file = open(output_path, "w")
        output_file2 = open(output_path2, "w")
        with open(input_path) as f:
            logger.info("Processing input data...")
            next(f)
            for line in f:
                # if (cont == 1000):
                #     break
                tokens = line.rstrip().split("\t")
                id = tokens[1]
                gene = tokens[3]
                #If the exonization is not repeated and the gene is presented in the gtf, generate randomizations
                # Generate a random number of positions.
                if(id not in exonizations_dict and gene in gene_start_end):
                    print(str(cont))
                    cont += 1
                    chr = id.split(";")[0]
                    start = id.split(";")[1]
                    end = id.split(";")[2]
                    strand = id.split(";")[3]
                    length = int(end) - int(start)
                    #Check if any of the positions are overlapping nor the exons in the GTF neither the exonization
                    coordinates = gene_start_end[gene]
                    cont2 = 0
                    while cont2<n:
                        #Generate a random start
                        random_start = random.randint(int(coordinates[0]), int(coordinates[1]))
                        random_end = random_start + length
                        #Check if this coordinates overlap any exonic region
                        if(overlap_any_region((random_start,random_end),gene_exons[gene])):
                            continue
                        cont2 += 1
                        #Output the gtf
                        new_name = "\"Exonization_" + str(cont) + "_Random_" + str(cont2) + "\""
                        output_file.write(chr+"\trandom_exonization\texon\t"+str(random_start)+"\t"+str(random_end)+"\t.\t"+strand+
                                          "\t.\tgene_id \""+gene+"\"; transcript_id "+new_name+"; exon_number \""+
                                          str(cont2)+"\"; gene_name \""+gene+"\"; transcript_name "+new_name+"; exon_id "+
                                          new_name+";\n")
                        # Output the bed
                        output_file2.write(chr+"\t"+str(random_start)+"\t"+str(random_end)+"\t"+new_name+"\t"+strand+"\t0\n")
                    #Generate also the coordinates for the exonization
                    # Output the gtf
                    new_name = "\"Exonization_" + str(cont) + "\""
                    output_file.write(chr + "\texonization\texon\t" + str(start) + "\t" + str(end) + "\t.\t" + strand +
                                      "\t.\tgene_id " + new_name + "; transcript_id " + new_name + "; exon_number " +
                                      str(cont) + "; gene_name " + new_name + "; transcript_name " + new_name + "; exon_id " +
                                      new_name + ";\n")
                    # Output the bed
                    output_file2.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + new_name + "\t" + strand +"\t0\n")
                    # Save this exonizations, for not generating repeated exonizations
                    exonizations_dict[id] = new_name
                else:
                    pass

        output_file.close()
        output_file2.close()

        #Load the file again and sort it by chromosome, start and end
        file = pd.read_table(output_path2, delimiter="\t")
        file.columns = ["chr","start","end","id","strand","score"]
        file["chr_num"] = file["chr"].apply(get_hg_chromosome_id)
        # Sort by sample id
        file_sorted = file.sort_values(["chr_num","start","end"])
        del file_sorted["chr_num"]
        file_sorted.to_csv(output_path2, sep="\t", index=False, header=False)

        logger.info("Saved "+output_path)
        logger.info("Saved "+output_path2)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
