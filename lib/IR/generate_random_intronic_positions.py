"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

generate_random_intronic_positions: generate a number of random position by intron retention
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
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

def get_coordinates(values):
    '''
    Return the start more upstream and the end more downstream
    '''
    start, end = 9999999999999999999, 0
    for tuple in values:
        if(tuple[1]<start):
            start=tuple[1]
        if (tuple[0] > end):
            end = tuple[0]
    return (start,end)


def generate_random_intronic_positions(input_path, gtf_path, n, output_path, output_path2):
# def generate_random_intronic_positions():

    try:
        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # n = int(sys.argv[3])
        # output_path = sys.argv[4]
        # output_path2 = sys.argv[5]
        #
        # input_path = "/projects_rg/test_ePydoor/IR_expressed_genes.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # n = 100
        # output_path = "/projects_rg/test_ePydoor/random_introns.gtf"
        # output_path2 = "/projects_rg/test_ePydoor/random_introns.bed"

        # Read the GTF data. Save the info of the exons by gene
        logger.info("Loading GTF...")
        gene_exons = {}
        with open(gtf_path) as f:
            for line in f:
                tokens = line.rstrip().split("\t")
                start = int(tokens[3])
                end = int(tokens[4])
                strand = tokens[6]
                info = tokens[8]
                # This line obtains the Ensembl gene
                gene = info.split(";")[0].split("\"")[1]
                # This line obtains the Gene_name
                # gene = info.split(";")[3].split("\"")[1]
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
                #Get the start more upstream and the end more downstream
                start,end = get_coordinates(values)
                #Get the end of the first exon and the start of the last
                # end = gene_exons[key][0][1]
                # start = gene_exons[key][len(values)-1][0]
                if (key not in gene_start_end):
                    # gene_start_end[key] = (end, start)
                    gene_start_end[key] = (start , end)
                else:
                    logger.info("Repeated gene "+key+"!!!")

        #Read the input data
        logger.info("Loading input data...")
        # Load the table using pandas
        input = pd.read_table(input_path, delimiter="\t")
        # Get the unique values of the Event_id,gene
        unique_values = input.drop_duplicates(['Event_id','Gene_id'])

        introns_dict = {}
        cont = 0
        output_file = open(output_path, "w")
        output_file2 = open(output_path2, "w")
        # with open(input_path) as f:
        #     logger.info("Processing input data...")
        #     header = next(f).rstrip().split("\t")
        #     Event_id_pos = header.index("Event_id")
        #     Gene_id_pos = header.index("Gene_id")
        #     for line in f:
        #         tokens = line.rstrip().split("\t")
        #         id = tokens[Event_id_pos]
        #         gene = tokens[Gene_id_pos]

        for i in range(0,len(unique_values.index)):
            logger.info("i: "+str(i))
            id = unique_values.Event_id.iloc[i]
            gene = unique_values.Gene_id.iloc[i]
            #If there are more than one gene, take the one with the greatest interval
            flag_found_gene = False
            if(len(gene.split(","))==1):
                if(gene in gene_start_end):
                    flag_found_gene = True
            else:
                coordinates_max = (0,0)
                for x in gene.split(","):
                    if(x in gene_start_end):
                        coordinates = gene_start_end[x]
                        if(int(coordinates_max[1]-coordinates_max[0])<int(coordinates[1]-coordinates[0])):
                            flag_found_gene = True
                            coordinates_max = coordinates
                            gene = x
            #If the intron is not repeated and the gene is presented in the gtf, generate randomizations
            # Generate a random number of positions.
            if(id not in introns_dict and flag_found_gene):
                cont += 1
                # print(str(cont))
                # print(str(id))
                # if(id=="chrY:9311665-9323642(-):kma_introns"):
                #     print("jobuoñbo")
                chr = id.split(":")[0]
                start = id.split(":")[1].split("(")[0].split("-")[0]
                end = id.split(":")[1].split("(")[0].split("-")[1]
                strand = id.split(":")[1][-2]
                length = int(end) - int(start)
                if(length<0):
                    raise Exception("Interval incorrect for id "+id)
                #Check if any of the positions are overlapping nor the exons in the GTF neither the exonization
                coordinates = gene_start_end[gene]
                cont2 = 0
                cont3 = 0
                flag_escape = False
                while cont2<n:
                    # If the program reach to a reasonable number of attempts, skip this line
                    if(cont3>100000):
                        flag_escape = True
                        break
                    cont3 += 1
                    #Generate a random start
                    random_start = random.randint(int(coordinates[0]), int(coordinates[1]))
                    random_end = random_start + length
                    #Check if this coordinates overlap any exonic region
                    if(overlap_any_region((random_start,random_end),gene_exons[gene])):
                        continue
                    cont2 += 1
                    #Output the gtf
                    new_name = "\"Intron_" + str(cont) + "_Random_" + str(cont2) + "\""
                    output_file.write(chr+"\trandom_intron\texon\t"+str(random_start)+"\t"+str(random_end)+"\t.\t"+strand+
                                      "\t.\tgene_id \""+gene+"\"; transcript_id "+new_name+"; exon_number \""+
                                      str(cont2)+"\"; gene_name \""+gene+"\"; transcript_name "+new_name+"; exon_id "+
                                      new_name+";\n")
                    # Output the bed
                    output_file2.write(chr+"\t"+str(random_start)+"\t"+str(random_end)+"\t"+new_name+"\t"+strand+"\t0\n")
                if(flag_escape):
                    logger.info("Not possible to obtain random intronic regions for "+gene)
                    continue
                #Generate also the coordinates for the intron
                # Output the gtf
                new_name = "\"Intron_" + str(cont) + "\""
                output_file.write(chr + "\tintron\texon\t" + str(start) + "\t" + str(end) + "\t.\t" + strand +
                                  "\t.\tgene_id " + new_name + "; transcript_id " + new_name + "; exon_number " +
                                  str(cont) + "; gene_name " + new_name + "; transcript_name " + new_name + "; exon_id " +
                                  new_name + ";\n")
                # Output the bed
                output_file2.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + new_name + "\t" + strand +"\t0\n")
                # Save this exonizations, for not generating repeated exonizations
                introns_dict[id] = new_name
            else:
                if(gene not in gene_start_end):
                    # logger.info("Gene "+gene+" not in gtf")
                    pass

        output_file.close()
        output_file2.close()

        #Load the file again and sort it by chromosome, start and end
        logger.info("Sort the file by chromosome...")
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

    except KeyboardInterrupt:
        logger.error('Interrupted program')
        logger.info(str(cont) + ": Number of attempts: " + str(cont3))
        sys.exit(1)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)

# if __name__ == '__main__':
#     generate_random_intronic_positions()