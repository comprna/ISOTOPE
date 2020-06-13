"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

extract_neoskipping_junctions: identify the junctions that could generate an neoskipping exon

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

description = \
    "Description:\n\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-i", "--input", required=True,
                    help="Input file")
parser.add_argument("-g", "--gtf", required=True,
                    help="Gtf file")
parser.add_argument("-t", "--threshold", required=False, type=float, default=0.0,
                    help="Number minimum of junction reads")
parser.add_argument("-o", "--output", required=True,
                    help="Output file")

def extract_neoskipping_junctions(input_path, gtf_path, threshold, output_path):

    try:
        logger.info("Starting execution")

        #Load the gtf file : Gene -> associated_exons
        gene_exons, gene_strand, gene_chr = {},{},{}
        cont = 0
        with open(gtf_path) as f:
            logger.info("Loading gtf exons...")
            for line in f:
                cont += 1
                # print(str(cont))
                # if (cont == 1000):
                #     break
                tokens = line.rstrip().split("\t")
                chr = tokens[0]
                start = int(tokens[3])
                end = int(tokens[4])
                strand = tokens[6]
                gene = tokens[8].split(";")[0].split("\"")[1]
                #Save the strand
                if(gene not in gene_strand):
                    gene_strand[gene] = strand
                else:
                    if(gene_strand[gene]!=strand):
                        raise Exception("Different strands associated to the same gene")
                # Save the chr
                if (gene not in gene_chr):
                    gene_chr[gene] = chr
                else:
                    if (gene_chr[gene] != chr):
                        raise Exception("Different chr associated to the same gene")
                #Save the exons
                if(gene not in gene_exons):
                    gene_exons[gene] = [[start, end]]
                else:
                    if([start, end] not in gene_exons[gene]):
                        gene_exons[gene].append([start, end])

        # Sort the list of exons
        # Depending on the strand, we have to sort it in a different way
        for aux_gene,exons_list in gene_exons.items():
            if (gene_strand[aux_gene] == "+"):
                exons_list_sorted = sorted(exons_list, key=lambda x: (x[1], x[0]))
            else:
                exons_list_sorted = sorted(exons_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_exons[aux_gene] = exons_list_sorted

        # Load the input file. Save all the junctions (type 2) to each associated gene: Gene -> junctions_type==2
        gene_junctions_1, gene_junctions_2, junction_type_dict = {}, {}, {}
        junction_reads = {}
        cont = 0
        with open(input_path) as f:
            header = next(f).rstrip().split("\t")[8:]
            logger.info("Loading junctions...")
            for line in f:
                cont += 1
                # print(str(cont))
                # if (cont == 1000):
                #     break
                tokens = line.rstrip().split("\t")
                junction_type = tokens[7]
                junction_id = tokens[0]
                start = int(junction_id.split(";")[1])
                end = int(junction_id.split(";")[2])
                gene_list = tokens[6].split(",")
                #Save also the reads associated
                reads = list(map(float,tokens[8:]))
                for x in gene_list:
                    gene_junction_id = x + ";" + tokens[0]
                    # Save the junctions only if the gene is in the gtf
                    if(junction_type=="1"):
                        if(x in gene_exons):
                            if(x not in gene_junctions_1):
                                gene_junctions_1[x] = [[start, end]]
                            else:
                                if([start, end] not in gene_junctions_1[x]):
                                    gene_junctions_1[x].append([start, end])
                            # Save this reads in the dictionary
                            if (gene_junction_id not in junction_reads):
                                junction_reads[gene_junction_id] = reads
                            else:
                                raise Exception("Junction id " + gene_junction_id + " repeated")
                        # Save the type of the junction
                        new_id = x+"_"+str(start)+"_"+str(end)
                        junction_type_dict[new_id] = junction_type
                    elif(junction_type=="2"):
                        if(x in gene_exons):
                            if(x not in gene_junctions_2):
                                gene_junctions_2[x] = [[start, end]]
                            else:
                                if([start, end] not in gene_junctions_2[x]):
                                    gene_junctions_2[x].append([start, end])
                            # Save this reads in the dictionary
                            if (gene_junction_id not in junction_reads):
                                junction_reads[gene_junction_id] = reads
                            else:
                                raise Exception("Junction id " + gene_junction_id + " repeated")
                        # Save the type of the junction
                        new_id = x+"_"+str(start)+"_"+str(end)
                        junction_type_dict[new_id] = junction_type
                    else:
                        pass

        # Sort the list of exons
        # Depending on the strand, we have to sort it in a different way
        for aux_gene, junction_list in gene_junctions_1.items():
            if (gene_strand[aux_gene] == "+"):
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]))
            else:
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_junctions_1[aux_gene] = junction_list_sorted

        for aux_gene, junction_list in gene_junctions_2.items():
            if (gene_strand[aux_gene] == "+"):
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]))
            else:
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_junctions_2[aux_gene] = junction_list_sorted

        # Go over the junctions lists and the exons and check if there is neoskipping
        logger.info("Processing files...")
        junctions_inside = {}
        #By gene, get all junctions==2 associated
        for gene,junction_list_2 in gene_junctions_2.items():
            if(gene in gene_exons):
                if(gene in gene_junctions_1):
                    junction_list_1 = gene_junctions_1[gene]
                    # if(strand=="+"):
                    #By junction, get all the junctions==1 associated to the same gene
                    for cont, junction_2 in enumerate(junction_list_2):
                        # Go over the junction list
                        for cont2,junction_1 in enumerate(junction_list_1):
                            # Save all the junctions that are falling inside of the neoskipping
                            if(int(junction_2[0])<=int(junction_1[0]) and int(junction_2[1])+1>=int(junction_1[1])):
                                junction_id_1_aux = gene + ";" + gene_chr[gene] + ";" + str(junction_1[0]) + ";" + str(junction_1[1]) + ";" + \
                                              gene_strand[gene]
                                junction_id_2_aux = gene + ";" + gene_chr[gene] + ";" + str(junction_2[0]) + ";" + str(junction_2[1]) + ";" + \
                                              gene_strand[gene]
                                if(junction_id_2_aux not in junctions_inside):
                                    junctions_inside[junction_id_2_aux] = [junction_id_1_aux]
                                else:
                                    junctions_inside[junction_id_2_aux].append(junction_id_1_aux)
                else:
                    pass
            else:
                # logger.info("Gene "+gene+" not in gtf")
                pass

        # Output the values obtaining the reads associated. If the reads is over the threshold, output the sample
        logger.info("Output values...")
        output_file = open(output_path, 'w')
        output_file.write("Sample_id\tGene_id\tNeoskipping_junction\tSkipped_junctions\tNeoskipping_ReadCounts\tmean_skipped_ReadCounts\n")
        cont = 0
        for neoskipping_junction,skipped_junctions in junctions_inside.items():
            cont += 1
            # print(str(cont))
            # Get the gene associated
            gene = neoskipping_junction.split(";")[0]
            junction_aux = neoskipping_junction.split(";")[1:]
            # Get the reads associated to the junction
            if(neoskipping_junction in junction_reads):
                neoskipped_reads = junction_reads[neoskipping_junction]
                # Look for reads above the threshold
                for i in range(0,len(neoskipped_reads)):
                    reads1 = float(neoskipped_reads[i])
                    if(reads1>=threshold):
                        total_reads, total_junctions = [], []
                        # Get also the reads associated to the skipped junctions for this sample
                        for x in skipped_junctions:
                            junction_aux2 = ";".join(x.split(";")[1:])
                            total_junctions.append(junction_aux2)
                            if (x in junction_reads):
                                reads2 = float((junction_reads[x])[i])
                                total_reads.append(reads2)
                            else:
                                logger.info("Junction " + x + " not in junction_reads")
                        # Get the mean
                        mean_reads = np.mean(total_reads)
                        #Output the values
                        sample_id = header[i]
                        output_file.write(sample_id+"\t"+gene+"\t"+";".join(junction_aux)+"\t"+
                                          ";".join(total_junctions)+"\t"+str(reads1)+"\t"+str(mean_reads)+"\n")

            else:
                logger.info("Junction "+neoskipping_junction+" not in junction_reads")

        output_file.close()

        # Get the fold-times there are reads in the neoskipping junction comapred with the others
        file = pd.read_table(output_path, delimiter="\t" )
        file["Fold"] = file.apply(lambda x: float(x[4])/float(x[5]) if float(x[5])!=0 else "Inf",axis=1)
        file.to_csv(output_path, sep="\t", index=False)

        logger.info("Created "+output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)