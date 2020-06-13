"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

compare_reads_random_junctions.py: for applying some filtering on the list of A5_A3 junctions, we are gonna
compare the readcounts for each junction against other new junctions associated to the same gene

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import random
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.sandbox.stats.multicomp import multipletests


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


def get_junctions(junction, junction_list, start):
    '''
    Get all the junctions that are overlapping the pos
    '''
    output_junction_list = []
    if (start):
        for x in junction_list:
            if (junction[0] == x[0] and (junction[0] != x[0] or junction[1] != x[1])):
                output_junction_list.append(x)
    else:
        for x in junction_list:
            if (junction[1] == x[1] and (junction[0] != x[0] or junction[1] != x[1])):
                output_junction_list.append(x)
    return output_junction_list


def get_associated_exon(pos, exons_list):
    '''
    Return the exons where the pos is falling
    '''
    for x in exons_list:
        if (x[0] <= pos <= x[1]):
            return x
    return ""


def get_info_splice_sites(line):
    '''
    Get the info of the coordinates of the new produced motif site
    '''
    can_junc = line["Canonical_Junction_id"].split(";")
    alt_junc = line["Alt_Junction_id"].split(";")
    if (can_junc[3] != alt_junc[3]):
        raise Exception("Different strand associated to the same exonization")
    else:
        strand = can_junc[3]
    if (strand == "+"):
        if (can_junc[2] == alt_junc[2]):
            # New donor site
            start = int(alt_junc[1]) + 1
            end = int(alt_junc[1]) + 3
            info = "New_donor"
        elif (can_junc[1] == alt_junc[1]):
            # New acceptor site
            start = int(alt_junc[2]) - 3
            end = int(alt_junc[2]) - 1
            info = "New_acceptor"
    else:  # strand=="-"
        if (can_junc[2] == alt_junc[2]):
            # New acceptor site
            start = int(alt_junc[1])  +1
            end = int(alt_junc[1]) + 3
            info = "New_acceptor"
        elif (can_junc[1] == alt_junc[1]):
            # New donor site
            start = int(alt_junc[2]) - 3
            end = int(alt_junc[2]) - 1
            info = "New_donor"
    return info, strand, start, end


def compare_reads_random_junctions(input_path, readCounts_path, gtf_path, output_path):

    try:
        logger.info("Starting execution")

        # Load the gtf file
        gene_exons, gene_strand, gene_chr = {}, {}, {}
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
                # Save the strand
                if (gene not in gene_strand):
                    gene_strand[gene] = strand
                else:
                    if (gene_strand[gene] != strand):
                        raise Exception("Different strands associated to the same gene")
                # Save the chr
                if (gene not in gene_chr):
                    gene_chr[gene] = chr
                else:
                    if (gene_chr[gene] != chr):
                        raise Exception("Different chr associated to the same gene")
                # Save the exons
                if (gene not in gene_exons):
                    gene_exons[gene] = [[start, end]]
                else:
                    if ([start, end] not in gene_exons[gene]):
                        gene_exons[gene].append([start, end])

        # Sort the list of exons
        # Depending on the strand, we have to sort it in a different way
        for aux_gene, exons_list in gene_exons.items():
            if (gene_strand[aux_gene] == "+"):
                exons_list_sorted = sorted(exons_list, key=lambda x: (x[1], x[0]))
            else:
                exons_list_sorted = sorted(exons_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_exons[aux_gene] = exons_list_sorted

        # Load the input file. Save all the junctions (except 5) to each associated gene and the associated reads
        gene_junctions_cannonical, gene_junctions_alternative, junction_type_dict, junction_reads = {}, {}, {}, {}
        cont = 0
        with open(readCounts_path) as f:
            logger.info("Loading readCounts...")
            header1 = next(f).rstrip().split("\t")
            header1 = header1[8:]
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
                reads_list = tokens[8:]
                for x in gene_list:
                    # Save the junctions only if the gene is in the gtf
                    if (x in gene_exons):
                        if (junction_type == "1"):
                            if (x not in gene_junctions_cannonical):
                                gene_junctions_cannonical[x] = [[start, end]]
                            else:
                                if ([start, end] not in gene_junctions_cannonical[x]):
                                    gene_junctions_cannonical[x].append([start, end])
                        elif (junction_type == "3" or junction_type == "4"):
                            if (x not in gene_junctions_alternative):
                                gene_junctions_alternative[x] = [[start, end]]
                            else:
                                if ([start, end] not in gene_junctions_alternative[x]):
                                    gene_junctions_alternative[x].append([start, end])
                        # Save the type of the junction
                        new_id = x + "_" + str(start) + "_" + str(end)
                        junction_type_dict[new_id] = junction_type
                        junction_reads[new_id] = reads_list
                    else:
                        # logger.info("Gene " + gene + " not in gtf")
                        pass

        # Sort the list of exons
        # Depending on the strand, we have to sort it in a different way
        for aux_gene, junction_list in gene_junctions_cannonical.items():
            if (gene_strand[aux_gene] == "+"):
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]))
            else:
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_junctions_cannonical[aux_gene] = junction_list_sorted

        for aux_gene, junction_list in gene_junctions_alternative.items():
            if (gene_strand[aux_gene] == "+"):
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]))
            else:
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_junctions_alternative[aux_gene] = junction_list_sorted

        # Go over the junctions lists, consulting all the other de novo junctions associated to the same gene
        pvalue_list, total_junctions_list = [], []
        cont = 0
        with open(input_path) as f:
            logger.info("Processing input file...")
            header2 = next(f).rstrip().split("\t")
            alt_junction_pos = header2.index("Alt_Junction_id")
            sample_id_pos = header2.index("Sample_id")
            gene_pos = header2.index("Gene")
            readCounts_pos = header2.index("ReadCounts")
            for line in f:
                cont += 1
                # logger.info(str(cont))
                tokens = line.rstrip().split("\t")
                alt_junction = tokens[alt_junction_pos]
                sample_id = tokens[sample_id_pos].strip()
                #Get the relative position of the sample in the header, for obtaining the associated reads
                pos_sample_id = header1.index(sample_id)
                gene = tokens[gene_pos]
                reads_junction = float(tokens[readCounts_pos].strip())
                #Recover all the_junctions associated to the gene
                junction_list_alt = gene_junctions_alternative[gene]
                total_junctions = len(junction_list_alt)
                total_junctions_list.append(total_junctions)
                # logger.info("Number of junctions "+str(total_junctions))
                #If there is at least 10 alternative associated junctions, take 50% of this for comparing against our number of read counts
                if(total_junctions>=10):
                    n_junctions = int(total_junctions/2)
                    if(n_junctions<10):
                        n_junctions = 10
                    cont2 = 0
                    list_selected = []
                    while cont2 < n_junctions:
                        # Generate a random start
                        random_number = random.randint(0, n_junctions-1)
                        # If this junction has not been chosen before, extract it
                        if(random_number not in list_selected):
                            list_selected.append(random_number)
                            cont2 += 1
                        else:
                            continue
                    #One we have all the selected junctions, recover the readcounts for each of them
                    reads_random = []
                    for x in list_selected:
                        id = junction_list_alt[x]
                        formatted_id = gene + "_" + str(id[0]) + "_" + str(id[1])
                        reads = junction_reads[formatted_id]
                        #Get the reads associated to our sample
                        reads_random.append(float(reads[pos_sample_id].strip()))
                    #Compare the reads from the original junction against the reads_random
                    ecdf = ECDF(reads_random)
                    pvalue = (1.0 - ecdf(reads_junction))
                    pvalue_list.append(pvalue)
                else:   #less than 10 junctions
                    pvalue_list.append(1)

        input_file = pd.read_table(input_path, delimiter="\t")
        input_file["pvalue_comparison_junctions"] = pvalue_list
        #Get the FDR
        _, pvals_corrected, _, _ = multipletests(pvalue_list, method='fdr_bh', alpha=0.05)
        input_file["FDR_comparison_junctions"] = pvals_corrected
        input_file["total_junctions_in_same_gene"] = total_junctions_list
        #Get only the significant cases according to the FDR
        input_file_f = input_file.loc[input_file["FDR_comparison_junctions"] <= 0.05]

        # Save the file
        input_file_f.to_csv(output_path, sep="\t", index=False)

        logger.info("Saved " + output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
