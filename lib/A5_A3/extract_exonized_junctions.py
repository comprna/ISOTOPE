"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

extract_exonized_junctions: identify the junctions that could generate an exonization

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


description = \
    "Description:\n\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-i", "--input", required=True,
                    help="Input file")
parser.add_argument("-g", "--gtf", required=True,
                    help="Gtf file")
parser.add_argument("-l", "--length", required=False, type=int, default=0,
                    help="Max length for the exonizations")
parser.add_argument("-o", "--output", required=True,
                    help="Output file")

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


def extract_exonized_junctions(input_path, gtf_path, max_length, output_path):

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

        # Load the input file. Save all the junctions (except 5) to each associated gene
        gene_junctions_cannonical, gene_junctions_alternative, junction_type_dict = {}, {}, {}
        cont = 0
        with open(input_path) as f:
            logger.info("Loading junctions...")
            next(f)
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

        # Go over the junctions lists and the exons and check if there are A5_A3
        logger.info("Processing files...")
        exon_junction3, exon_junction4, exon_gene, exon_lengths = {}, {}, {}, {}
        flag3, flag4 = False, False
        output_file = open(output_path, "w")
        output_file.write(
            "Gene\tCanonical_Junction_id\tAlt_Junction_id\tCanonical_Exon\tAlt_Exon_id\tstrand\tOffset\tNew_Exon_length\n")
        for gene, junction_list in gene_junctions_cannonical.items():
            # print("Gene "+gene)
            if (gene in gene_exons):
                strand = gene_strand[gene]
                exons_list = gene_exons[gene]
                chr = gene_chr[gene]
                if (strand == "+"):
                    # print("Strand +")
                    for cont, junction in enumerate(junction_list):
                        # print("junction " + ",".join(map(str,junction)))
                        # Get all the junctions falling in the same start position
                        if (gene in gene_junctions_alternative):
                            junction_list_alt = gene_junctions_alternative[gene]
                            start_associated_junctions = get_junctions(junction, junction_list_alt, True)
                            # print("start")
                            if (len(start_associated_junctions) > 1):
                                # Get the exons where the junction end is falling (for obtaining the total length of the aberrant exon)
                                exon = get_associated_exon(junction[1], exons_list)
                                if (exon != ""):
                                    for x in start_associated_junctions:
                                        id_can = chr + ";" + str(junction[0]) + ";" + str(junction[1]) + ";" + strand
                                        id_alt = chr + ";" + str(x[0]) + ";" + str(x[1]) + ";" + strand
                                        exon_can = chr + ";" + str(exon[0]) + ";" + str(exon[1]) + ";" + strand
                                        exon_alt = chr + ";" + str(x[1]) + ";" + str(exon[1]) + ";" + strand
                                        # If the junction is falling upstream
                                        if (x[1] < exon[0]):
                                            new_length = (exon[0] - x[1]) + (exon[1] - exon[0])
                                            offset = exon[0] - x[1]
                                            # If the resultant exon is shorter than n, then output the aberrant A5
                                            if (new_length <= max_length):
                                                output_file.write(
                                                    gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                    "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")
                                        # If the junction is falling inside the exon
                                        elif (x[1] < exon[1]):
                                            new_length = (exon[0] - x[1]) + (exon[1] - exon[0])
                                            offset = exon[0] - x[1]
                                            output_file.write(
                                                gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")

                            # Get all the junctions falling in the same end position
                            end_associated_junctions = get_junctions(junction, junction_list_alt, False)
                            # print("end")
                            if (len(end_associated_junctions) > 1):
                                # Get the exons where the junction end is falling (for obtaining the total length of the aberrant exon)
                                exon = get_associated_exon(junction[0], exons_list)
                                if (exon != ""):
                                    for x in end_associated_junctions:
                                        id_can = chr + ";" + str(junction[0]) + ";" + str(junction[1]) + ";" + strand
                                        id_alt = chr + ";" + str(x[0]) + ";" + str(x[1]) + ";" + strand
                                        exon_can = chr + ";" + str(exon[0]) + ";" + str(exon[1]) + ";" + strand
                                        exon_alt = chr + ";" + str(exon[0]) + ";" + str(x[0]+1) + ";" + strand
                                        # If the junctions is falling upstream
                                        if (exon[1] < x[0]):
                                            new_length = (x[0] - exon[1]) + (exon[1] - exon[0])
                                            offset = x[0] - exon[1]
                                            # If the resultant exon is shorter than n, then output the aberrant A5
                                            if (new_length <= max_length):
                                                output_file.write(
                                                    gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                    "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")
                                        # If the junction is falling inside the exon
                                        elif (exon[0] < x[0]):
                                            new_length = (x[0] - exon[1]) + (exon[1] - exon[0])
                                            offset = x[0] - exon[1]
                                            output_file.write(
                                                gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")
                        else:
                            # print("No alt junctions for this gene "+gene)
                            pass

                else:  # strand=="-"
                    # print("Strand -")
                    for cont, junction in enumerate(junction_list):
                        # print("junction " + ",".join(map(str,junction)))
                        # Get all the junctions falling in the same start position
                        if (gene in gene_junctions_alternative):
                            junction_list_alt = gene_junctions_alternative[gene]
                            start_associated_junctions = get_junctions(junction, junction_list_alt, False)
                            # print("start")
                            if (len(start_associated_junctions) > 1):
                                # Get the exons where the junction end is falling (for obtaining the total length of the aberrant exon)
                                exon = get_associated_exon(junction[0], exons_list)
                                if (exon != ""):
                                    for x in start_associated_junctions:
                                        id_can = chr + ";" + str(junction[0]) + ";" + str(junction[1]) + ";" + strand
                                        id_alt = chr + ";" + str(x[0]) + ";" + str(x[1]) + ";" + strand
                                        exon_can = chr + ";" + str(exon[0]) + ";" + str(exon[1]) + ";" + strand
                                        # exon_alt = chr + ";" + str(x[1]) + ";" + str(exon[1]) + ";" + strand
                                        exon_alt = chr + ";" + str(exon[0]) + ";" + str(x[0] + 1) + ";" + strand
                                        # If the junction is falling upstream
                                        if (x[0] > exon[1]):
                                            offset = x[0] - exon[1]
                                            new_length = offset + (exon[1] - exon[0])
                                            # If the resultant exon is shorter than n, then output the aberrant A5
                                            if (new_length <= max_length):
                                                output_file.write(
                                                    gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                    "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")
                                        # If the junction is falling inside the exon
                                        elif (x[0] > exon[0]):
                                            offset = x[0] - exon[1]
                                            new_length = offset + (exon[1] - exon[0])
                                            output_file.write(
                                                gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" + exon_alt +
                                                "\t" + strand + "\t" + str(offset) + "\t" + str(new_length) + "\n")

                            # Get all the junctions falling in the same end position
                            end_associated_junctions = get_junctions(junction, junction_list_alt, True)
                            # print("end")
                            if (len(end_associated_junctions) > 1):
                                # Get the exons where the junction end is falling (for obtaining the total length of the aberrant exon)
                                exon = get_associated_exon(junction[1], exons_list)
                                if (exon != ""):
                                    for x in end_associated_junctions:
                                        id_can = chr + ";" + str(junction[0]) + ";" + str(junction[1]) + ";" + strand
                                        id_alt = chr + ";" + str(x[0]) + ";" + str(x[1]) + ";" + strand
                                        exon_can = chr + ";" + str(exon[0]) + ";" + str(exon[1]) + ";" + strand
                                        exon_alt = chr + ";" + str(x[1]) + ";" + str(exon[1]) + ";" + strand
                                        # If the junctions is falling upstream
                                        if (exon[0] > x[1]):
                                            offset = exon[0] - x[1]
                                            new_length = offset + (exon[1] - exon[0])
                                            # If the resultant exon is shorter than n, then output the aberrant A5
                                            if (new_length <= max_length):
                                                output_file.write(
                                                    gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" +
                                                    exon_alt + "\t" + strand + "\t" + str(offset) + "\t" + str(
                                                        new_length) + "\n")
                                        # If the junction is falling inside the exon
                                        elif (exon[1] > x[1]):
                                            offset = exon[0] - x[1]
                                            new_length = offset + (exon[1] - exon[0])
                                            output_file.write(
                                                gene + "\t" + id_can + "\t" + id_alt + "\t" + exon_can + "\t" +
                                                exon_alt + "\t" + strand + "\t" + str(offset) + "\t" + str(
                                                    new_length) + "\n")
                        else:
                            # print("No alt junctions for this gene "+gene)
                            pass

        output_file.close()

        # Get the motifs of the splice sites
        # Load the file we just created
        logger.info("Obtaining new motifs...")
        logger.info("loading file...")
        file = pd.read_table(output_path, delimiter="\t")
        # Transform the junction cordinates to bed file
        chr = file["Canonical_Junction_id"].apply(lambda x: x.split(";")[0])
        # Depending on the strand and the junctions we will return the splice site of the acceptor site or the donor
        logger.info("get_info_splice_sites...")
        aux = file.apply(get_info_splice_sites, axis=1)
        # logger.info("info...")
        info = aux.apply(lambda x: x[0])
        # logger.info("strand...")
        strand = aux.apply(lambda x: x[1])
        start = aux.apply(lambda x: x[2])
        end = aux.apply(lambda x: x[3])
        # logger.info("id...")
        # id = file.apply(lambda x: x["Gene"]+"|"+x["Canonical_Junction_id"]+"|"+x["Alt_Junction_id"]+"|"+info,axis=1)

        logger.info("creating bed...")
        bed = [("chr", chr), ("start", start), ("end", end), ("strand", strand)]
        bed_file = pd.DataFrame.from_items(bed)
        # bed_file['id'] = id
        bed_file['score'] = 0
        # Save this variables as bed file
        path1 = "/".join(output_path.split("/")[:-1])
        logger.info("Saving bed files...")
        bed_file.to_csv(path1 + "/aux_bed.bed", sep="\t", index=False, header=False)
        # Run getFasta from MoSEA (it needs python2 and BEDTools)
        logger.info("Running Mosea...")
        command1 = "module load Python/2.7.11; module load BEDTools; python2 /genomics/users/juanluis/Software/MoSEA-master/mosea.py getfasta --bedfile " + \
                   path1 + "/aux_bed.bed --genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa " \
                           "--output " + path1 + "/aux_bed.fasta"
        os.system(command1)

        # Load the motifs and associate it to the exons. Make it a table
        motif, position = [], []
        with open(path1 + "/aux_bed.fasta") as f:
            for line in f:
                if (line[0] == ">"):
                    pass
                else:
                    motif.append(line.rstrip())

        # Incorporate this info to the final df
        file["motif"] = motif
        file["splice_site_type"] = info
        file.to_csv(output_path, sep="\t", index=False)

        logger.info("Saved " + output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
