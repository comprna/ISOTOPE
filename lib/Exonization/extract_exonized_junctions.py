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

# description = \
#     "Description:\n\n"
#
# parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
#                         add_help=True)
# parser.add_argument("-i", "--input", required=True,
#                     help="Input file")
# parser.add_argument("-g", "--gtf", required=True,
#                     help="Gtf file")
# parser.add_argument("-l", "--length", required=False, type=int, default=0,
#                     help="Max length for the exonizations")
# parser.add_argument("-o", "--output", required=True,
#                     help="Output file")

def extract_exonized_junctions(input_path, gtf_path, max_length, output_path):

    try:
        logger.info("Starting execution")

        #Load the gtf file
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

        # Load the input file. Save all the junctions (type 3 or 4) to each associated gene
        # Check if any of the junction (type 3 or 4) overlapps with any exon
        gene_junctions, junction_type_dict = {}, {}
        cont = 0
        with open(input_path) as f:
            logger.info("Loading junctions...")
            for line in f:
                cont += 1
                # print(str(cont))
                # if (cont == 1000):
                #     break
                tokens = line.rstrip().split("\t")
                junction_type = tokens[7]
                if(junction_type=="3" or junction_type=="4"):
                    junction_id = tokens[0]
                    # chr = int(junction_id.split(";")[0])
                    start = int(junction_id.split(";")[1])
                    end = int(junction_id.split(";")[2])
                    # strand = junction_id.split(";")[3]
                    gene_list = tokens[6].split(",")
                    for x in gene_list:
                        # Save the junctions only if the gene is in the gtf
                        if(x in gene_exons):
                            if(x not in gene_junctions):
                                gene_junctions[x] = [[start, end]]
                            else:
                                if([start, end] not in gene_junctions[x]):
                                    gene_junctions[x].append([start, end])
                            # Save the type of the junction
                            new_id = x+"_"+str(start)+"_"+str(end)
                            junction_type_dict[new_id] = junction_type
                        else:
                            # logger.info("Gene " + gene + " not in gtf")
                            pass

        # Sort the list of exons
        # Depending on the strand, we have to sort it in a different way
        for aux_gene, junction_list in gene_junctions.items():
            if (gene_strand[aux_gene] == "+"):
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]))
            else:
                junction_list_sorted = sorted(junction_list, key=lambda x: (x[1], x[0]), reverse=True)
            gene_junctions[aux_gene] = junction_list_sorted

        # Go over the junctions lists and the exons and check if there is exonization
        #For each junction type 3, check all the exons. If there is a overlapping, check for junctions type 4
        logger.info("Processing files...")
        exon_junction3, exon_junction4, exon_gene, exon_lengths = {},{},{},{}
        flag3, flag4 = False, False
        for gene,junction_list in gene_junctions.items():
            if(gene in gene_exons):
                strand = gene_strand[gene]
                exons_list = gene_exons[gene]
                if(strand=="+"):
                    for cont, junction in enumerate(junction_list):
                        # Go over the junction list
                        new_id = gene+"_"+str(junction[0])+"_"+str(junction[1])
                        for cont2,exon in enumerate(exons_list):
                            # Go over the exon list
                            if(int(junction[0])+1==int(exon[1])):
                            # if (int(junction[1]) + 1 == int(exon[0])):
                                # print("Passed 1st condition")
                                #Found 3. Iterate over the rest of the downstream junctions, looking for a 4
                                for cont3,junction2 in enumerate(junction_list[cont+1:]):
                                    # Go over the junction list
                                    new_id2 = gene + "_" + str(junction2[0]) + "_" + str(junction2[1])
                                    for cont4, exon2 in enumerate(exons_list[cont2+1:]):
                                        #Check if there is overlapping and if the junctions extremes are not giving a negative exon
                                        # if (int(junction2[1]) == int(exon2[0]) and (junction[1] - junction2[0])>0):
                                        if (int(junction2[1]) == int(exon2[0]) and (junction2[0] - junction[1]) > 0):
                                            # print("Passed 2nd condition")
                                            #Found 4. Output this result
                                            #Save the junctions associated to this new exon
                                            new_exon_id = gene_chr[gene]+";"+str(junction[1])+";"+str(junction2[0])+";"+gene_strand[gene]
                                            new_junction3_id = gene_chr[gene] + ";" + str(junction[0]) + ";" + str(
                                                junction[1]) + ";" + gene_strand[gene]
                                            new_junction4_id = gene_chr[gene] + ";" + str(junction2[0]) + ";" + str(
                                                junction2[1]) + ";" + gene_strand[gene]
                                            exon_gene[new_exon_id] = gene
                                            length = abs(junction[1] - junction2[0])
                                            exon_lengths[new_exon_id] = length
                                            if(new_exon_id not in exon_junction3):
                                                exon_junction3[new_exon_id] = [new_junction3_id]
                                            else:
                                                #Only add it if its not repeated
                                                if(new_junction3_id not in exon_junction3[new_exon_id]):
                                                    exon_junction3[new_exon_id].append(new_junction3_id)
                                            if(new_exon_id not in exon_junction4):
                                                exon_junction4[new_exon_id] = [new_junction4_id]
                                            else:
                                                #Only add it if its not repeated
                                                if(new_junction4_id not in exon_junction4[new_exon_id]):
                                                    exon_junction4[new_exon_id].append(new_junction4_id)
                else:
                    for cont, junction in enumerate(junction_list):
                        # Go over the junction list
                        new_id = gene+"_"+str(junction[0])+"_"+str(junction[1])
                        for cont2,exon in enumerate(exons_list):
                            # Go over the exon list
                            if(int(junction[1])+1==int(exon[0])):
                            # if (int(junction[1]) + 1 == int(exon[0])):
                                # print("Passed 1st condition")
                                #Found 3. Iterate over the rest of the downstream junctions, looking for a 4
                                for cont3,junction2 in enumerate(junction_list[cont+1:]):
                                    # Go over the junction list
                                    new_id2 = gene + "_" + str(junction2[0]) + "_" + str(junction2[1])
                                    for cont4, exon2 in enumerate(exons_list[cont2+1:]):
                                        #Check if there is overlapping and if the junctions extremes are not giving a negative exon
                                        if (int(junction2[0]) == int(exon2[1]) and (junction[0] - junction2[1])>0):
                                            # print("Passed 2nd condition")
                                            #Found 4. Output this result
                                            #Save the junctions associated to this new exon
                                            new_exon_id = gene_chr[gene]+";"+str(junction2[1])+";"+str(junction[0])+";"+gene_strand[gene]
                                            new_junction3_id = gene_chr[gene] + ";" + str(junction[0]) + ";" + str(
                                                junction[1]) + ";" + gene_strand[gene]
                                            new_junction4_id = gene_chr[gene] + ";" + str(junction2[0]) + ";" + str(
                                                junction2[1]) + ";" + gene_strand[gene]
                                            exon_gene[new_exon_id] = gene
                                            length = abs(junction[0] - junction2[1])
                                            exon_lengths[new_exon_id] = length
                                            if(new_exon_id not in exon_junction3):
                                                exon_junction3[new_exon_id] = [new_junction3_id]
                                            else:
                                                #Only add it if its not repeated
                                                if(new_junction3_id not in exon_junction3[new_exon_id]):
                                                    exon_junction3[new_exon_id].append(new_junction3_id)
                                            if(new_exon_id not in exon_junction4):
                                                exon_junction4[new_exon_id] = [new_junction4_id]
                                            else:
                                                #Only add it if its not repeated
                                                if(new_junction4_id not in exon_junction4[new_exon_id]):
                                                    exon_junction4[new_exon_id].append(new_junction4_id)
            else:
                # logger.info("Gene "+gene+" not in gtf")
                pass

        # Output the values
        logger.info("Output values...")
        output_file = open(output_path, 'w')
        output_file.write("Gene\tJunction_id3\tJunction_id4\tNew_exon\tExon_length\n")
        for x,junctions3_list in exon_junction3.items():
            if(x in exon_junction4):
                junctions4_list = exon_junction4[x]
            else:
                raise Exception("Exon id:"+ x + "not in junction4 list")
            if(max_length>0):
                if(exon_lengths[x]<=max_length):
                    output_file.write(exon_gene[x]+"\t"+",".join(junctions3_list)+"\t"+",".join(junctions4_list)+"\t"+
                                      x +"\t" + str(exon_lengths[x])+"\n")
            else:
                output_file.write(exon_gene[x] + "\t" + ",".join(junctions3_list) + "\t" + ",".join(junctions4_list) + "\t" +
                    x + "\t" + str(exon_lengths[x]) + "\n")
        output_file.close()

        logger.info("Created "+output_path)

        # Remove from the created file those exonizations that are falling in exonic regions
        # Create a bed file from the previous file
        logger.info("Removing exonizations overlapping exons...")
        logger.info("Creating bed file with the exonizations...")
        exonizations = pd.read_table(output_path, delimiter="\t")
        chr = exonizations['New_exon'].apply(lambda x: x.split(";")[0])
        start = exonizations['New_exon'].apply(lambda x: x.split(";")[1])
        end = exonizations['New_exon'].apply(lambda x: x.split(";")[2])
        strand = exonizations['New_exon'].apply(lambda x: x.split(";")[3])
        # Save this variables as bed file
        path1 = "/".join(output_path.split("/")[:-1])
        bed = [("chr", chr), ("start", start), ("end", end), ("id", exonizations['New_exon']),
               ("strand", strand)]
        bed_file = pd.DataFrame.from_items(bed)
        bed_file['score'] = 0
        bed_file.to_csv(path1 + "/exonizations.bed", sep="\t", index=False, header=False)

        # Get the information of the gene and transform it ot bed
        outFile = open(path1 + "/gtf.bed", 'w')
        with open(gtf_path) as f:
            logger.info("Creating bed file from gtf...")
            for line in f:
                # Skip the comments and keep just the exons
                if (not(re.search("#", line)) and re.search("exon",line)):
                    tokens = line.rstrip().split("\t")
                    # Extract just the number of the exon
                    exon = str(tokens[8].split(";")[2].split("\"")[1])
                    transcript_id = str(tokens[8].split(";")[1])[16:31]
                    gene_id = str(tokens[8].split(";")[3].split("\"")[1])
                    outFile.write(tokens[0] + "\t" + str(int(tokens[3]) - 1) + "\t" + tokens[
                        4] + "\t" + gene_id + ":" + transcript_id +
                                  ":" + "exon_" + exon + ":" + str(int(tokens[3]) - 1) + "_" + tokens[4] + ":" + tokens[
                                      6] + "\t" + str(0) + "\t" + tokens[6] + "\n")

        outFile.close()

        # Run interesectBed for obtaining the new exons that are not in coding regions
        logger.info("Running intersectBed...")
        command = "module load BEDTools; intersectBed -v -a " + path1 + "/exonizations.bed -b " + path1 + "/gtf.bed > " + \
                  path1 + "/intersection.bed; module unload module load BEDTools"
        os.system(command)

        # Take from the output file the exons obtained with intersectBed
        logger.info("Filtering exons from final file...")
        intersection = pd.read_table(path1 + "/intersection.bed", delimiter="\t", header=None)
        intersection.columns = ['chr', 'start', 'end', 'exon_id','strand', 'score']
        pos = exonizations["New_exon"].isin(intersection["exon_id"]).as_matrix()
        exonizations_filtered = exonizations.iloc[pos]  # .tolist()
        exonizations_filtered.to_csv(output_path, sep="\t", index=False, header=True)

        # Get the motifs of the splice sites
        # Transform the exons cordinates to 2 bed files: one with the cordinates of the 5' and other with the 3'
        chr = exonizations_filtered["New_exon"].apply(lambda x: x.split(";")[0])
        start = exonizations_filtered["New_exon"].apply(lambda x: str(int(x.split(";")[1])-1))
        start2 = exonizations_filtered["New_exon"].apply(lambda x: str(int(x.split(";")[1])-3))
        end = exonizations_filtered["New_exon"].apply(lambda x: str(int(x.split(";")[2])+1))
        end2 = exonizations_filtered["New_exon"].apply(lambda x: str(int(x.split(";")[2])+3))
        strand = exonizations_filtered["New_exon"].apply(lambda x: x.split(";")[3])

        bed5 = [("chr", chr), ("start", start2), ("end", start), ("strand", strand)]
        bed5_file = pd.DataFrame.from_items(bed5)
        bed5_file['id'] = exonizations_filtered["New_exon"]
        bed5_file['score'] = 0
        #Resort the columns
        cols = bed5_file.columns.tolist()
        cols = cols[0:3] + cols[4:6] + [cols[3]]
        bed5_file = bed5_file[cols]

        bed3 = [("chr", chr), ("start", end), ("end", end2), ("strand", strand)]
        bed3_file = pd.DataFrame.from_items(bed3)
        bed3_file['id'] = exonizations_filtered["New_exon"]
        bed3_file['score'] = 0
        #Resort the columns
        cols = bed3_file.columns.tolist()
        cols = cols[0:3] + cols[4:6] + [cols[3]]
        bed3_file = bed3_file[cols]

        # Save this variables as bed file
        path1 = "/".join(output_path.split("/")[:-1])
        logger.info("Saving bed files...")
        bed5_file.to_csv(path1 + "/bed5.bed", sep="\t", index=False, header=False)
        bed3_file.to_csv(path1 + "/bed3.bed", sep="\t", index=False, header=False)

        # Run getFasta from MoSEA (it needs python2 and BEDTools)
        command1 = "module load Python/2.7.11; module load BEDTools; python2 /genomics/users/juanluis/Software/MoSEA-master/mosea.py getfasta --bedfile " + \
                   path1 + "/bed5.bed --genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa " \
                           "--output " + path1 + "/bed5.fasta"
        os.system(command1)
        command2 = "module load Python/2.7.11; module load BEDTools; python2 /genomics/users/juanluis/Software/MoSEA-master/mosea.py getfasta --bedfile " + \
                   path1 + "/bed3.bed --genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa " \
                           "--output " + path1 + "/bed3.fasta"
        os.system(command2)

        logger.info("Obtaining the splice sites...")
        #Load the motifs and associate it to the exons. Make it a table
        motif5, position5 = [], []
        with open(path1 + "/bed5.fasta") as f:
            for line in f:
                if(line[0]==">"):
                    # strand = line.rstrip()[1]
                    # aux = line.rstrip()[4:-3]
                    # chr = aux.split(":")[0]
                    # start = aux.split(":")[1].split("-")[0]
                    # end = aux.split(":")[1].split("-")[1]
                    # id = chr + ";" + start + ";" + end + ";" + strand
                    id = line.rstrip()[1:-3]
                    position5.append(id)
                else:
                    motif5.append(line.rstrip())

        motif3, position3 = [], []
        with open(path1 + "/bed3.fasta") as f:
            for line in f:
                if(line[0]==">"):
                    # strand = line.rstrip()[1]
                    # aux = line.rstrip()[4:-3]
                    # chr = aux.split(":")[0]
                    # start = aux.split(":")[1].split("-")[0]
                    # end = aux.split(":")[1].split("-")[1]
                    # id = chr + ";" + start + ";" + end + ";" + strand
                    id = line.rstrip()[1:-3]
                    position3.append(id)
                else:
                    motif3.append(line.rstrip())

        #Associate the motifs to the original exonization file
        exonizations_filtered["splice_site5"] = motif5
        exonizations_filtered["splice_site3"] = motif3
        # exonizations_filtered["position5"] = position5
        # exonizations_filtered["position3"] = position3

        exonizations_filtered.to_csv(output_path, sep="\t", index=False, header=True)

        # Remove auxiliary files
        # os.remove(path1 + "/exonizations.bed")
        # os.remove(path1 + "/gtf.bed")
        # os.remove(path1 + "/intersection.bed")

        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)

