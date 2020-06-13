"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

associate_gene_ids.py: given a list of introns, we will check in which genes
are falling
arg[1]: list of introns
arg[2]: gtf
arg[3]: list of introns - with the corresponding genes
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys
import pandas as pd
import os

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


def IR_kma_associate_gene_ids(introns_path, gtf_path, output_path):


    try:
        logger.info("Starting execution")

        # 1. Load the introns file
        logger.info("Reading introns file...")
        introns_file = pd.read_table(introns_path, delimiter="\t")
        chr = introns_file['Event_id'].apply(lambda x: x.split(":")[0])
        start = introns_file['Event_id'].apply(lambda x: x.split(":")[1].split("-")[0])
        end = introns_file['Event_id'].apply(lambda x: x.split(":")[1].split("-")[1].split("(")[0])
        strand = introns_file['Event_id'].apply(lambda x: x.split("(")[1].split(")")[0])
        input_bedtools = pd.DataFrame(
            {'chr': chr.tolist(), 'start': start.tolist(), 'end': end.tolist(), 'id': introns_file['Event_id'].tolist(),
             'reads': 0, 'strand': strand.tolist()})

        columns = ['chr','start','end','id','reads','strand']
        input_bedtools2 = input_bedtools[columns]

        # 2. Save input_bedtools for running bedtools
        output_aux = "/".join(output_path.split("/")[:-1])+"/input_bedtools.bed"
        input_bedtools2.to_csv(output_aux, sep="\t", index=False, header=False)
        logger.info("Running bedtools...")
        command = "module load BEDTools; intersectBed -wao -a " + output_aux + " -b " + gtf_path + " > " + \
                  "/".join(output_path.split("/")[:-1]) + "/input_bedtools2.bed; module unload BEDTools;"
        os.system(command)

        # 3. Read the created file, saving all the genes associated to each intron
        logger.info("Formatting bedtools output...")
        dict_introns1, dict_introns2, dict_intron_strand = {}, {}, {}
        dict_introns_flag_5ss, dict_introns_flag_3ss = {}, {}
        # cont = 0
        for line in open("/".join(output_path.split("/")[:-1]) + "/input_bedtools2.bed"):
            tokens = line.rstrip().split("\t")
            intron_id = tokens[3]
            # print(str(cont))
            # cont += 1
            #Check if the intron overlaps on any exon
            if(int(tokens[15])!=0):  # (int(tokens[1])+40==int(tokens[9]) and int(tokens[2])-39)==int(tokens[10])):
                if(intron_id in dict_introns_flag_5ss and intron_id in dict_introns_flag_3ss):
                    continue
                elif(intron_id not in dict_introns_flag_5ss and int(tokens[1])+40==int(tokens[10])):
                    dict_introns_flag_5ss[intron_id] = True
                elif(intron_id not in dict_introns_flag_3ss and int(tokens[2])-39==int(tokens[9])):
                    dict_introns_flag_3ss[intron_id] = True

                if(intron_id in dict_introns_flag_5ss and intron_id in dict_introns_flag_3ss):
                    gene_name = tokens[14].split("gene_name")[1].split(";")[0].replace("\"", "")[1:]
                    gene_id = tokens[14].split("gene_id")[1].split(";")[0].replace("\"", "")[1:]
                    #Get the gene and the strands
                    if(intron_id not in dict_introns1):
                        dict_introns1[intron_id] = [gene_name]
                    elif(gene_name not in dict_introns1[intron_id]):
                        dict_introns1[intron_id].append(gene_name)
                    if(intron_id not in dict_introns2):
                        dict_introns2[intron_id] = [gene_id]
                    elif(gene_id not in dict_introns2[intron_id]):
                        dict_introns2[intron_id].append(gene_id)
                    if(intron_id not in dict_intron_strand):
                        dict_intron_strand[intron_id] = tokens[12]
                    elif(dict_intron_strand[intron_id] != tokens[12]):
                        logger.error("Different strands associated to the same intron")
                    # raise Exception ("Different strands associated to the same intron")

        # 4. Read the initial file assigning the genes to each intron id
        # Read the events file
        logger.info("Processing final file...")
        outFile = open(output_path, 'w')
        # outFile.write("IR\tGene_name\tGene_id\tSample_id\tZscore\tp_value\tFDR\tTPM_mut\tTPM_nonmut\n")
        with open(introns_path) as e:
            header = next(e).rstrip().split("\t")
            outFile.write(header[0]+"\tGene_name\tGene_id\t"+"\t".join(header[1:])+"\n")
            for line in e:
                tokens = line.rstrip().split("\t")
                intron_id = tokens[0]
                strand = intron_id.split("(")[1].split(")")[0]
                if(intron_id in dict_introns1):
                    gene_names = ",".join(dict_introns1[intron_id])
                else:
                    gene_names = "No gene"
                if(intron_id in dict_introns2):
                    gene_ids = ",".join(dict_introns2[intron_id])
                else:
                    gene_ids = "No gene"
                #Extract the strand
                if(intron_id in dict_intron_strand):
                    if(strand != dict_intron_strand[intron_id]):
                        intron_id = tokens[0].split("(")[0]+"("+dict_intron_strand[intron_id]+"):"+tokens[0].split(":")[-1]
                outFile.write(intron_id + "\t" + gene_names + "\t" + gene_ids + "\t" + "\t".join(tokens[1:]) + "\n")

        # 5. Remove auxiliary files before closing
        # os.remove("/".join(output_path.split("/")[:-1]) + "/input_bedtools.bed")
        # os.remove("/".join(output_path.split("/")[:-1]) + "/input_bedtools2.bed")

        logger.info("Created "+output_path)
        outFile.close()
        logger.info("Done. Exiting program.")

    except Exception as error:
        print('ERROR: ' + repr(error))
        print("Aborting execution")
        sys.exit(1)
