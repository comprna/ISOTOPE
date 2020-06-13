"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

select_fasta_candidates: from our list of exonizations with a change in peptide sequence, we are gonna extract
from the file with the fasta sequences of the reference and the exonization, the ones in this list (for running NetMHC)
We're gonna output also the fasta sequences one file per exon.
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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


def select_fasta_candidates(input_list_path, input_sequence_path, output_sequence_path, output_sequence_pieces_path):

    try:
        logger.info("Starting execution")

        # input_list_path = sys.argv[1]
        # input_sequence_path = sys.argv[2]
        # output_sequence_path = sys.argv[3]
        # output_sequence_pieces_path = sys.argv[4]

        # input_list_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/all_exonizations_ORF_filtered_peptide_change.tab"
        ## input_list_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_exonized_junctions_ORF_filtered.tab"
        # input_sequence_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_peptide_sequence.fa"
        # output_sequence_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_peptide_sequence_filtered.fa"
        # output_sequence_pieces_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonization_fasta_files"

        # 1. Load the list of fasta sequences
        logger.info("Loading fasta sequences...")
        sequences_dict = {}
        flag = False
        with open(input_sequence_path) as f:
            for line in f:
                if (not flag and re.search(">", line)):
                    transcript_id_ref = line.rstrip()[1:]
                elif(not flag and not re.search(">", line)):
                    sequence_ref = line.rstrip()
                    flag = True
                elif(flag and re.search(">", line)):
                    transcript_id_ex = line.rstrip()[1:].split("_")[0]
                    id_exonization = line.rstrip()[1:].split("|")[1]
                elif(flag and not re.search(">", line)):
                    sequence_ex = line.rstrip()
                    flag = False
                    #Save the sequences
                    if(transcript_id_ref==transcript_id_ex):
                        sequences_dict[id_exonization] = (sequence_ref,sequence_ex)
                    else:
                        raise Exception("Different transcript ids")

        # # 2. Read the input list with the candidates, outputing the sequences that change peptide
        # # Output in the whole file and by exonization
        # outFile = open(output_sequence_path,"w")
        # cont = 0
        # with open(input_list_path) as f:
        #     next(f)
        #     for line in f:
        #         cont += 1
        #         tokens = line.rstrip().split("\t")
        #         if(tokens[14]=="True"):
        #             if(cont==1):
        #                 new_sample_id = tokens[0].rstrip()
        #                 outFile_individual = open(output_sequence_pieces_path + "/" + new_sample_id + ".fa", "w")
        #             elif(tokens[0].rstrip() != new_sample_id):
        #                 new_sample_id = tokens[0].rstrip()
        #                 #Open a new file and close the previous one
        #                 outFile_individual.close()
        #                 outFile_individual = open(output_sequence_pieces_path+"/"+new_sample_id+".fa", "w")
        #             if(tokens[1] in sequences_dict):
        #                 sequence_ref, sequence_ex = sequences_dict[tokens[1]]
        #                 outFile.write(">"+tokens[1]+"_reference\n")
        #                 outFile.write(sequence_ref+"\n")
        #                 outFile.write(">" + tokens[1] + "_exonization\n")
        #                 outFile.write(sequence_ex + "\n")
        #                 outFile_individual.write(">r_"+new_sample_id+"_"+tokens[15]+"\n")
        #                 outFile_individual.write(sequence_ref+"\n")
        #                 outFile_individual.write(">e_" +new_sample_id+"_"+tokens[15]+"\n")
        #                 outFile_individual.write(sequence_ex + "\n")
        #             else:
        #                 raise Exception(tokens[1]+ " is not in the dictionary!!")

        # 2. Read the input list with the candidates, outputing the sequences that change peptide
        # Output in the whole file and by exonization
        logger.info("Filtering output...")
        cont = 1
        outFile = open(output_sequence_path,"w")
        with open(input_list_path) as f:
            header = next(f).rstrip().split("\t")
            for line in f:
                tokens = line.rstrip().split("\t")
                New_Exon_pos = header.index("New_exon")
                Sample_id_pos = header.index("Sample_id")
                Index_pos = header.index("Index")
                if(tokens[New_Exon_pos] in sequences_dict):
                    new_sample_id = tokens[Sample_id_pos].rstrip()
                    outFile_individual = open(output_sequence_pieces_path + "/" + tokens[Index_pos] + ".fa", "w")
                    sequence_ref, sequence_ex = sequences_dict[tokens[New_Exon_pos]]
                    outFile.write(">"+tokens[New_Exon_pos]+"_reference_"+new_sample_id+"\n")
                    outFile.write(sequence_ref+"\n")
                    outFile.write(">" + tokens[New_Exon_pos] + "_exonization_"+new_sample_id+"\n")
                    outFile.write(sequence_ex + "\n")
                    outFile_individual.write(">r_"+new_sample_id+"_"+tokens[Index_pos]+"\n")
                    outFile_individual.write(sequence_ref+"\n")
                    outFile_individual.write(">e_" +new_sample_id+"_"+tokens[Index_pos]+"\n")
                    outFile_individual.write(sequence_ex + "\n")
                    outFile_individual.close()
                    cont += 1
                else:
                    raise Exception(tokens[Sample_id_pos]+ " is not in the dictionary!!")

        outFile.close()

        logger.info("Saved "+output_sequence_path)
        logger.info("Saved also same file in pieces in "+output_sequence_pieces_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)

