"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

select_fasta_candidates: from our list of neoskipping with a change in peptide sequence, we are gonna extract
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

        # 2. Read the input list with the candidates, outputing the sequences that change peptide and are not falling into NMD
        # Output in the whole file and by neoskipping
        logger.info("Filtering output...")
        cont = 0
        outFile = open(output_sequence_path,"w")
        with open(input_list_path) as f:
            header = next(f).rstrip().split("\t")
            New_Exon_pos = header.index("Neoskipping_junction")
            Sample_id_pos = header.index("Sample_id")
            Peptide_change_pos = header.index("Peptide_change")
            NMD_pos = header.index("NMD")
            Index_pos = header.index("Index")
            for line in f:
                cont += 1
                tokens = line.rstrip().split("\t")
                # if(len(tokens)==19):
                # if(tokens[Peptide_change_pos]=="True" and tokens[NMD_pos]=="False"):
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
                else:
                    raise Exception(tokens[New_Exon_pos]+ " is not in the dictionary!!")

        outFile.close()

        logger.info("Saved "+output_sequence_path)
        logger.info("Saved also same file in pieces in "+output_sequence_pieces_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
