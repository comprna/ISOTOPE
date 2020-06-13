"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

select_fasta_candidates: from our list of with a change in peptide sequence, we are gonna extract
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
        logger.info("Processing fasta sequences file...")
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
                    id_exonization = "|".join(line.rstrip()[1:].split("|")[1:])
                elif(flag and not re.search(">", line)):
                    sequence_ex = line.rstrip()
                    flag = False
                    #Save the sequences
                    if(transcript_id_ref==transcript_id_ex):
                        sequences_dict[id_exonization] = (sequence_ref,sequence_ex)
                    else:
                        raise Exception("Different transcript ids")

        # 2. Read the input list with the candidates, outputing the sequences that change peptide
        # Output in the whole file and by exonization
        cont = 1
        outFile = open(output_sequence_path,"w")
        with open(input_list_path) as f:
            header = next(f).rstrip().split("\t")
            for line in f:
                tokens = line.rstrip().split("\t")
                Canonical_Exon_pos = header.index("Canonical_Exon")
                Alt_Exon_pos = header.index("Alt_Exon_id")
                Sample_id_pos = header.index("Sample_id")
                id = tokens[Canonical_Exon_pos]+"|"+tokens[Alt_Exon_pos]
                if(id in sequences_dict):
                    new_sample_id = tokens[Sample_id_pos].rstrip()
                    outFile_individual = open(output_sequence_pieces_path + "/" + str(cont) + ".fa", "w")
                    sequence_ref, sequence_ex = sequences_dict[id]
                    outFile.write(">"+id+"_reference_"+new_sample_id+"\n")
                    outFile.write(sequence_ref+"\n")
                    outFile.write(">" + id+ "_exonization_"+new_sample_id+"\n")
                    outFile.write(sequence_ex + "\n")
                    outFile_individual.write(">r_"+new_sample_id+"_"+str(cont)+"\n")
                    outFile_individual.write(sequence_ref+"\n")
                    outFile_individual.write(">e_" +new_sample_id+"_"+str(cont)+"\n")
                    outFile_individual.write(sequence_ex + "\n")
                    outFile_individual.close()
                else:   #There is no sequence generated for this case. Skip it
                    logger.info(id+ " is not in the dictionary")
                cont += 1

        outFile.close()

        logger.info("Saved "+output_sequence_path)
        logger.info("Saved also same file in pieces in "+output_sequence_pieces_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)

