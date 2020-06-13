"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

format_to_SPADA: from the list obtained, format the output for being run with SPADA
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


def format_to_SPADA(input_path1, input_path2, input_path3, input_path4, output_path1, output_path2, output_path3):

    try:
        logger.info("Starting execution")

        # input_path1 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_ORF.tab"
        # input_path2 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_ORF_sequences.tab"
        # input_path3 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_Interpro.tab"
        # input_path4 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_IUPred.tab"
        # output_path1 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_formatted_SPADA.tab"
        # output_path2 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_formatted_SPADA.fasta"
        # output_path3 = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v3/all_neoskipping_formatted_SPADA_features.tab"

        # 1. Load the neoskipping with the gene associated
        peptide_ref, peptide_change = {}, {}
        cont = 0
        with open(input_path2) as f:
            logger.info("Loading sequences file...")
            next(f)
            for line in f:
                cont += 1
                # print(str(cont))
                tokens = line.rstrip().split("\t")
                if(len(tokens)!=1):
                    index = tokens[0]
                    pep_ref = tokens[3]
                    pep_change = tokens[4]
                    if(index not in peptide_ref):
                        peptide_ref[index] = pep_ref
                    else:
                        raise Exception("Repeated index " + str(index))
                    if (index not in pep_change):
                        peptide_change[index] = pep_change
                    else:
                        raise Exception("Repeated index " + str(index))

        # 2. Go over the file formatting the info
        output_file1 = open(output_path1,"w")
        # output_file.write("Sample_id\tCancer_type\tGene\tMutation_id\tTranscript_id\tReference\tAlternative\n")
        output_file1.write("GeneID\tControl_transcript\tCase_transcript\tSamples\n")
        output_file2 = open(output_path2,"w")
        output_file3 = open(output_path3,"w")
        output_file3.write("transcript\tfeatureType\tfeature_id\tstart\tend\n")
        with open(input_path1) as f:
            logger.info("Formatting the input file...")
            header = next(f).rstrip().split("\t")
            sample_ID_pos = header.index("Sample_id")
            transcript_pos = header.index("Transcript_id")
            gene_pos = header.index("Gene_id")
            exonization_pos = header.index("Neoskipping_junction")
            index_pos = header.index("Index")
            peptide_change_pos = header.index("Peptide_change")
            NMD_pos = header.index("NMD")
            for line in f:
                tokens = line.rstrip().split("\t")
                if(len(tokens)==18):
                    #Save the lines that changes the peptide and is not falling in NMD
                    if(tokens[peptide_change_pos]=="True" and tokens[NMD_pos]=="False"):
                        sample_ID = tokens[sample_ID_pos]
                        transcript = tokens[transcript_pos]
                        gene = tokens[gene_pos]
                        exonization = tokens[exonization_pos]
                        index = tokens[index_pos]
                        if(index in peptide_ref and index in peptide_change):
                            #Remove the * from the sequences
                            reference = peptide_ref[index].replace("*","")
                            alternative = peptide_change[index].replace("*","")
                        else:
                            raise Exception("Index " + str(index) + " missing")
                        #Save the values
                        # output_file1.write(sample_ID+"\tSCLC\t"+gene+"\t"+exonization+"\t"+transcript+"\t"+reference+"\t"+alternative+"\n")
                        output_file1.write(gene+"\t"+transcript+"\t"+exonization+"\t"+sample_ID+"\n")
                        output_file2.write(">"+transcript+"\n")
                        output_file2.write(reference+"\n")
                        output_file2.write(">"+exonization+"\n")
                        output_file2.write(alternative+"\n")

        # 3. Format also the IUPred and Interpro predictions
        with open(input_path3) as f:
            logger.info("Formatting the Interpro predictions...")
            for line in f:
                tokens = line.rstrip().split("\t")
                transcript = tokens[0]
                featureType = tokens[3]
                feature_id = tokens[4]
                start = tokens[6]
                end = tokens[7]
                if(featureType=="Pfam" or featureType=="ProSiteProfiles"):
                    #Save this line to the output
                    output_file3.write(transcript+"\t"+featureType+"\t"+feature_id+"\t"+start+"\t"+end+"\n")

        with open(input_path4) as f:
            logger.info("Formatting the IUPred predictions...")
            next(f)
            for line in f:
                #Save all lines
                output_file3.write(line)

        output_file1.close()
        output_file2.close()
        output_file3.close()
        logger.info("Saved "+output_path1)
        logger.info("Saved "+output_path2)
        logger.info("Saved "+output_path3)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)

