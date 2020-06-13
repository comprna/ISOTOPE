"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_neoskipping_CHESS: from our list opf neoskipping, we are gonna check if this events appears in the CHESS
annotation. CHESS is a db that obtains all expressed transcripts in GTEX data. We have processed the info from CHESS
with SUPPA and obtained the events. If we detect any SE in CHESS overlapping any exonization, we will remove it
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


def filter_IR_CHESS(IR_path, CHESS_RI_path, output_path):

    try:
        logger.info("Starting execution")

        # IR_path = sys.argv[1]
        # CHESS_RI_path = sys.argv[2]
        # output_path = sys.argv[3]

        # IR_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered.tab"
        # CHESS_RI_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered2.tab"

        # 1. Load the RI CHESS events, save the flanking exons
        logger.info("Loading CHESS db...")
        CHESS_event = {}
        with open(CHESS_RI_path) as f:
            header = next(f).rstrip()
            for line in f:
                tokens = line.rstrip().rsplit("\t")
                event_id = tokens[2]
                #Extract the chr and coordinates
                chr = event_id.split(":")[1]
                e1 = event_id.split(":")[3].split("-")[0]
                s2 = event_id.split(":")[3].split("-")[1]
                strand = event_id.split(":")[5]
                event_id_f = chr + ";" + e1 + ";" + s2 + ";" + strand
                #Save it in a dictionary
                if(event_id_f not in CHESS_event):
                    CHESS_event[event_id_f] = 0
                else:
                    pass
                    # logger.info("Repeated event "+event_id_f)

        # 2. Load the IR, checking which ones are in CHESS
        logger.info("Checking IR...")
        cont = 0
        outfile = open(output_path,"w")
        with open(IR_path) as f:
            header = next(f).rstrip().rsplit("\t")
            outfile.write("\t".join(header)+"\n")
            IR_pos = header.index("Event_id")
            for line in f:
                cont += 1
                # print(str(cont))
                tokens = line.rstrip().rsplit("\t")
                IR = tokens[IR_pos]
                #Get the cordinates and the chr
                chr_ir = IR.split(":")[0]
                start_ir = int(IR.split(":")[1].split("-")[0])+40
                end_ir = int(IR.split(":")[1].split("-")[1][:-3])-39
                strand = IR.split("(")[1].split(")")[0]
                event_id_f = chr_ir + ";" + str(start_ir) + ";" + str(end_ir) + ";" + strand
                # Look for it in the CHESS db.
                if(event_id_f in CHESS_event):
                    #Dont plot the line and go with the next line
                    logger.info(str(cont) + ": found event in CHESS")
                    continue
                else:
                    outfile.write(line)
                # # Check if event
                # for key, values in CHESS_event.items():
                #     chr = key.split(";")[0]
                #     start = int(key.split(";")[1])
                #     end = int(key.split(";")[2])
                #     if(chr==chr_ir and start_ir<start and end<end_ir):
                #         logger.info("Found event in CHESS "+event_id_f+"|"+key)
                #         exit(0)
                #         break

        outfile.close()
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
