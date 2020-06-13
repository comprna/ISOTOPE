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


def filter_neoskipping_CHESS(neoskipping_path, CHESS_SE_path, output_path):

    try:
        logger.info("Starting execution")

        # 1. Load the SE CHESS events, save the flanking exons
        logger.info("Loading CHESS db...")
        CHESS_event = {}
        with open(CHESS_SE_path) as f:
            header = next(f).rstrip()
            for line in f:
                tokens = line.rstrip().rsplit("\t")
                event_id = tokens[2]
                #Extract the chr and coordinates
                chr = event_id.split(":")[1]
                flanking_5 = event_id.split(":")[2].split("-")[0]
                flanking_3 = event_id.split(":")[3].split("-")[1]
                strand = event_id.split(":")[4]
                event_id_f = chr + ";" + flanking_5 + ";" + flanking_3 + ";" + strand
                #Save it in a dictionary
                if(event_id_f not in CHESS_event):
                    CHESS_event[event_id_f] = 0
                else:
                    pass
                    # logger.info("Repeated event "+event_id_f)

        # 2. Load the neoskipping, checking which ones are in CHESS
        logger.info("Checking neoskipping...")
        cont = 0
        outfile = open(output_path,"w")
        with open(neoskipping_path) as f:
            header = next(f).rstrip().rsplit("\t")
            outfile.write("\t".join(header)+"\n")
            Neoskipping_junction_pos = header.index("Neoskipping_junction")
            # Junction_id3_pos = header.index("Junction_id3")
            # Junction_id4_pos = header.index("Junction_id4")
            for line in f:
                cont += 1
                # print(str(cont))
                tokens = line.rstrip().rsplit("\t")
                flanking_junction = tokens[Neoskipping_junction_pos]
                #Get the cordinates and the chr
                chr = flanking_junction.split(";")[0]
                start_flanking_junction = int(flanking_junction.split(";")[1])+1
                end_flanking_junction = int(flanking_junction.split(";")[2])
                strand = flanking_junction.split(";")[3]
                # Junction_id3 = tokens[Junction_id3_pos]
                #Take the end
                # end_Junction_id3 = int(Junction_id3.split(";")[2])
                # Junction_id4 = tokens[Junction_id4_pos]
                # # Take the start
                # start_Junction_id4 = int(Junction_id4.split(";")[1])
                # Assemble the event
                event_id_f = chr + ";" + str(start_flanking_junction) + ";" + str(end_flanking_junction) + ";" + strand
                # Look for it in the CHESS db.
                if(event_id_f in CHESS_event):
                    #Dont plot the line and go with the next line
                    logger.info(str(cont) + ": found event in CHESS")
                    continue
                else:
                    outfile.write(line)

                # Check if event
                # for key, values in CHESS_event.items():
                #     # start_central = int(key.split(":")[1].split("-")[1])
                #     # end_central = int(key.split(":")[2].split("-")[0])
                #     start = int(key.split(";")[1])
                #     end = int(key.split(";")[2])
                #     if(start-5<int(start_flanking_junction)<start+5 and end-5<int(end_flanking_junction)<end+5):
                #         logger.info("Found event in CHESS "+event_id_f)
                #         break

        outfile.close()
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)