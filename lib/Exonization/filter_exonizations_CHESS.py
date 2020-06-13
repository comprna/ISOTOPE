"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_exonizations_CHESS: from our list opf exonizations, we are gonna check if this events appears in the CHESS
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


def filter_exonizations_CHESS(exonizations_path, CHESS_SE_path, output_path):

    try:
        logger.info("Starting execution")

        # 1. Load the SE CHESS events, save the central exon
        logger.info("Loading CHESS db...")
        CHESS_event = {}
        with open(CHESS_SE_path) as f:
            header = next(f).rstrip()
            for line in f:
                tokens = line.rstrip().rsplit("\t")
                event_id = tokens[2]
                #Extract the chr and coordinates
                chr = event_id.split(":")[1]
                start = event_id.split(":")[2].split("-")[1]
                end = event_id.split(":")[3].split("-")[0]
                strand = event_id.split(":")[4]
                event_id_f = chr + ";" + start + ";" + end + ";" + strand
                #Save it in a dictionary
                if(event_id_f not in CHESS_event):
                    CHESS_event[event_id_f] = 0
                else:
                    pass
                    # logger.info("Repeated event "+event_id_f)

        # 2. Load the exonizations, checking which ones are in CHESS
        logger.info("Checking exonizations...")
        cont = 0
        outfile = open(output_path,"w")
        with open(exonizations_path) as f:
            header = next(f).rstrip().rsplit("\t")
            outfile.write("\t".join(header)+"\n")
            New_exon_pos = header.index("New_exon")
            Junction_id3_pos = header.index("Junction_id3")
            Junction_id4_pos = header.index("Junction_id4")
            for line in f:
                cont += 1
                # print(str(cont))
                tokens = line.rstrip().rsplit("\t")
                central_exon = tokens[New_exon_pos]
                #Get the cordinates and the chr
                chr = central_exon.split(";")[0]
                start_central_exon = int(central_exon.split(";")[1])
                end_central_exon = int(central_exon.split(";")[2])+1
                strand = central_exon.split(";")[3]
                # Junction_id3 = tokens[Junction_id3_pos]
                # #Take the end
                # end_Junction_id3 = int(Junction_id3.split(";")[2])
                # Junction_id4 = tokens[Junction_id4_pos]
                # # Take the start
                # start_Junction_id4 = int(Junction_id4.split(";")[1])
                # Assemble the event
                event_id_f = chr + ";" + str(start_central_exon) + ";" + str(end_central_exon) + ";" + strand
                # Look for it in the CHESS db.
                if(event_id_f in CHESS_event):
                    #Dont plot the line and go with the next line
                    logger.info(str(cont) + ": found event in CHESS")
                    continue
                else:
                    outfile.write(line)

                    # Check if event
                # for key, values in CHESS_event.items():
                #     start_central = int(key.split(":")[1].split("-")[1])
                #     end_central = int(key.split(":")[2].split("-")[0])
                #     end_junction = int(key.split(":")[1].split("-")[0])
                #     start_junction = int(key.split(":")[2].split("-")[1])
                #     if(start_central-5<int(start_central_exon)<start_central+5 and
                #                        end_central-5<int(end_central_exon)<end_central+5 and
                #                        end_junction - 5 < int(end_Junction_id3) < end_junction + 5 and
                #                        start_junction - 5 < int(start_Junction_id4) < start_junction + 5):
                #         logger.info("Found event in CHESS "+event_id_f)

        outfile.close()
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
