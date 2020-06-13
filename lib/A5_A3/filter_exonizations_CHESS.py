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


def filter_exonizations_CHESS(exonization_path, CHESS_A5_path, CHESS_A3_path, output_path):

    try:
        logger.info("Starting execution")

        # 1. Load the A5/A3 CHESS events,
        logger.info("Loading CHESS db...")
        CHESS_event_A5 = {}
        with open(CHESS_A5_path) as f:
            header = next(f).rstrip()
            for line in f:
                tokens = line.rstrip().rsplit("\t")
                event_id = tokens[2]
                #Extract the chr and coordinates
                chr = event_id.split(":")[1]
                strand = event_id.split(":")[4]
                #Depending on the strand, the coordiantes change
                if(strand=="+"):
                    e2 = event_id.split(":")[2].split("-")[0]
                    s3 = event_id.split(":")[2].split("-")[1]
                    e1 = event_id.split(":")[3].split("-")[0]
                    event_id_f = chr + ";" + e1 + ";" + e2 + ";" + s3 + ";" + strand
                else:
                    s3 = event_id.split(":")[3].split("-")[1]
                    s2 = event_id.split(":")[2].split("-")[1]
                    e1 = event_id.split(":")[2].split("-")[0]
                    event_id_f = chr + ";" + e1 + ";" + s2 + ";" + s3 + ";" + strand
                #Save it in a dictionary
                if(event_id_f not in CHESS_event_A5):
                    CHESS_event_A5[event_id_f] = 0
                else:
                    pass
                    # logger.info("Repeated event "+event_id_f)

        CHESS_event_A3 = {}
        with open(CHESS_A3_path) as f:
            header = next(f).rstrip()
            for line in f:
                tokens = line.rstrip().rsplit("\t")
                event_id = tokens[2]
                #Extract the chr and coordinates
                chr = event_id.split(":")[1]
                strand = event_id.split(":")[4]
                #Depending on the strand, the coordiantes change
                if(strand=="+"):
                    e1 = event_id.split(":")[2].split("-")[0]
                    s2 = event_id.split(":")[2].split("-")[1]
                    s3 = event_id.split(":")[3].split("-")[1]
                    event_id_f = chr + ";" + e1 + ";" + s2 + ";" + s3 + ";" + strand
                else:
                    s3 = event_id.split(":")[3].split("-")[1]
                    e2 = event_id.split(":")[2].split("-")[0]
                    e1 = event_id.split(":")[3].split("-")[0]
                    event_id_f = chr + ";" + e1 + ";" + e2 + ";" + s3 + ";" + strand
                #Save it in a dictionary
                if(event_id_f not in CHESS_event_A3):
                    CHESS_event_A3[event_id_f] = 0
                else:
                    pass
                    # logger.info("Repeated event "+event_id_f)

        # 2. Load the exonization, checking which ones are in CHESS
        logger.info("Checking exonization...")
        cont = 0
        outfile = open(output_path,"w")
        with open(exonization_path) as f:
            header = next(f).rstrip().rsplit("\t")
            outfile.write("\t".join(header)+"\n")
            Canonical_Junction_pos = header.index("Canonical_Junction_id")
            Alt_Junction_pos = header.index("Alt_Junction_id")
            splice_site_pos = header.index("splice_site_type")
            strand_pos = header.index("strand")
            for line in f:
                cont += 1
                logger.info(str(cont)+"...")
                tokens = line.rstrip().rsplit("\t")
                canonical_junction = tokens[Canonical_Junction_pos]
                alt_junction = tokens[Alt_Junction_pos]
                strand = tokens[strand_pos]
                splice_site = tokens[splice_site_pos]
                if(splice_site=="New_acceptor"):    #A3
                    if (strand == "+"):
                        chr = canonical_junction.split(";")[0]
                        c1 = int(canonical_junction.split(";")[2])
                        c2 = int(alt_junction.split(";")[1]) + 1
                        c3 = int(alt_junction.split(";")[2])
                        if (c1 < c3):
                            # Assemble the event
                            event_id_f = chr + ";" + str(c2) + ";" + str(c1) + ";" + str(c3) + ";" + strand
                        else:
                            # Assemble the event
                            event_id_f = chr + ";" + str(c2) + ";" + str(c3) + ";" + str(c1) + ";" + strand
                    else:
                        chr = canonical_junction.split(";")[0]
                        c1 = int(canonical_junction.split(";")[1]) + 1
                        c2 = int(canonical_junction.split(";")[2])
                        c3 = int(alt_junction.split(";")[1]) + 1
                        if (c1 < c3):
                            # Assemble the event
                            event_id_f = chr + ";" + str(c1) + ";" + str(c3) + ";" + str(c2) + ";" + strand
                        else:
                            # Assemble the event
                            event_id_f = chr + ";" + str(c3) + ";" + str(c1) + ";" + str(c2) + ";" + strand
                    # Look for it in the CHESS db.
                    if (event_id_f in CHESS_event_A3):
                        # Dont plot the line and go with the next line
                        logger.info(str(cont) + ": found A3 event in CHESS: " + event_id_f)
                        continue
                    else:
                        outfile.write(line)
                    # # Check if event is in CHESS_event_A3
                    # c1 = event_id_f.split(";")[1]
                    # c2 = event_id_f.split(";")[2]
                    # c3 = event_id_f.split(";")[3]
                    # for key, values in CHESS_event_A3.items():
                    #     c1_aux = int(key.split(";")[1])
                    #     c2_aux = int(key.split(";")[2])
                    #     c3_aux = int(key.split(";")[3])
                    #     if (c1_aux - 5 < int(c1) < c1_aux + 5 and c2_aux - 5 < int(c2) < c2_aux + 5
                    #         and c3_aux - 5 < int(c3) < c3_aux + 5):
                    #         logger.info("Found A3 event in CHESS " + key + "|" + event_id_f)
                    #         break
                else:  # A5
                    if(strand=="+"):
                        chr = canonical_junction.split(";")[0]
                        c1 = int(canonical_junction.split(";")[1])+1
                        c2 = int(canonical_junction.split(";")[2])
                        c3 = int(alt_junction.split(";")[1])+1
                        if(c1<c3):
                            # Assemble the event
                            event_id_f = chr + ";" + str(c1) + ";" + str(c3) + ";" + str(c2) + ";" + strand
                        else:
                            # Assemble the event
                            event_id_f = chr + ";" + str(c3) + ";" + str(c1) + ";" + str(c2) + ";" + strand
                    else:
                        chr = canonical_junction.split(";")[0]
                        c1 = int(canonical_junction.split(";")[2])
                        c2 = int(alt_junction.split(";")[1])+1
                        c3 = int(alt_junction.split(";")[2])
                        if(c1<c3):
                            # Assemble the event
                            event_id_f = chr + ";" + str(c2) + ";" + str(c1) + ";" + str(c3) + ";" + strand
                        else:
                            # Assemble the event
                            event_id_f = chr + ";" + str(c2) + ";" + str(c3) + ";" + str(c1) + ";" + strand
                    # Look for it in the CHESS db.
                    if (event_id_f in CHESS_event_A5):
                        # Dont plot the line and go with the next line
                        logger.info(str(cont) + ": found A3 event in CHESS: " + event_id_f)
                        continue
                    else:
                        outfile.write(line)
                    # # Check if event is in CHESS_event_A5
                    # c1 = event_id_f.split(";")[1]
                    # c2 = event_id_f.split(";")[2]
                    # c3 = event_id_f.split(";")[3]
                    # for key, values in CHESS_event_A5.items():
                    #     c1_aux = int(key.split(";")[1])
                    #     c2_aux = int(key.split(";")[2])
                    #     c3_aux = int(key.split(";")[3])
                    #     if (c1_aux - 5 < int(c1) < c1_aux + 5 and c2_aux - 5 < int(c2) < c2_aux + 5
                    #         and c3_aux - 5 < int(c3) < c3_aux + 5):
                    #         logger.info("Found A5 event in CHESS " + key + "|" + event_id_f)
                    #         break

        outfile.close()
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
