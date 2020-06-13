"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

check_mutations_nearby: check if in the neoskipping there are mutations nearby
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

def get_associated_samples(exon,dict):
    if(exon in dict):
        return dict[exon]
    else:
        return ""

def check_coincidences(line):
    samples = line["mut_samples"]
    aux2 = line["Sample_id"].rstrip()
    for x in samples:
        aux = x.split(":")[0]
        if(aux2==aux):
            return True
    return False

def get_offset_A5(line):
    if(line[5]=="+"):
        offset_A5 = line[8]-line["original_start"]
    else:
        offset_A5 = line[8]-line["original_end"]
    return offset_A5

def get_offset_A3(line):
    if(line[5]=="+"):
        offset_A3 = line[8]-line["original_end"]
    else:
        offset_A3 = line[8]-line["original_start"]
    return offset_A3

def check_mutations_nearby(neoskipping_path, mutations_path, window, output_path):

    try:
        logger.info("Starting execution")

        #Load the neoskipping
        neoskipping = pd.read_table(neoskipping_path, delimiter="\t", )
        # Create a bed file
        new_exon_unique = neoskipping.Neoskipping_junction.unique()
        chr = list(map(lambda x: x.split(";")[0],new_exon_unique))
        #Add to the start and end the window shift
        start = list(map(lambda x: int(x.split(";")[1])-window,new_exon_unique))
        end = list(map(lambda x: int(x.split(";")[2])+window,new_exon_unique))
        strand = list(map(lambda x: x.split(";")[3],new_exon_unique))
        # Save this variables as bed file
        path1 = "/".join(output_path.split("/")[:-1])
        bed = [("chr", chr), ("start", start), ("end", end), ("id", new_exon_unique),("score", 0),("strand", strand)]
        bed_file = pd.DataFrame.from_items(bed)
        # bed_file['score'] = 0
        bed_file.to_csv(path1 + "/neoskipping.bed", sep="\t", index=False, header=False)

        # Run interesectBed for obtaining the new exons that are not in coding regions
        logger.info("Running intersectBed...")
        command = "module load BEDTools; intersectBed -wao -a " + path1 + "/neoskipping.bed -b " + mutations_path + " > " + \
                  path1 + "/intersection_mutations.bed; module unload module load BEDTools"
        os.system(command)

        # Take from the output file the exons obtained with intersectBed
        logger.info("Obtain neoskipping overlapping with mutations...")
        intersection = pd.read_table(path1 + "/intersection_mutations.bed", delimiter="\t", header=None)
        # intersection.columns = ['chr', 'start', 'end', 'exon_id','strand', 'score']
        # Get the rows with a 1 in the last column
        intersection_f = intersection.loc[intersection.iloc[:,21] == 1]
        original_start = intersection_f[3].apply(lambda x: int(x.split(";")[1]))
        original_end = intersection_f[3].apply(lambda x: int(x.split(";")[2]))
        intersection_f["original_start"] = original_start
        intersection_f["original_end"] = original_end
        offset_A5 = intersection_f.apply(lambda x: get_offset_A5(x),axis=1)
        offset_A3 = intersection_f.apply(lambda x: get_offset_A3(x),axis=1)
        #Extract the samples per exon that has associated mutation
        #Generate also bed tracks for the visualization of the mutations and the neoskipping
        exon_samples, exon_offset_A5, exon_offset_A3 = {},{},{}
        mut_id = intersection_f.apply(lambda x: x[10]+":"+x[11]+">"+x[12], axis=1)
        intersection_f["mut_id"] = mut_id
        for i in range(0,len(intersection_f.index)):
            exon = intersection_f.iloc[i, 3]
            mut_id_aux = intersection_f.iloc[i].loc["mut_id"]
            if(exon not in exon_samples):
                exon_samples[exon] = [mut_id_aux]
                exon_offset_A5[exon] = [offset_A5.iloc[i]]
                exon_offset_A3[exon] = [offset_A3.iloc[i]]
            else:
                exon_samples[exon].append(mut_id_aux)
                exon_offset_A5[exon].append(offset_A5.iloc[i])
                exon_offset_A3[exon].append(offset_A3.iloc[i])

        #Generate bedtracks of the mutations
        logger.info("Generating bedtracks...")
        # event_info_bed = pd.DataFrame({'chr': intersection_f.iloc[:,6].tolist(), 'start': intersection_f.iloc[:,7].tolist(),
        #                                'end': intersection_f.iloc[:,8].tolist(), 'id': intersection_f['mut_id'].tolist(),
        #                                'score': 0, 'strand': intersection_f.iloc[:,4].tolist()})
        bed2 = [('chr', intersection_f.iloc[:,6].tolist()), ('start', intersection_f.iloc[:,7].tolist()),
                ('end', intersection_f.iloc[:,8].tolist()), ('id', intersection_f['mut_id'].tolist()),
                 ('score', 0), ('strand', intersection_f.iloc[:,5].tolist())]
        event_info_bed = pd.DataFrame.from_items(bed2)
        bedtrack_output_file = open(path1 + "/Neoskipping_track.bed", 'w')
        bedtrack_output_file.write("track name=Mutations description=\"Mutations\" color=138,0,0\n")
        event_info_bed.to_csv(bedtrack_output_file, sep="\t", index=False, header=False, mode='a')
        #Generate bedtracks of the neoskipping
        start2 = list(map(lambda x: int(x.split(";")[1])-1,new_exon_unique))
        end2 = list(map(lambda x: int(x.split(";")[2])+1,new_exon_unique))
        bed2 = [("chr", chr), ("start", start2), ("end", end2), ("id", new_exon_unique),("score", 0),("strand", strand)]
        bed_file2 = pd.DataFrame.from_items(bed2)
        bedtrack_output_file.write("track name=Neoskipping description=\"Neoskipping\" color=0,0,0\n")
        bed_file2.to_csv(bedtrack_output_file, sep="\t", index=False, header=False, mode='a')
        bedtrack_output_file.close()
        logger.info("Saved " + path1 + "/Neoskipping_track.bed")

        #Associate to each exon in the original file, the mutations
        samples = neoskipping['Neoskipping_junction'].apply(lambda x: get_associated_samples(x,exon_samples))
        aux_offset_A5 = neoskipping['Neoskipping_junction'].apply(lambda x: get_associated_samples(x,exon_offset_A5))
        aux_offset_A3 = neoskipping['Neoskipping_junction'].apply(lambda x: get_associated_samples(x,exon_offset_A3))
        neoskipping["mut_samples"] = samples
        neoskipping["mut_offset_A5"] = aux_offset_A5
        neoskipping["mut_offset_A3"] = aux_offset_A3
        #Check if any of the samples associated to the exonization has a mutation
        coincidence = neoskipping.apply(lambda x: check_coincidences(x), axis=1)
        neoskipping["mut_coincidence"] = coincidence

        neoskipping.to_csv(output_path, sep="\t", index=False, header=True)

        logger.info("Saved "+output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)