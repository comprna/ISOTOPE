"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

check_mutations_nearby: check if in the exonizations there are mutations nearby
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

def get_offset(line):
    '''
    Offset between the canoncial and the alternative junction
    '''
    can_junction = line[3].split("|")[0]
    start_can = int(can_junction.split(";")[1])
    end_can = int(can_junction.split(";")[2])
    alt_junction = line[3].split("|")[1]
    start_alt = int(alt_junction.split(";")[1])
    end_alt = int(alt_junction.split(";")[2])
    mut_pos = int(line[8])
    if(start_can==start_alt):
        offset = mut_pos-end_alt
    else:
        offset = mut_pos-start_alt
    return offset

def get_positions(line):
    type = line["splice_site_type"]
    strand = line["strand"]
    alt_start = int(line["Alt_Exon_id"].split(";")[1])
    alt_end = int(line["Alt_Exon_id"].split(";")[2])
    if((strand=="+" and type=="New_acceptor") or (strand=="-" and type=="New_donor")):
        start = alt_start -2
        end = alt_start
    else:
        start = alt_end +1
        end = alt_end +3
    return (start,end)

def check_mutations_nearby(exonizations_path, mutations_path, window, output_path):

    try:
        logger.info("Starting execution")

        #Load the exonizations
        exonizations = pd.read_table(exonizations_path, delimiter="\t", )
        # Create a bed file
        chr = exonizations["Canonical_Junction_id"].apply(lambda x: x.split(";")[0])
        #Add to the start and end the window shift
        positions = exonizations.apply(get_positions,axis=1)
        start = positions.apply(lambda x: int(x[0])-window)
        end = positions.apply(lambda x: int(x[1])+window)
        id = exonizations.apply(lambda x: x["Canonical_Junction_id"]+"|"+x["Alt_Junction_id"],axis=1)
        exonizations["id"] = id
        # Save this variables as bed file
        path1 = "/".join(output_path.split("/")[:-1])
        bed = [("chr", chr), ("start", start), ("end", end), ("id", id),("score", 0),("strand", exonizations["strand"])]
        bed_file = pd.DataFrame.from_items(bed)
        # bed_file['score'] = 0
        bed_file.to_csv(path1 + "/exonizations.bed", sep="\t", index=False, header=False)

        # Run interesectBed for obtaining the new exons that are not in coding regions
        logger.info("Running intersectBed...")
        command = "module load BEDTools; intersectBed -wao -a " + path1 + "/exonizations.bed -b " + mutations_path + " > " + \
                  path1 + "/intersection_mutations.bed; module unload module load BEDTools"
        os.system(command)

        # Take from the output file the exons obtained with intersectBed
        logger.info("Obtain exonizations overlapping with mutations...")
        intersection = pd.read_table(path1 + "/intersection_mutations.bed", delimiter="\t", header=None)
        # intersection.columns = ['chr', 'start', 'end', 'exon_id','strand', 'score']
        # Get the rows with a 1 in the last column
        intersection_f = intersection.loc[intersection.iloc[:,21] == 1]

        offset = intersection_f.apply(get_offset,axis=1)

        # original_start = intersection_f[3].apply(lambda x: int(x.split(";")[1]))
        # original_end = intersection_f[3].apply(lambda x: int(x.split(";")[2]))
        intersection_f["offset"] = offset
        # intersection_f["original_end"] = original_end
        # offset_A5 = intersection_f.apply(lambda x: get_offset_A5(x),axis=1)
        # offset_A3 = intersection_f.apply(lambda x: get_offset_A3(x),axis=1)
        #Extract the samples per exon that has associated mutation
        #Generate also bed tracks for the visualization of the mutations and the exonizations
        exon_samples, exon_offset = {},{}
        mut_id = intersection_f.apply(lambda x: x[10]+":"+x[11]+">"+x[12], axis=1)
        intersection_f["mut_id"] = mut_id
        for i in range(0,len(intersection_f.index)):
            exon = intersection_f.iloc[i, 3]
            mut_id_aux = intersection_f.iloc[i].loc["mut_id"]
            if(exon not in exon_samples):
                exon_samples[exon] = [mut_id_aux]
                exon_offset[exon] = [offset.iloc[i]]
            else:
                exon_samples[exon].append(mut_id_aux)
                exon_offset[exon] = [offset.iloc[i]]

        #Generate bedtracks of the mutations
        logger.info("Generating bedtracks...")
        # event_info_bed = pd.DataFrame({'chr': intersection_f.iloc[:,6].tolist(), 'start': intersection_f.iloc[:,7].tolist(),
        #                                'end': intersection_f.iloc[:,8].tolist(), 'id': intersection_f['mut_id'].tolist(),
        #                                'score': 0, 'strand': intersection_f.iloc[:,4].tolist()})
        bed2 = [('chr', intersection_f.iloc[:,6].tolist()), ('start', intersection_f.iloc[:,7].tolist()),
                ('end', intersection_f.iloc[:,8].tolist()), ('id', intersection_f['mut_id'].tolist()),
                 ('score', 0), ('strand', intersection_f.iloc[:,5].tolist())]
        event_info_bed = pd.DataFrame.from_items(bed2)
        #Get unique rows from the df
        event_info_bed2 = event_info_bed.drop_duplicates()
        bedtrack_output_file = open(path1 + "/A5_A3_junctions_track.bed", 'w')
        bedtrack_output_file.write("track name=Mutations description=\"Mutations\" color=138,0,0\n")
        event_info_bed2.to_csv(bedtrack_output_file, sep="\t", index=False, header=False, mode='a')
        #Generate bedtracks of the alternative junctions
        id = intersection_f.apply(lambda x: x.iloc[6]+";"+str(int(x.iloc[1])+window)+";"+str(int(x.iloc[2]) - window), axis=1)
        alt_junction = intersection_f.iloc[:,3].apply(lambda x: x.split("|")[1])
        start2 = intersection_f.iloc[:,1].apply(lambda x: int(x)+window)
        end2 = intersection_f.iloc[:,2].apply(lambda x: int(x)-window)
        bed2 = [("chr", intersection_f.iloc[:,6].tolist()), ("start", start2.tolist()), ("end", end2.tolist()),
                ("id", id.tolist()),("score", 0),("strand", intersection_f.iloc[:,5].tolist())]
        bed_file2 = pd.DataFrame.from_items(bed2)
        #Get unique rows from the df
        bed_file3 = bed_file2.drop_duplicates()
        bedtrack_output_file.write("track name=A5_A3_junctions description=\"A5_A3_junctions\" color=0,0,0\n")
        bed_file3.to_csv(bedtrack_output_file, sep="\t", index=False, header=False, mode='a')
        bedtrack_output_file.close()
        logger.info("Saved " + path1 + "/A5_A3_junctions_track.bed")

        #Associate to each exon in the original file, the mutations
        samples = exonizations['id'].apply(lambda x: get_associated_samples(x,exon_samples))
        aux_offset = exonizations['id'].apply(lambda x: get_associated_samples(x,exon_offset))
        # aux_offset_A3 = exonizations['New_exon'].apply(lambda x: get_associated_samples(x,exon_offset_A3))
        exonizations["mut_samples"] = samples
        exonizations["mut_offset"] = aux_offset
        # exonizations["mut_offset_A3"] = aux_offset_A3
        #Check if any of the samples associated to the exonization has a mutation
        coincidence = exonizations.apply(lambda x: check_coincidences(x), axis=1)
        #Remove the id column
        del exonizations['id']
        exonizations["mut_coincidence"] = coincidence

        exonizations.to_csv(output_path, sep="\t", index=False, header=True)

        logger.info("Saved "+output_path)
        logger.info("Done. Exiting program.")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
