"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_peptide_sequence: from our list of A5_A3 junctions:
    - Get the peptide sequence from the reference and the alternative transcript
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


def check_overlappings(df1, df2):
    '''
    Return True if the introduced exon is not overlapping any other exon
    False in other case
    '''
    s2 = int(df2["start"].tolist().pop())
    e2 = int(df2["end"].tolist().pop())
    # Check if the exon in df2 is overlapping any other exon in df1
    for i in range(0, len(df1.index)):
        s1 = int(df1.iloc[i, :]["start"])
        e1 = int(df1.iloc[i, :]["end"])
        if ((s1 < s2 and s2 < e1) or (s1 < e2 and e2 < e1)):
            return False
    return True

def check_exonization(id, exons):
    '''
    Returns True if the exonization exists, False on any other case
    '''
    can_exon = id.split("|")[0]
    alt_exon = id.split("|")[1]
    can_exon_start = can_exon.split(";")[1]
    can_exon_end = can_exon.split(";")[2]
    can_exon_strand = can_exon.split(";")[3]
    alt_exon_start = alt_exon.split(";")[1]
    alt_exon_end = alt_exon.split(";")[2]
    alt_exon_strand = alt_exon.split(";")[3]

    # 5.1 Go over all the exons checking where the canonical exon is located
    # Save the position where the A3_A5_junction is located
    for i in range(0, len(exons.index)):
        start = exons.iloc[i, 3]
        end = exons.iloc[i, 4]
        strand = exons.iloc[i, 6]
        if (strand != can_exon_strand):
            return False
        if (int(can_exon_start) == int(start) and int(can_exon_end) == int(end)):
            junction_position = i
            # Substitue the exon by the alternative exon
            line = exons.iloc[i]
            line["start"] = alt_exon_start
            line["end"] = alt_exon_end
            line["rest_information"] = "Alt_exon"
            df = pd.Series.to_frame(line).transpose()
            # If there is only 1 exon
            if (len(exons.index) == 1):
                return True
            # If is the last exon
            elif (i == len(exons.index) - 1):
                # Check if there is no overlap between the new exon and the rest
                if (check_overlappings(exons.iloc[:i], df)):
                    return True
            # If is the first exon
            elif (i == 0):
                # Check if there is no overlap between the new exon and the rest
                if (check_overlappings(exons.iloc[i + 1:], df)):
                    return True
            else:
                # Check if there is no overlap between the new exon and the rest
                if (check_overlappings(pd.concat(
                        [exons.iloc[:i], exons.iloc[i + 1:]]).reset_index(drop=True), df)):
                    return True

def get_peptide_sequence(exonizations_path, transcript_expression_path, gtf_path, codons_gtf_path, output_peptide_path,
                         output_sequence_path, output_path2, output_path3, output_path4, output_path5, mosea,
                         fast_genome, orfs_scripts, interpro, IUPred, remove_temp_files):

    try:
        logger.info("Starting execution")

        # 0. Load the exonizations with the gene associated
        exonizations_gene, exonizations_event_type = {}, {}
        with open(exonizations_path) as f:
            logger.info("Processing A5_A3_junctions file...")
            next(f)
            for line in f:
                tokens = line.rstrip().split("\t")
                # gene = tokens[3]
                # exonization = tokens[1]
                gene = tokens[1]
                can_exon = tokens[4]
                alt_exon = tokens[5]
                id = can_exon + "|" + alt_exon
                if(id not in exonizations_gene):
                    exonizations_gene[id] = gene
                else:
                    # logger.info("Repeated id  " + str(id))
                    pass
                splice_site_type = tokens[11]
                if (id not in exonizations_event_type):
                    exonizations_event_type[id] = splice_site_type
                else:
                    # logger.info("Repeated id  " + str(id))
                    pass

        # 1. Load the expression associated to each transcript
        logger.info("Load the expression associated to each transcript...")
        transcript_expression = {}
        with open(transcript_expression_path) as f:
            for line in f:
                tokens = line.rstrip().split("\t")
                transcript = tokens[0]
                tpm = tokens[1]
                # Save the values
                if (transcript not in transcript_expression):
                    transcript_expression[transcript] = tpm
                else:
                    logger.info("Repeated transcript " + transcript + " in transcript_expression")

        # 2. Get the association gene - transcript from the gtf
        gene_transcript = {}
        with open(gtf_path) as f:
            logger.info("Loading genes - transcripts...")
            for line in f:
                if (re.search("#", line)):
                    pass
                else:
                    tokens = line.rstrip().split("\t")
                    gene_id = re.sub("\"", "", tokens[8].split(";")[0].split("gene_id")[1]).strip()
                    transcript_id = re.sub("\"", "",
                                           tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                    if (gene_id not in gene_transcript):
                        gene_transcript[gene_id] = [transcript_id]
                    else:
                        if (transcript_id not in gene_transcript[gene_id]):
                            gene_transcript[gene_id].append(transcript_id)

        # 3. Get the start and end codon from each transcript
        transcript_start_codon, transcript_stop_codon = {}, {}
        with open(codons_gtf_path) as f:
            logger.info("Loading start and end codons...")
            for line in f:
                if (re.search("#", line)):
                    pass
                else:
                    tokens = line.rstrip().split("\t")
                    if (re.search("start_codon", line)):
                        start = tokens[3]
                        end = tokens[4]
                        transcript_id = re.sub("\"", "",
                                               tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                        if (transcript_id not in transcript_start_codon):
                            transcript_start_codon[transcript_id] = (start, end)
                        else:
                            # logger.info("Repeated start codon for " + str(transcript))
                            pass
                    elif (re.search("stop_codon", line)):
                        start = tokens[3]
                        end = tokens[4]
                        transcript_id = re.sub("\"", "",
                                               tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                        if (transcript_id not in transcript_stop_codon):
                            transcript_stop_codon[transcript_id] = (start, end)
                        else:
                            # logger.info("Repeated end codon for " + str(transcript))
                            pass
                    else:
                        pass

        # 4. Load the gtf as a pandas dataframe
        logger.info("Loading gtf file...")
        gtf = pd.read_table(gtf_path, delimiter="\t", header=None)
        gtf.columns = ['chr', 'type1', 'type2', 'start', 'end', 'dot', 'strand', 'dot2', 'rest_information']
        gtf["transcript_id"] = gtf["rest_information"].apply(lambda x: x.split(";")[1].split("\"")[1])

        # 5. Get the peptidic sequences of the reference and the transcript with the A5_A3 junctions replacement
        peptide_change, frame_shift, created_sequences, NMD, Stalling = {}, {}, {}, {}, {}
        index_DNA_ref, index_DNA_ex, index_AA_ref, index_AA_ex = {}, {}, {}, {}
        exonization_transcript = {}
        path1 = "/".join(output_peptide_path.split("/")[:-1])
        outFile_peptide = open(output_peptide_path, 'w')
        outFile_sequence = open(output_sequence_path, 'w')
        outFile_peptide_Interpro = open(path1 + "/A5_A3_peptide_sequence_Interpro.temp", 'w')
        outFile_IUPred = open(output_path5, 'w')
        outFile_IUPred.write("transcript\tfeatureType\tfeature_id\tstart\tend\n")
        cont1 = 0
        with open(exonizations_path) as f:
            logger.info("Processing exonizations file...")
            next(f)
            for line in f:
                tokens = line.rstrip().split("\t")
                gene = tokens[1]
                can_exon = tokens[4]
                alt_exon = tokens[5]
                id = can_exon + "|" + alt_exon
                can_exon_start = can_exon.split(";")[1]
                can_exon_end = can_exon.split(";")[2]
                can_exon_strand = can_exon.split(";")[3]
                alt_exon_start = alt_exon.split(";")[1]
                alt_exon_end = alt_exon.split(";")[2]
                alt_exon_strand = alt_exon.split(";")[3]
                event_type = exonizations_event_type[id]
                flag_exit = False
                cont1 += 1
                logger.info(str(cont1))
                # if(cont1==2):
                #     break
                # Intitialize the dictionaries
                peptide_change[id] = False
                NMD[id] = False
                # frame_shift[id] = False
                Stalling[id] = False
                created_sequences[id] = False
                # Get the transcripts associated to the gene
                if (gene in gene_transcript):
                    associated_transcripts = gene_transcript[gene]
                else:
                    logger.info("Gene " + gene + " not in gtf")
                    continue

                TPM_associated = 0
                transcript_id = "None"
                for transcript in associated_transcripts:

                    # Get the exons associated to this transcript
                    if (can_exon_strand == "+"):
                        exons_associated = (gtf.loc[gtf['transcript_id'] == transcript]).sort_values('start')
                    else:
                        exons_associated = (gtf.loc[gtf['transcript_id'] == transcript]).sort_values('start',
                                                                                                          ascending=False)

                    # Check if the neoskiipping is included on this transcript
                    if (check_exonization(id, exons_associated) and transcript in transcript_expression):
                        # Get the TPM expression and the id. We will take the transcript with the greatest expression
                        if (float(transcript_expression[transcript]) > TPM_associated):
                            TPM_associated = float(transcript_expression[transcript])
                            transcript_id = transcript

                            # If no transcript has been chosen, that means that there is no transcript in the annotation with this exonization
                if (transcript_id == "None"):
                    logger.info("A5_A3 event " + id + " not in the annotation")
                    continue

                # Associate this transcript to tje exonization
                exonization_transcript[id] = transcript_id

                # Get the exons associated to this transcript
                if (can_exon_strand == "+"):
                    exons_associated = (gtf.loc[gtf['transcript_id'] == transcript_id]).sort_values('start')
                else:
                    exons_associated = (gtf.loc[gtf['transcript_id'] == transcript_id]).sort_values('start',ascending=False)

                # 5.1 Go over all the exons checking where the canonical exon is located
                # Save the position where the A3_A5_junction is located
                for i in range(0, len(exons_associated.index)):
                    start = exons_associated.iloc[i, 3]
                    end = exons_associated.iloc[i, 4]
                    strand = exons_associated.iloc[i, 6]
                    if (strand != can_exon_strand):
                        raise Exception("Different strands associated to the same gene!!!")
                    if (int(can_exon_start) == int(start) and int(can_exon_end) == int(end)):
                        junction_position = i
                        # Substitue the exon by the alternative exon
                        line = exons_associated.iloc[i]
                        line["start"] = alt_exon_start
                        line["end"] = alt_exon_end
                        line["rest_information"] = "Alt_exon"
                        df = pd.Series.to_frame(line).transpose()
                        # If there is only 1 exon
                        if (len(exons_associated.index) == 1):
                            flag_exit = True
                            exons_associated_with_exonization = pd.concat([df]).reset_index(drop=True)
                            break
                        # If is the last exon
                        elif (i == len(exons_associated.index) - 1):
                            # Check if there is no overlap between the new exon and the rest
                            if (check_overlappings(exons_associated.iloc[:i], df)):
                                flag_exit = True
                                exons_associated_with_exonization = pd.concat(
                                    [exons_associated.iloc[:i], df]).reset_index(drop=True)
                                break
                        # If is the first exon
                        elif (i == 0):
                            # Check if there is no overlap between the new exon and the rest
                            if (check_overlappings(exons_associated.iloc[i + 1:], df)):
                                flag_exit = True
                                exons_associated_with_exonization = pd.concat(
                                    [df, exons_associated.iloc[i + 1:]]).reset_index(drop=True)
                                break
                        else:
                            # Check if there is no overlap between the new exon and the rest
                            if (check_overlappings(pd.concat(
                                    [exons_associated.iloc[:i], exons_associated.iloc[i + 1:]]).reset_index(drop=True),
                                                   df)):
                                flag_exit = True
                                exons_associated_with_exonization = pd.concat(
                                    [exons_associated.iloc[:i], df, exons_associated.iloc[i + 1:]]).reset_index(drop=True)
                                break

                if (flag_exit):
                    created_sequences[id] = True
                    # 5.2 Get the dna sequence associated with this exons and the ones without the exonization, using MosEA
                    # 5.2.1. Format the exonization exons it in a bed format
                    # remember to substract 1 to the start position
                    path1 = "/".join(output_peptide_path.split("/")[:-1])
                    exons_associated_with_exonization['start'] = exons_associated_with_exonization['start'].apply(
                        lambda x: str(int(x) - 1))
                    # Include in the id the start and end of the exon
                    id_formatted = exons_associated_with_exonization.apply(
                        lambda x: id + ":" + x['chr'] + ":" +
                                  str(x['start']) + "-" + str(x['end']), axis=1)
                    bed = [("chr", exons_associated_with_exonization['chr']),
                           ("start", exons_associated_with_exonization['start']),
                           ("end", exons_associated_with_exonization['end']), ("id", id_formatted),
                           ("strand", exons_associated_with_exonization['strand'])]
                    bed_file = pd.DataFrame.from_items(bed)
                    bed_file['score'] = 0
                    bed_file.to_csv(path1 + "/aux_exonization_A5_A3.bed", sep="\t", index=False, header=False)
                    # Format the reference transcript in a bed format
                    exons_associated['start'] = exons_associated['start'].apply(lambda x: str(int(x) - 1))
                    id_formatted = exons_associated.apply(lambda x: id + ":" + x['chr'] + ":" +
                                                                             str(x['start']) + "-" + str(x['end']),
                                                                   axis=1)
                    bed = [("chr", exons_associated['chr']), ("start", exons_associated['start']),
                           ("end", exons_associated['end']), ("id", id_formatted),
                           ("strand", exons_associated['strand'])]
                    bed_file = pd.DataFrame.from_items(bed)
                    bed_file['score'] = 0
                    bed_file.to_csv(path1 + "/aux_reference_A5_A3.bed", sep="\t", index=False, header=False)

                    # 5.2.2. Get the sequence from Mosea
                    # logger.info("Obtaining fasta exonizations sequence...")
                    command1 = "module load Python/2.7.11; module load BEDTools; python " + mosea + " getfasta --bedfile " + \
                               path1 + "/aux_exonization_A5_A3.bed --genome " + fast_genome + " --output " + path1 + \
                               "/aux_exonization_A5_A3.fa" + "; module unload Python/2.7.11"
                    # print(command1)
                    os.system(command1)

                    # logger.info("Obtaining fasta reference sequence...")
                    command2 = "module load Python/2.7.11; module load BEDTools; python " + mosea + " getfasta --bedfile " + \
                               path1 + "/aux_reference_A5_A3.bed --genome " + fast_genome + " --output " + path1 + \
                               "/aux_reference_A5_A3.fa" + "; module unload Python/2.7.11"
                    # print(command2)
                    os.system(command2)

                    # 5.3 Get the peptidic sequence
                    # logger.info("Obtaining peptidic exonizations sequence...")

                    # 5.3.1. Get the reference ORF using the information from the GTF with the start and stop codons
                    # Read the DNA reference file
                    # If the transcript is not in the gtf, there wont be change
                    if (transcript_id not in transcript_start_codon or transcript_id not in transcript_stop_codon):
                        continue
                    else:
                        # 5.3.1.1. Check if the start codon exists in the aberrant transcript. Also, obtain the similar sequence
                        cont_same_exons = 0
                        cont_start_codon = 0
                        counter = 0
                        if (can_exon_strand == "+"):
                            start_codon = transcript_start_codon[transcript_id][0]
                            stop_codon = transcript_stop_codon[transcript_id][1]
                        else:
                            start_codon = transcript_start_codon[transcript_id][1]
                            stop_codon = transcript_stop_codon[transcript_id][0]
                        sequence_similar = ""
                        flag_diff, flag_interval, flag_save_sequence, flag_exit, flag_break = False, False, False, False, False
                        with open(path1 + "/aux_exonization_A5_A3.fa") as f1, open(path1 + "/aux_reference_A5_A3.fa") as f2:
                            for x, y in zip(f1, f2):
                                if (re.search(">", x) and re.search(">", y)):
                                    if (x == y):
                                        coordinates = y.split(":")[2]
                                        start_coordinates = coordinates.split("-")[0]
                                        end_coordinates = coordinates.split("-")[1][:-4]
                                        if (int(start_coordinates) <= int(start_codon) <= int(end_coordinates)):
                                            flag_save_sequence = True
                                            flag_diff = True
                                            if (can_exon_strand == "+"):
                                                offset = (int(start_codon), int(end_coordinates))
                                                sequence_line = "y"
                                            else:
                                                offset = (int(start_coordinates), int(start_codon))
                                                sequence_line = "y"
                                        elif (flag_save_sequence):
                                            flag_diff = True
                                            offset = (int(start_coordinates), int(end_coordinates))
                                            sequence_line = "y"
                                            # counter += 1
                                    else:
                                        # Get the offset relative to the shorter exon
                                        aberrant_coordinates = x.split(":")[2]
                                        aberrant_start_coordinates = aberrant_coordinates.split("-")[0]
                                        aberrant_end_coordinates = aberrant_coordinates.split("-")[1][:-4]
                                        # aberrant_length = int(aberrant_start_coordinates )- int(aberrant_end_coordinates)
                                        aberrant_length = int(aberrant_end_coordinates) - int(aberrant_start_coordinates)
                                        ref_coordinates = y.split(":")[2]
                                        ref_start_coordinates = ref_coordinates.split("-")[0]
                                        ref_end_coordinates = ref_coordinates.split("-")[1][:-4]
                                        # ref_length = int(ref_start_coordinates )- int(ref_end_coordinates)
                                        ref_length = int(ref_end_coordinates) - int(ref_start_coordinates)
                                        if (event_type == "New_donor"):  # A5 event
                                            if (flag_save_sequence):  # the start codon was already found. Save the sequence until the discrepancy point
                                                flag_save_sequence = True
                                                flag_exit = True
                                                flag_interval = True
                                                if (aberrant_length < ref_length):
                                                    offset = (int(aberrant_end_coordinates), int(aberrant_start_coordinates))
                                                    # The sequence to be save comes from the aberrant exon
                                                    sequence_line = "x"
                                                else:
                                                    offset = (int(ref_end_coordinates), int(ref_start_coordinates))
                                                    # The sequence to be save comes from the reference exon
                                                    sequence_line = "y"
                                            else:  # the start codon has not been found yet. Check if it's in both exons
                                                if (int(aberrant_start_coordinates) <= int(start_codon) <= int(
                                                        aberrant_end_coordinates) and
                                                                int(ref_start_coordinates) <= int(start_codon) <= int(
                                                            ref_end_coordinates)):
                                                    if (aberrant_length < ref_length):
                                                        flag_save_sequence = True
                                                        flag_exit = True
                                                        flag_diff = True
                                                        # The sequence to be save comes from the aberrant exon
                                                        sequence_line = "x"
                                                        offset = (int(start_codon), int(aberrant_end_coordinates))
                                                    else:
                                                        flag_save_sequence = True
                                                        flag_exit = True
                                                        flag_diff = True
                                                        # The sequence to be save comes from the reference exon
                                                        sequence_line = "y"
                                                        offset = (int(start_codon), int(ref_end_coordinates))
                                                elif (int(ref_start_coordinates) <= int(start_codon) <= int(
                                                        ref_end_coordinates)):
                                                    # It's only in the reference. No similarity sequence. Necessary to recalculate a start_codon for the aberrant
                                                    break
                                                else:  # The start codon is after the junction. There wont be change in both ORFs
                                                    flag_break = True
                                                    break

                                        else:  # A3 event
                                            if (flag_save_sequence):  # the start codon was already found. It will be necessary to caluclate a new sc
                                                break
                                                # if (can_exon_strand == "+"):
                                                #     if (aberrant_length < ref_length):
                                                #         offset = int(aberrant_start_coordinates) - int(ref_start_coordinates)
                                                #     else:
                                                #         offset = int(ref_start_coordinates) - int(aberrant_start_coordinates)
                                                # else:
                                                #     if (aberrant_length < ref_length):
                                                #         offset = int(ref_end_coordinates) - int(aberrant_end_coordinates)
                                                #     else:
                                                #         offset = int(aberrant_end_coordinates) - int(ref_end_coordinates)
                                            else:  # the start codon has not been found yet. Check if it's in both exons
                                                if (int(aberrant_start_coordinates) <= int(start_codon) <= int(
                                                        aberrant_end_coordinates) and
                                                                int(ref_start_coordinates) <= int(start_codon) <= int(
                                                            ref_end_coordinates)):  # The start codon is after the junction. There wont be change in both ORFs
                                                    flag_break = True
                                                    break
                                                elif (int(ref_start_coordinates) <= int(start_codon) <= int(
                                                        ref_end_coordinates)):
                                                    # It's only in the reference. No similarity sequence. Necessary to recalculate a start_codon for the aberrant
                                                    break
                                                else:  # The start codon is after the junction. There wont be change in both ORFs
                                                    flag_break = True
                                                    break

                                else:
                                    if (flag_save_sequence):
                                        if (sequence_line == "y"):
                                            sequence_line2 = y
                                        else:
                                            sequence_line2 = x
                                        if (flag_diff and can_exon_strand == "+"):
                                            sequence_similar += sequence_line2.rstrip()[-(offset[1] - offset[0] + 1):]
                                            flag_diff = False
                                        elif (flag_diff and can_exon_strand == "-"):
                                            sequence_similar = sequence_line2.rstrip()[
                                                               :(offset[1] - offset[0])] + sequence_similar
                                            flag_diff = False
                                        elif (flag_interval and can_exon_strand == "+"):
                                            # sequence_similar += y.rstrip()[offset[0] : offset[1]]
                                            sequence_similar += sequence_line2.rstrip()[:(offset[0] - offset[1])]
                                            flag_interval = False
                                        elif (flag_interval and can_exon_strand == "-"):
                                            # sequence_similar = y.rstrip()[offset[0] : offset[1]] + sequence_similar
                                            sequence_similar = sequence_line2.rstrip()[
                                                               -(offset[0] - offset[1] + 1):] + sequence_similar
                                            flag_interval = False
                                        else:
                                            logger.info("Case not contemplated " + id)
                                            pass
                                    if (flag_exit):
                                        break

                        if (flag_break):
                            continue

                        # # If the start codon is after the A5_A3 junctions, the ORFs will be similar and there wont be any peptide change
                        # if (cont_start_codon > cont_same_exons):
                        #     peptide_change[id] = False
                        #     frame_shift[id] = False
                        #     NMD[id] = False
                        #     Stalling[id] = False
                        #     continue
                        else:
                            # 5.3.1.2. Get the reference sequence given by the start and stop codons from the GTF
                            cont2 = 0
                            flag_start, flag_end, flag_same_exons = False, False, True
                            sequence_total_REF = ""
                            if (can_exon_strand == "+"):
                                start_codon = transcript_start_codon[transcript_id][0]
                                stop_codon = transcript_stop_codon[transcript_id][1]
                                with open(path1 + "/aux_reference_A5_A3.fa") as f:
                                    for line in f:
                                        # If the exons are the same in the refernece and the exonization, we will store the
                                        # sequence_similar between the ref and the ex in an additional variable
                                        if (flag_same_exons and cont2 > cont_same_exons):
                                            # sequence_similar = sequence_total_REF
                                            flag_same_exons = False
                                        # If its header, pass the line
                                        if (re.search(">", line)):
                                            cont2 += 1
                                            coordinates = line.split(":")[2]
                                            start_coordinates = coordinates.split("-")[0]
                                            end_coordinates = coordinates.split("-")[1][:-4]
                                            offset1 = -1
                                            offset2 = -1
                                            pass
                                        else:
                                            sequence = line.rstrip()
                                            # Find for the start codon
                                            if (not flag_start and not flag_end and int(start_coordinates) <= int(
                                                    start_codon) <= int(
                                                    end_coordinates)):
                                                # Get the relative position in the sequence
                                                offset1 = int(start_codon) - int(start_coordinates)
                                                flag_start = True
                                            if (flag_start and not flag_end and int(start_coordinates) <= int(
                                                    stop_codon) <= int(
                                                    end_coordinates)):
                                                # Get the relative position in the sequence
                                                offset2 = int(stop_codon) - int(start_coordinates)
                                                flag_end = True
                                            # Store the sequence inside the ORF
                                            if (offset1 != -1 and offset2 != -1):
                                                sequence_total_REF = sequence_total_REF + sequence[offset1 - 1:offset2]
                                            elif (offset1 == -1 and offset2 != -1):
                                                sequence_total_REF = sequence_total_REF + sequence[:offset2]
                                            elif (offset1 != -1 and offset2 == -1):
                                                sequence_total_REF = sequence_total_REF + sequence[offset1 - 1:]
                                            elif (flag_start and offset1 == -1 and offset2 == -1):
                                                sequence_total_REF = sequence_total_REF + sequence
                                            else:
                                                pass
                                            # If both flags are True, stop the iteration
                                            if (flag_start and flag_end):
                                                break
                            else:  # can_exon_strand=="-"
                                start_codon = transcript_start_codon[transcript_id][1]
                                stop_codon = transcript_stop_codon[transcript_id][0]
                                with open(path1 + "/aux_reference_A5_A3.fa") as f:
                                    for line in f:
                                        # If the exons are the same in the refernece and the exonization, we will store the sequence
                                        # in an additional variable
                                        if (flag_same_exons and cont2 > cont_same_exons):
                                            # sequence_similar = sequence_total_REF
                                            flag_same_exons = False
                                        # If its header, pass the line
                                        if (re.search(">", line)):
                                            cont2 += 1
                                            coordinates = line.split(":")[2]
                                            start_coordinates = coordinates.split("-")[0]
                                            end_coordinates = coordinates.split("-")[1][:-4]
                                            offset1 = -1
                                            offset2 = -1
                                            pass
                                        else:
                                            sequence = line.rstrip()
                                            # Find for the start codon
                                            if (not flag_start and not flag_end and int(start_coordinates) <= int(
                                                    start_codon) <= int(
                                                    end_coordinates)):
                                                # Get the relative position in the sequence
                                                # offset1 = int(start_codon) - int(start_coordinates)
                                                offset1 = int(end_coordinates) - int(start_codon)
                                                flag_start = True
                                            if (flag_start and not flag_end and int(start_coordinates) <= int(
                                                    stop_codon) <= int(
                                                    end_coordinates)):
                                                # Get the relative position in the sequence
                                                # offset2 = int(stop_codon) - int(start_coordinates)
                                                offset2 = int(end_coordinates) - int(stop_codon)
                                                flag_end = True
                                            # Store the sequence inside the ORF
                                            if (offset1 != -1 and offset2 != -1):
                                                # sequence_total_REF = sequence_total_REF + sequence[offset1 - 1:offset2]
                                                sequence_total_REF = sequence_total_REF + sequence[-offset2:-offset1 - 1]
                                            elif (offset1 == -1 and offset2 != -1):
                                                # sequence_total_REF = sequence_total_REF + sequence[:offset2]
                                                sequence_total_REF = sequence[-offset2 - 1:] + sequence_total_REF
                                            elif (offset1 != -1 and offset2 == -1):
                                                # sequence_total_REF = sequence_total_REF + sequence[offset1 - 1:]
                                                sequence_total_REF = sequence[:-offset1] + sequence_total_REF
                                            elif (flag_start and offset1 == -1 and offset2 == -1):
                                                sequence_total_REF = sequence + sequence_total_REF
                                            else:
                                                pass
                                            # If both flags are True, stop the iteration
                                            if (flag_start and flag_end):
                                                break

                            # 5.3.2. Get the ORFs associated to the exonization DNA sequence
                            # If this is strand negative, we have to reverse the sequences before
                            sequence_total_EX = ""
                            with open(path1 + "/aux_exonization_A5_A3.fa") as f:
                                for line in f:
                                    # If its header, pass the line
                                    if (re.search(">", line)):
                                        pass
                                    else:
                                        sequence = line.rstrip()
                                        if (can_exon_strand == "+"):
                                            sequence_total_EX = sequence_total_EX + sequence
                                        else:
                                            my_seq = Seq(sequence)
                                            rev_compl_sequence = my_seq.reverse_complement()
                                            sequence_total_EX = sequence_total_EX + rev_compl_sequence

                            outFile_aux = open(path1 + "/aux_sequence_total_EX_A5_A3.fa", "w")
                            outFile_aux.write(">" + transcript_id + "_exonized|" + id + "\n")
                            outFile_aux.write(str(sequence_total_EX) + "\n")
                            outFile_aux.close()

                            # 5.3.2.1. Run extract_orfs.py for obtaining all possible ORFs in the sequence
                            # logger.info("Obtaining ORFs...")
                            command1 = "module load Python/2.7.11; python " + orfs_scripts + " " + path1 + \
                                       "/aux_sequence_total_EX_A5_A3.fa" + " 50 > " + path1 + "/aux_sequence_total_EX_ORF_A5_A3.fa" \
                                       + " ; module unload Python/2.7.11"
                            # print(command1)
                            os.system(command1)

                            # 5.3.2.2. Get the ORF of the exonizations. 2 options:
                            #   - 5.3.2.2.1: The start codon of the reference it exists also in the aberrant. Then take the
                            # one with the shortest length that starts with the similar sequence. Check if both peptides has the same ORF
                            #   - 5.3.2.2.2: The start codon doesn't exist in the exonization. Then get the smaller ORF downstream of
                            # the previous start codon
                            ORF_EX = ""
                            if (sequence_similar != ""):  # 5.3.2.2.1
                                # Check the file from the end
                                # If the gene is in reverse, get the rev_compl from the sequence_similar
                                if (can_exon_strand == "-"):
                                    my_seq = Seq(sequence_similar)
                                    sequence_similar = my_seq.reverse_complement()
                                flag_found2 = False
                                Stalling[id] = False
                                for line in reversed(list(open(path1 + "/aux_sequence_total_EX_ORF_A5_A3.fa"))):
                                    if (re.search(">", line)):
                                        pass
                                    else:
                                        # if (re.search(str(sequence_similar), line.rstrip())):
                                        if (line.rstrip().startswith(str(sequence_similar))):
                                            ORF_EX = line.rstrip()
                                            flag_found2 = True
                                            break
                                if (not flag_found2):
                                    # No stop codon
                                    Stalling[id] = True
                                    continue
                                # If there is no ORF from the previous step, go to 5.3.2.2.2
                                if (ORF_EX == ""):
                                    sequence_similar = ""
                            if (sequence_similar == ""):  # 5.3.2.2.2
                                # 5.3.2.2.2.1. Get the relative position of the start codon
                                relative_pos = 0
                                with open(path1 + "/aux_exonization_A5_A3.fa") as f:
                                    for line in f:
                                        # If it is not header, pass the line
                                        if (re.search(">", line)):
                                            coordinates = line.split(":")[2]
                                            start_coordinates = int(coordinates.split("-")[0])
                                            end_coordinates = int(coordinates.split("-")[1][:-4])
                                            exon_length = end_coordinates - start_coordinates
                                            if (can_exon_strand == "+" and event_type == "New_donor"):  # A5
                                                if (start_coordinates > int(start_codon)):
                                                    # If the start codon is in the previous exon,stop
                                                    break
                                                else:
                                                    relative_pos += exon_length
                                            elif (can_exon_strand == "-" and event_type == "New_donor"):  # A5
                                                if (end_coordinates < int(
                                                        start_codon)):  # If the start codon is in the previous exon,stop
                                                    break
                                                else:
                                                    relative_pos += exon_length
                                            elif (can_exon_strand == "+" and event_type == "New_acceptor"):  # A3
                                                if (int(start_codon) < start_coordinates):
                                                    # If the start codon is in the ref but not the ab, the newsc will be at the begining of the new exon, stop
                                                    break
                                                else:
                                                    relative_pos += exon_length
                                            elif (can_exon_strand == "-" and event_type == "New_acceptor"):  # A3
                                                if (int(start_codon) > end_coordinates):
                                                    # If the start codon is in the ref but not the ab, the newsc will be at the begining of the new exon, stop
                                                    break
                                                else:
                                                    relative_pos += exon_length
                                        else:
                                            pass

                                # 5.3.2.2.2.2. Get the first ORF downstream from the start codon
                                ORF_upstream_val = 999999999
                                ORF_length = 999999999
                                flag_seq = False
                                Stalling[id] = False
                                with open(path1 + "/aux_sequence_total_EX_ORF_A5_A3.fa") as f:
                                    for line in f:
                                        # If its header, take the ORF_interval. If it's downstream from the start_codon, save it
                                        if (re.search(">", line)):
                                            lengths = line.split(":")[1]
                                            if (int(lengths.split("-")[0]) > relative_pos):
                                                if (ORF_upstream_val > int(lengths.split("-")[0])):
                                                    ORF_upstream_val = int(lengths.split("-")[0])
                                                    ORF_length = int(lengths.split("-")[1]) - int(lengths.split("-")[0])
                                                    flag_seq = True
                                                # If it's in the exact same position, take the smaller
                                                elif (ORF_upstream_val == int(lengths.split("-")[0]) and
                                                              ORF_length > int(lengths.split("-")[1]) - int(
                                                              lengths.split("-")[0])):
                                                    ORF_length = int(lengths.split("-")[1]) - int(lengths.split("-")[0])
                                                    flag_seq = True
                                        else:
                                            if (flag_seq):
                                                ORF_EX = line.rstrip()
                                                flag_seq = False

                                # If we haven't found any ORF, there will be stalling. Break
                                if (ORF_upstream_val == 999999999):
                                    Stalling[id] = True
                                    continue

                                # 5.3.2.2.2.3. We need to recover also the exact position of the start codon of the ORF_EX
                                distance = ORF_upstream_val
                                with open(path1 + "/aux_exonization_A5_A3.fa") as f:
                                    for line in f:
                                        # If it is not header, pass the line
                                        if (re.search(">", line)):
                                            coordinates = line.split(":")[2]
                                            start_coordinates = int(coordinates.split("-")[0])
                                            end_coordinates = int(coordinates.split("-")[1][:-4])
                                            aux_length = end_coordinates - start_coordinates
                                            if (int(distance) > aux_length):
                                                distance -= aux_length
                                            else:
                                                if (can_exon_strand == "+"):
                                                    start_codon = start_coordinates + distance
                                                else:
                                                    start_codon = end_coordinates - distance
                                                break
                                        else:
                                            pass

                            # 5.3.3. Get the translation from the ORFs (reference and exonization)
                            ORF_EX_f = ORF_EX.replace("T", "U")
                            messenger_rna = Seq(ORF_EX_f, IUPAC.unambiguous_rna)
                            peptide_exonizations = messenger_rna.translate()

                            # If the gene is in reverse, get the rev_compl from the sequence_total_REF
                            if (can_exon_strand == "-"):
                                my_seq = Seq(sequence_total_REF)
                                sequence_total_REF = my_seq.reverse_complement()
                            ORF_REF_f = str(sequence_total_REF).replace("T", "U")
                            messenger_rna = Seq(ORF_REF_f, IUPAC.unambiguous_rna)
                            peptide_reference = messenger_rna.translate()

                            # 5.4. Save both DNA and peptidic sequences to the output
                            outFile_peptide.write(">" + transcript_id + "\n")
                            outFile_peptide.write(str(peptide_reference) + "\n")
                            outFile_peptide.write(">" + transcript_id + "_exonized|" + id + "\n")
                            outFile_peptide.write(str(peptide_exonizations) + "\n")

                            # Save also the sequence fasta output. First the reference.
                            outFile_sequence.write(">" + transcript_id + "\n")
                            outFile_sequence.write(str(sequence_total_REF) + "\n")
                            outFile_sequence.write(">" + transcript_id + "_exonized|" + id + "\n")
                            outFile_sequence.write(str(ORF_EX) + "\n")

                            # Save the sequences in a separate structure for outputing in another file
                            index_DNA_ref[id] = str(sequence_total_REF)
                            index_DNA_ex[id] = str(ORF_EX)
                            index_AA_ref[id] = str(peptide_reference)
                            index_AA_ex[id] = str(peptide_exonizations)

                            # Output only the peptide sequence of the exonizations (for the Intrepro prediction)
                            outFile_peptide_Interpro.write(">" + id + "\n")
                            outFile_peptide_Interpro.write(str(peptide_exonizations).replace("*", "") + "\n")

                            # Output only the peptide sequence of the exonizations (for the IUPred prediction)
                            outFile_peptide_IUPred = open(path1 + "/A5_A3_peptide_sequence_IUPred.temp", 'w')
                            outFile_peptide_IUPred.write(">" + id + "\n")
                            outFile_peptide_IUPred.write(str(peptide_exonizations).replace("*", "") + "\n")
                            outFile_peptide_IUPred.close()

                            # Run IUPred for obtaining the disordered regions
                            command4 = "module load Python; python " + IUPred + "/iupred2a.py " + path1 + \
                                       "/A5_A3_peptide_sequence_IUPred.temp long > " + path1 + "/A5_A3_peptide_sequence_IUPred.temp.out; " \
                                                                                               "module unload Python;"
                            os.system(command4)

                            # Process the output of IUPred
                            flag_new_interval = True
                            length = 0
                            with open(path1 + "/A5_A3_peptide_sequence_IUPred.temp.out") as f:
                                for line in f:
                                    if (re.search("#", line)):
                                        continue
                                    else:
                                        tokens = line.rstrip().split("\t")
                                        prediction_value = float(tokens[2])
                                        AA_position = int(tokens[0])
                                        if (prediction_value > 0.5 and flag_new_interval):
                                            flag_new_interval = False
                                            start = AA_position
                                            length = 1
                                        elif (not flag_new_interval):
                                            if (prediction_value > 0.5):
                                                length += 1
                                                continue
                                            else:
                                                flag_new_interval = True
                                                end = AA_position - 1
                                                # Close the interval and save it if the length of the interval is greater than 5
                                                if (length >= 5):
                                                    outFile_IUPred.write(id + "\tIDR\tDisordered_region\t" + str(start) +
                                                                         "\t" + str(end) + "\n")

                            # 5.5. If there is a peptide change, check if the exonized sequence will go to NMD
                            peptide_change[id] = (not peptide_reference == peptide_exonizations)
                            if (peptide_reference == peptide_exonizations):
                                NMD[id] = False
                            else:
                                # Count the number of the exons in the file. Also check if the start codon is in on any
                                # of the exons
                                if (can_exon_strand == "+"):
                                    start_codon = transcript_start_codon[transcript_id][0]
                                else:
                                    start_codon = transcript_start_codon[transcript_id][1]
                                n_exons = 0
                                flag_start_codon = False
                                with open(path1 + "/aux_exonization_A5_A3.fa") as f:
                                    for line in f:
                                        # Count the lines with a header
                                        if (re.search(">", line)):
                                            n_exons += 1
                                            coordinates = line.split(":")[2]
                                            start_coordinates = coordinates.split("-")[0]
                                            end_coordinates = coordinates.split("-")[1][:-4]
                                            if (int(start_coordinates) <= int(start_codon) <= int(end_coordinates)):
                                                flag_start_codon = True

                                # Read the DNA exonization file, checking where the stop codon is falling
                                cont3 = 0
                                flag_start, flag_end = False, False
                                ORF_counter_positions = len(ORF_EX)
                                if (can_exon_strand == "+"):
                                    start_codon = transcript_start_codon[transcript_id][0]
                                else:
                                    start_codon = transcript_start_codon[transcript_id][1]
                                with open(path1 + "/aux_exonization_A5_A3.fa") as f:
                                    for line in f:
                                        # If its header, pass the line
                                        if (re.search(">", line)):
                                            cont3 += 1
                                            coordinates = line.split(":")[2]
                                            start_coordinates = coordinates.split("-")[0]
                                            end_coordinates = coordinates.split("-")[1][:-4]
                                            offset1 = -1
                                            offset2 = -1
                                            pass
                                        else:
                                            sequence = line.rstrip()
                                            # Find for the start codon
                                            if (not flag_start and not flag_end and int(start_coordinates) <= int(
                                                    start_codon) <= int(
                                                end_coordinates)):
                                                # Get the number of bp of the ORF that are in this exon
                                                if (can_exon_strand == "+"):
                                                    offset1 = int(end_coordinates) - int(start_codon)
                                                else:
                                                    offset1 = int(start_codon) - int(start_coordinates)
                                                flag_start = True
                                            # Substract the number of positions left in the ORF
                                            if (flag_start and offset1 != -1):
                                                ORF_counter_positions -= offset1
                                            elif (flag_start and offset1 == -1):
                                                ORF_counter_positions -= len(sequence)
                                            else:
                                                pass
                                            # If ORF_counter_positions < 0, stop iterating
                                            # Check where the stop codon is falling
                                            if (ORF_counter_positions - 1 <= 0):
                                                # If there are more than 50 bp remaining and is not the last exon, the transcript goes to NMD
                                                if (abs(ORF_counter_positions - 1) >= 50 and n_exons != cont3):
                                                    NMD[id] = True
                                                else:
                                                    NMD[id] = False
                                                break

                                    if (ORF_counter_positions - 1 > 0):
                                        logger.info("ORF_counter_positions over 0!! for " + id)
                                        # raise Exception("ORF_counter_positions over 0!! for "+id)

                else:
                    pass

        outFile_peptide.close()
        outFile_sequence.close()
        outFile_peptide_Interpro.close()
        outFile_IUPred.close()

        # 6. Add the columns to the initial list
        logger.info("Processing output file...")
        outFile2 = open(output_path2, "w")
        outFile3 = open(output_path3, "w")
        index = 0
        with open(exonizations_path) as f:
            outFile2.write("Index\t" + next(
                f).rstrip() + "\tTranscript_id\tTranscript_TPM\tCreated_sequences?\tPeptide_change\tNMD\tStalling\n")
            outFile3.write("Index\tORF_reference\tORF_changed\tPeptide_reference\tPeptide_changed\n")
            for line in f:
                tokens = line.rstrip().split("\t")
                gene = tokens[1]
                can_exon = tokens[4]
                alt_exon = tokens[5]
                id = can_exon + "|" + alt_exon
                if(id in exonization_transcript):
                    transcript = exonization_transcript[id]
                    if(transcript in transcript_expression):
                        tpm = float(transcript_expression[transcript])
                    else:
                        tpm = 0
                else:
                    transcript = ""
                    tpm = 0
                if (id in index_DNA_ref and id and index_DNA_ex and id in index_AA_ref and id in index_AA_ex):
                    outFile3.write(
                        str(index) + "\t" + index_DNA_ref[id] + "\t" + index_DNA_ex[id] + "\t" + index_AA_ref[id] +
                        "\t" + index_AA_ex[id] + "\n")
                else:
                    outFile3.write(str(index) + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\n")
                if (id in peptide_change and id in NMD and id in created_sequences and id in Stalling):
                    outFile2.write(str(index) + "\t" + line.rstrip() + "\t" + transcript + "\t" + str(tpm) + "\t" + str(
                        created_sequences[id]) + "\t" +
                                   str(peptide_change[id]) + "\t" + str(NMD[id]) + "\t" + str(Stalling[id]) + "\n")
                    index += 1
                else:
                    raise Exception("A5_A3 junction " + id + " not in dictionaries")

        outFile2.close()
        outFile3.close()

        # 7. Run Interpro
        logger.info("Run Interpro...")
        command3 = interpro + " -i " + path1 + "/A5_A3_peptide_sequence_Interpro.temp -f tsv -o " + output_path4
        os.system(command3)

        logger.info("Saved " + output_sequence_path)
        logger.info("Saved " + output_peptide_path)
        logger.info("Saved " + output_path2)
        logger.info("Saved " + output_path3)
        logger.info("Saved " + output_path4)
        logger.info("Saved " + output_path5)
        logger.info("Done. Exiting program.")

        # Remove temporary files
        if (remove_temp_files):
            os.remove(path1 + "/aux_exonization_A5_A3.bed")
            os.remove(path1 + "/aux_reference_A5_A3.bed")
            os.remove(path1 + "/aux_exonization_A5_A3.fa")
            os.remove(path1 + "/aux_reference_A5_A3.fa")
            os.remove(path1 + "/aux_sequence_total_EX_A5_A3.fa")
            os.remove(path1 + "/aux_sequence_total_EX_ORF_A5_A3.fa")
            os.remove(path1 + "/A5_A3_peptide_sequence_Interpro.temp")
            os.remove(path1 + "/A5_A3_peptide_sequence_IUPred.temp")
            os.remove(path1 + "/A5_A3_peptide_sequence_IUPred.temp.out")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    exonizations_path = sys.argv[1]
    transcript_expression_path = sys.argv[2]
    gtf_path = sys.argv[3]
    codons_gtf_path = sys.argv[4]
    output_peptide_path = sys.argv[5]
    output_sequence_path = sys.argv[6]
    output_path2 = sys.argv[7]
    output_path3 = sys.argv[8]
    output_path4 = sys.argv[9]
    output_path5 = sys.argv[10]
    mosea = sys.argv[11]
    fast_genome = sys.argv[12]
    orfs_scripts = sys.argv[13]
    interpro = sys.argv[14]
    IUPred = sys.argv[15]
    remove_temp_files = sys.argv[16]
    get_peptide_sequence(exonizations_path, transcript_expression_path, gtf_path, codons_gtf_path, output_peptide_path,
                         output_sequence_path, output_path2, output_path3, output_path4, output_path5, mosea,
                         fast_genome, orfs_scripts, interpro, IUPred, remove_temp_files)
