"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_peptide_sequence: from our list of neoskipping:
    - Get the peptide sequence from the refernce and the alternative transcript
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


def check_neoskipping(neoskipping, exons):
    '''
    Returns True if the neoskipping exists, False on any other case
    '''
    neoskipping_start = neoskipping.split(";")[1]
    neoskipping_end = neoskipping.split(";")[2]
    neoskipping_strand = neoskipping.split(";")[3]

    flag_start, flag_end = False, False
    pos_start, pos_end = 0, 0

    # 5.1. Go over all the exons checking where the neoskipping is located
    start_prev = exons.iloc[0, 3]
    end_prev = exons.iloc[0, 4]
    strand_prev = exons.iloc[0, 6]
    if (len(exons.index) != 1):
        for i in range(1, len(exons.index)):
            start = exons.iloc[i, 3]
            end = exons.iloc[i, 4]
            strand = exons.iloc[i, 6]
            if (strand != neoskipping_strand):
                raise False
            if (strand == "+"):
                if (not flag_start and int(end_prev) == int(neoskipping_start) + 1):
                    pos_start = i
                    flag_start = True
                if (int(start) == int(neoskipping_end)):
                    pos_end = i
                    flag_end = True
                if (flag_start and flag_end):
                    if (pos_start == pos_end):
                        return False
                    else:
                        return True
            else:  # (strand=="-")
                if (int(start_prev) == int(neoskipping_end)):
                    pos_start = i
                    flag_start = True
                if (int(end) == int(neoskipping_start) + 1):
                    pos_end = i
                    flag_end = True
                if (flag_start and flag_end):
                    if (pos_start == pos_end):
                        return False
                    else:
                        return True
            start_prev = start
            end_prev = end
            strand_prev = strand

    else:  # Only 1 exon in the df, no neoskipping possible
        return False

def get_expression(sample_id,transcript_id,transcript_expression):
    if (transcript_id in transcript_expression[sample_id]):
        return transcript_expression[sample_id][transcript_id]
    else:
        return -1

def get_peptide_sequence(neoskipping_path, transcript_expression_path, gtf_path, codons_gtf_path, output_peptide_path,
                         output_sequence_path, output_path2, output_path3, output_path4, output_path5, mosea,
                         fast_genome, orfs_scripts, interpro, IUPred, remove_temp_files, python2):

    try:
        logger.info("Starting execution")

        # Create a dict per cell line
        logger.info("Load the expression associated to each transcript...")
        transcript_expression = {}
        with open(transcript_expression_path) as f:
            header = next(f).rstrip().split("\t")
            # Initialize the list of dictionaries
            for i in range(0, len(header)):
                transcript_expression[header[i]] = {}
                # transcript_expression.append({})
            for line in f:
                tokens = line.rstrip().split("\t")
                transcript = tokens[0]
                tpm = tokens[1:]
                for i in range(0, len(tpm)):
                    if (transcript not in transcript_expression[header[i]]):
                        transcript_expression[header[i]][transcript] = float(tpm[i])
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
                    transcript_id = re.sub("\"", "", tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                    if(gene_id not in gene_transcript):
                        gene_transcript[gene_id] = [transcript_id]
                    else:
                        if(transcript_id not in gene_transcript[gene_id]):
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
                        transcript_id = re.sub("\"","",tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                        if(transcript_id not in transcript_start_codon):
                            transcript_start_codon[transcript_id] = (start,end)
                        else:
                            # logger.info("Repeated start codon for " + str(transcript))
                            pass
                    elif (re.search("stop_codon", line)):
                        start = tokens[3]
                        end = tokens[4]
                        transcript_id = re.sub("\"","",tokens[8].split(";")[1].split("transcript_id")[1]).strip()
                        if(transcript_id not in transcript_stop_codon):
                            transcript_stop_codon[transcript_id] = (start,end)
                        else:
                            # logger.info("Repeated end codon for " + str(transcript))
                            pass
                    else:
                        pass

        # 4. Load the gtf as a pandas dataframe
        logger.info("Loading gtf file...")
        gtf = pd.read_table(gtf_path, delimiter="\t",header=None)
        gtf.columns = ['chr', 'type1', 'type2', 'start', 'end', 'dot', 'strand', 'dot2', 'rest_information']
        gtf["transcript_id"] = gtf["rest_information"].apply(lambda x: x.split(";")[1].split("\"")[1])

        # 5. Get the peptidic sequences of the reference and the transcript with the neoskipping
        peptide_change, frame_shift, NMD, Stalling = {}, {}, {}, {}
        index_DNA_ref, index_DNA_ex, index_AA_ref, index_AA_ex = {}, {}, {}, {}
        neoskipping_transcript = {}
        path1 = "/".join(output_peptide_path.split("/")[:-1])
        outFile_peptide = open(output_peptide_path, 'w')
        outFile_sequence = open(output_sequence_path, 'w')
        outFile_peptide_Interpro = open(path1 + "/Neoskippings_peptide_sequence_Interpro.temp", 'w')
        outFile_IUPred = open(output_path5, 'w')
        outFile_IUPred.write("transcript\tfeatureType\tfeature_id\tstart\tend\n")
        cont1 = 0
        with open(neoskipping_path) as f:
            logger.info("Processing neoskipping file...")
            header = next(f).rstrip().split("\t")
            Gene_id_pos = header.index("Gene_id")
            Sample_id_pos = header.index("Sample_id")
            Neoskipping_junction_pos = header.index("Neoskipping_junction")
            for line in f:
                tokens = line.rstrip().split("\t")
                gene = tokens[Gene_id_pos]
                neoskipping = tokens[Neoskipping_junction_pos]
                sample_id = tokens[Sample_id_pos].rstrip()
                cont1+=1
                logger.info(str(cont1))
                # if(cont1==213):
                #     print("guiuoh")
                gene_neoskipping = gene + ";" + neoskipping
                # Intitialize the dictionaries
                peptide_change[gene_neoskipping] = False
                NMD[gene_neoskipping] = False
                frame_shift[gene_neoskipping] = False
                Stalling[gene_neoskipping] = False
                # Get the transcripts associated to the gene
                if(gene in gene_transcript):
                    associated_transcripts = gene_transcript[gene]
                else:
                    logger.info("Gene "+gene+" not in gtf")
                    continue

                # Check in which transcripts the neoskipping is included
                neoskipping_start = neoskipping.split(";")[1]
                neoskipping_end = neoskipping.split(";")[2]
                neoskipping_strand = neoskipping.split(";")[3]

                TPM_associated = 0
                transcript_id = "None"
                for transcript in associated_transcripts:

                    #Get the exons associated to this transcript
                    if(neoskipping_strand=="+"):
                        exons_associated = (gtf.loc[gtf['transcript_id'] == transcript]).sort_values('start')
                    else:
                        exons_associated = (gtf.loc[gtf['transcript_id'] == transcript]).sort_values('start',ascending=False)

                    #Check if the neoskiipping is included on this transcript
                    if(check_neoskipping(neoskipping,exons_associated)):
                        TPM = get_expression(sample_id,transcript,transcript_expression)
                        if (TPM != -1):
                            # Get the TPM expression and the id. We will take the transcript with the greatest expression
                            if (TPM > TPM_associated):
                                TPM_associated = TPM
                                transcript_id = transcript



                #If no transcript has been chosen, that means that there is no transcript in the annotation with this neoskipping
                if(transcript_id == "None"):
                    logger.info("Neoskipping "+gene_neoskipping+" not in the annotation")
                    continue

                #Associate this transcript to tje neoskipping
                neoskipping_transcript[gene_neoskipping] = transcript_id

                # Get the exons associated to this transcript
                if (neoskipping_strand == "+"):
                    exons_associated = (gtf.loc[gtf['transcript_id'] == transcript_id]).sort_values('start')
                else:
                    exons_associated = (gtf.loc[gtf['transcript_id'] == transcript_id]).sort_values('start',ascending=False)

                # 5.1. Go over all the exons checking where the neoskipping is located
                flag_start, flag_end = False, False
                pos_start, pos_end = 0, 0
                start_prev = exons_associated.iloc[0, 3]
                end_prev = exons_associated.iloc[0, 4]
                strand_prev = exons_associated.iloc[0, 6]
                for i in range(1,len(exons_associated.index)):
                    start = exons_associated.iloc[i,3]
                    end = exons_associated.iloc[i,4]
                    strand = exons_associated.iloc[i, 6]
                    if(strand!=neoskipping_strand):
                        raise Exception("Different strands associated to the same gene!!!")
                    if(strand=="+"):
                        if(not flag_start and int(end_prev) == int(neoskipping_start)+1):
                            pos_start = i
                            flag_start = True
                        if(int(start)==int(neoskipping_end)):
                            pos_end = i
                            flag_end = True
                        if(flag_start and flag_end):
                            exons_associated_with_neoskipping = pd.concat([exons_associated.iloc[:pos_start],
                                                                           exons_associated.iloc[pos_end:]]).reset_index(drop=True)
                            break
                    else:  # (strand=="-")
                        # logger.info("Negative strand!!")
                        if (int(start_prev) == int(neoskipping_end)):
                            pos_start = i
                            flag_start = True
                        if (int(end) == int(neoskipping_start)+1):
                            pos_end = i
                            flag_end = True
                        if (flag_start and flag_end):
                            exons_associated_with_neoskipping = pd.concat(
                                [exons_associated.iloc[:pos_start], exons_associated.iloc[pos_end:]]).reset_index(
                                drop=True)
                            break
                    start_prev = start
                    end_prev = end

                # If both flags are true, means that we've found exons coinciding with the neoskipping junction.from
                # If dont, the main transcript dont have matching exons with the neoskipping and dont generate peptide sequence
                if(flag_start and flag_end):
                    # 5.2 Get the dna sequence associated with this exons and the ones without the neoskipping, using MosEA
                    # 5.2.1. Format the neoskipping exons it in a bed format
                    # remember to substract 1 to the start position
                    exons_associated_with_neoskipping['start'] = exons_associated_with_neoskipping['start'].apply(lambda x: str(int(x) - 1))
                    # Include in the id the start and end of the exon
                    id_formatted = exons_associated_with_neoskipping.apply(lambda x: neoskipping + ":" + x['chr'] + ":" +
                                                                            str(x['start']) + "-" + str(x['end']), axis=1)
                    bed = [("chr", exons_associated_with_neoskipping['chr']), ("start", exons_associated_with_neoskipping['start']),
                    ("end", exons_associated_with_neoskipping['end']), ("id", id_formatted),
                    ("strand", exons_associated_with_neoskipping['strand'])]
                    bed_file = pd.DataFrame.from_items(bed)
                    bed_file['score'] = 0
                    bed_file.to_csv(path1 + "/aux_neoskipping_Exoniz.bed", sep="\t", index=False, header=False)
                    # Format the reference transcript in a bed format
                    exons_associated['start'] = exons_associated['start'].apply(lambda x: str(int(x) - 1))
                    id_formatted = exons_associated.apply(lambda x: neoskipping + ":" + x['chr'] + ":" +
                                                                             str(x['start']) + "-" + str(x['end']),
                                                                   axis=1)
                    bed = [("chr", exons_associated['chr']), ("start", exons_associated['start']),
                    ("end", exons_associated['end']), ("id", id_formatted),
                    ("strand", exons_associated['strand'])]
                    bed_file = pd.DataFrame.from_items(bed)
                    bed_file['score'] = 0
                    bed_file.to_csv(path1 + "/aux_reference_Exoniz.bed", sep="\t", index=False, header=False)

                    # 5.2.2. Get the sequence from Mosea
                    # logger.info("Obtaining fasta neoskipping sequence...")
                    command1 = "module load " + python2 + " ; module load BEDTools; python " + mosea + " getfasta --bedfile " + \
                               path1 + "/aux_neoskipping_Exoniz.bed --genome " + fast_genome + " --output " + path1 + \
                               "/aux_neoskipping_Exoniz.fa" + "; module unload " + python2
                    # print(command1)
                    os.system(command1)

                    # logger.info("Obtaining fasta reference sequence...")
                    command2 = "module load " + python2 + " ; module load BEDTools; python " + mosea + " getfasta --bedfile " + \
                               path1 + "/aux_reference_Exoniz.bed --genome " + fast_genome + " --output " + path1 + \
                               "/aux_reference_Exoniz.fa" + "; module unload " + python2
                    # print(command2)
                    os.system(command2)

                    # 5.3 Get the peptidic sequence
                    # logger.info("Obtaining peptidic neoskipping sequence...")

                    # 5.3.1. Get the reference ORF using the information from the GTF with the start and stop codons
                    # Read the DNA reference file
                    # If the transcript is not in the gtf, there wont be no change
                    if (transcript_id not in transcript_start_codon or transcript_id not in transcript_stop_codon):
                        logger.info("Transcript " + transcript_id + " don't have codons in gtf")
                        continue
                    else:
                        # 5.3.1.1. Check how many initial exons are similar. We will do this in order to obtain the part of the ORF
                        # reference that should be in common with the ORF neoskipping
                        cont_same_exons = 0
                        counter = 0
                        if (neoskipping_strand == "+"):
                            start_codon = transcript_start_codon[transcript_id][0]
                        else:
                            start_codon = transcript_start_codon[transcript_id][1]
                        with open(path1 + "/aux_neoskipping_Exoniz.fa") as f1, open(path1 + "/aux_reference_Exoniz.fa") as f2:
                            for x, y in zip(f1, f2):
                                if (re.search(">", x)):
                                    # coordinates = x.split(":")[2]
                                    # Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                    coordinates = x.split(":")[-1]
                                    start_coordinates = coordinates.split("-")[0]
                                    end_coordinates = coordinates.split("-")[1][:-4]
                                    if (int(start_coordinates) <= int(start_codon) <= int(end_coordinates)):
                                        cont_start_codon = counter
                                    counter += 1
                                if (re.search(">", x) and re.search(">", y) and x == y):
                                    continue
                                # Save the sequence if they are the same exon
                                elif (not re.search(">", x) and not re.search(">", y) and x == y):
                                    cont_same_exons += 1
                                # If the exons are different, stop the loop
                                else:
                                    break

                    # If the start codon is after the neoskipping, the ORFs will be similar and there wont be any peptide change
                    # This will happen also if the stop codon is falling before the neoskipping
                    if (neoskipping_strand == "+"):
                        start_codon = transcript_start_codon[transcript_id][0]
                        stop_codon = transcript_stop_codon[transcript_id][1]
                    else:
                        start_codon = transcript_start_codon[transcript_id][1]
                        stop_codon = transcript_stop_codon[transcript_id][0]
                    if(strand == "+" and neoskipping_end < start_codon or strand == "-" and neoskipping_start > start_codon or
                       strand == "+" and neoskipping_start > stop_codon or strand == "-" and neoskipping_end < stop_codon):
                        continue
                    else:
                        # 5.3.1.2. Get the reference sequence given by the start and stop codons from the GTF
                        cont2 = 0
                        flag_start, flag_end, flag_same_exons = False, False, True
                        sequence_total_REF, sequence_similar = "", ""
                        if (neoskipping_strand == "+"):
                            with open(path1 + "/aux_reference_Exoniz.fa") as f:
                                for line in f:
                                    #If the exons are the same in the refernece and the neoskipping, we will store the sequence
                                    # in an additional variable
                                    if(flag_same_exons and cont2>cont_same_exons):
                                        sequence_similar = sequence_total_REF
                                        flag_same_exons = False
                                    # If its header, pass the line
                                    if (re.search(">", line)):
                                        cont2 += 1
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
                                        start_coordinates = coordinates.split("-")[0]
                                        end_coordinates = coordinates.split("-")[1][:-4]
                                        offset1 = -1
                                        offset2 = -1
                                        pass
                                    else:
                                        sequence = line.rstrip()
                                        #Find for the start codon
                                        if(not flag_start and not flag_end and int(start_coordinates)<=int(start_codon)<=int(end_coordinates)):
                                            # Get the relative position in the sequence
                                            offset1=int(start_codon)-int(start_coordinates)
                                            flag_start = True
                                        # Find for the stop codon
                                        if(flag_start and not flag_end and int(start_coordinates)<=int(stop_codon)<=int(end_coordinates)):
                                            # Get the relative position in the sequence
                                            offset2=int(stop_codon)-int(start_coordinates)
                                            flag_end = True
                                        #Store the sequence inside the ORF
                                        if(offset1!=-1 and offset2!=-1):
                                            sequence_total_REF = sequence_total_REF + sequence[offset1-1:offset2]
                                        elif(offset1==-1 and offset2!=-1):
                                            sequence_total_REF = sequence_total_REF + sequence[:offset2]
                                        elif(offset1!=-1 and offset2==-1):
                                            if(offset1!=0):
                                                sequence_total_REF = sequence_total_REF + sequence[offset1-1:]
                                            else:
                                                sequence_total_REF = sequence_total_REF + sequence[offset1-1:]
                                        elif(flag_start and offset1==-1 and offset2==-1):
                                            sequence_total_REF = sequence_total_REF + sequence
                                        else:
                                            pass
                                        #If both flags are True, stop the iteration
                                        if(flag_start and flag_end):
                                            break

                        else:       #neoskipping_strand == "-"
                            with open(path1 + "/aux_reference_Exoniz.fa") as f:
                                for line in f:
                                    # If the exons are the same in the refernece and the neoskipping, we will store the sequence
                                    # in an additional variable
                                    if (flag_same_exons and cont2 > cont_same_exons):
                                        sequence_similar = sequence_total_REF
                                        flag_same_exons = False
                                    # If its header, pass the line
                                    if (re.search(">", line)):
                                        cont2 += 1
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
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
                                        # Find for the stop codon
                                        if (flag_start and not flag_end and int(start_coordinates) <= int(stop_codon) <= int(
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
                                            if(offset1!=0):
                                                sequence_total_REF = sequence[:-offset1] + sequence_total_REF
                                            else:
                                                sequence_total_REF = sequence + sequence_total_REF
                                        elif (flag_start and offset1 == -1 and offset2 == -1):
                                            sequence_total_REF = sequence + sequence_total_REF
                                        else:
                                            pass
                                        # If both flags are True, stop the iteration
                                        if (flag_start and flag_end):
                                            break

                        # 5.3.2. Get the ORFs associated to the neoskipping DNA sequence
                        # If this is strand negative, we have to reverse the sequences before
                        sequence_total_EX = ""
                        with open(path1 + "/aux_neoskipping_Exoniz.fa") as f:
                            for line in f:
                                # If its header, pass the line
                                if (re.search(">", line)):
                                    pass
                                else:
                                    sequence = line.rstrip()
                                    if (neoskipping_strand == "+"):
                                        sequence_total_EX = sequence_total_EX + sequence
                                    else:
                                        my_seq = Seq(sequence)
                                        rev_compl_sequence = my_seq.reverse_complement()
                                        sequence_total_EX = sequence_total_EX + rev_compl_sequence

                        outFile_aux = open(path1 + "/aux_sequence_total_EX_Exoniz.fa","w")
                        outFile_aux.write(">"+transcript_id+"_exonized|"+neoskipping+"\n")
                        outFile_aux.write(str(sequence_total_EX)+"\n")
                        outFile_aux.close()

                        # 5.3.2.1. Run extract_orfs.py for obtaining all possible ORFs in the sequence
                        # logger.info("Obtaining ORFs...")
                        command1 = "module load " + python2 + " ; python " + orfs_scripts + " " + path1 + \
                                   "/aux_sequence_total_EX_Exoniz.fa" + " 50 > " + path1 + "/aux_sequence_total_EX_ORF_Exoniz.fa" \
                                   + " ; module unload " + python2
                        # print(command1)
                        os.system(command1)

                        # 5.3.2.2. Get the ORF of the exonizations. 2 options:
                        #   - 5.3.2.2.1: The start codon of the reference it exists also in the exonization. Then take the
                        # one with the shortest length that includes the sequence_similar. Check if both peptides has the same ORF
                        #   - 5.3.2.2.2: The start codon doesn't exist in the exonization. Then get the smaller ORF more upstream
                        ORF_EX = ""
                        if (sequence_similar != ""):  # 5.3.2.2.1
                            # Check the file from the end
                            # If the gene is in reverse, get the rev_compl from the sequence_similar
                            if (neoskipping_strand == "-"):
                                my_seq = Seq(sequence_similar)
                                sequence_similar = my_seq.reverse_complement()

                            for line in reversed(list(open(path1 + "/aux_sequence_total_EX_ORF_Exoniz.fa"))):
                                if (re.search(">", line)):
                                    pass
                                else:
                                    # if (re.search(str(sequence_similar), line.rstrip())):
                                    if (line.rstrip().startswith(str(sequence_similar))):
                                        ORF_EX = line.rstrip()
                                        break
                            # If there is no ORF from the previous step, go to 5.3.2.2.2
                            if (ORF_EX == ""):
                                sequence_similar = ""
                        if (sequence_similar == ""):  # 5.3.2.2.2
                            # 5.3.2.2.2.1. Get the relative position of the start codon
                            relative_pos = 0
                            with open(path1 + "/aux_neoskipping_Exoniz.fa") as f:
                                for line in f:
                                    # If it is not header, pass the line
                                    if (re.search(">", line)):
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
                                        start_coordinates = int(coordinates.split("-")[0])
                                        end_coordinates = int(coordinates.split("-")[1][:-4])
                                        exon_length = end_coordinates - start_coordinates
                                        #If we have skipped the start codon, stop iterating
                                        if (neoskipping_strand == "+"):
                                            if(start_coordinates > int(start_codon)):
                                                break
                                            else:
                                                relative_pos += exon_length
                                        else:   #strand "-"
                                            if (end_coordinates < int(start_codon)):
                                                break
                                            else:
                                                relative_pos += exon_length
                                    else:
                                        pass

                            # 5.3.2.2.2.2. Get the first ORF downstream from the start codon
                            ORF_upstream_val = 999999999
                            ORF_length = 999999999
                            flag_seq = False
                            Stalling[gene_neoskipping] = False
                            with open(path1 + "/aux_sequence_total_EX_ORF_Exoniz.fa") as f:
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
                                Stalling[gene_neoskipping] = True
                                continue

                            # 5.3.2.2.2.3. We need to recover also the exact position of the start codon of the ORF_EX
                            distance = ORF_upstream_val
                            with open(path1 + "/aux_neoskipping_Exoniz.fa") as f:
                                for line in f:
                                    # If it is not header, pass the line
                                    if (re.search(">", line)):
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
                                        start_coordinates = int(coordinates.split("-")[0])
                                        end_coordinates = int(coordinates.split("-")[1][:-4])
                                        aux_length = end_coordinates - start_coordinates
                                        if (int(distance) > aux_length):
                                            distance -= aux_length
                                        else:
                                            if (neoskipping_strand == "+"):
                                                start_codon = start_coordinates + distance + 1
                                            else:
                                                start_codon = end_coordinates - distance + 1
                                            break
                                    else:
                                        pass

                        # 5.3.3. Get the translation from the ORFs (reference and neoskipping)
                        ORF_EX_f = str(ORF_EX).replace("T", "U")
                        messenger_rna = Seq(ORF_EX_f, IUPAC.unambiguous_rna)
                        peptide_neoskipping = messenger_rna.translate()

                        # If the gene is in reverse, get the rev_compl from the sequence_total_REF
                        if (neoskipping_strand == "-"):
                            my_seq = Seq(sequence_total_REF)
                            sequence_total_REF = my_seq.reverse_complement()
                        ORF_REF_f = str(sequence_total_REF).replace("T", "U")
                        messenger_rna = Seq(ORF_REF_f, IUPAC.unambiguous_rna)
                        peptide_reference = messenger_rna.translate()

                        # 5.4. Save both DNA and peptidic sequences to the output
                        outFile_peptide.write(">"+transcript_id+"\n")
                        outFile_peptide.write(str(peptide_reference)+"\n")
                        outFile_peptide.write(">"+transcript_id+"_exonized|"+neoskipping+"\n")
                        outFile_peptide.write(str(peptide_neoskipping)+"\n")

                        outFile_sequence.write(">"+transcript_id+"\n")
                        outFile_sequence.write(str(sequence_total_REF)+"\n")
                        outFile_sequence.write(">"+transcript_id+"_exonized|"+neoskipping+"\n")
                        outFile_sequence.write(str(ORF_EX)+"\n")

                        # Save the sequences in a separate structure for outputing in another file
                        index_DNA_ref[gene_neoskipping] = str(sequence_total_REF)
                        index_DNA_ex[gene_neoskipping] = str(ORF_EX)
                        index_AA_ref[gene_neoskipping] = str(peptide_reference)
                        index_AA_ex[gene_neoskipping] = str(peptide_neoskipping)

                        # Output only the peptide sequence of the neoskipping (for the Intrepro prediction)
                        outFile_peptide_Interpro.write(">"+neoskipping+"\n")
                        outFile_peptide_Interpro.write(str(peptide_neoskipping).replace("*","")+"\n")

                        # Output only the peptide sequence of the neoskipping (for the IUPred prediction)
                        outFile_peptide_IUPred = open(path1 + "/Neoskippings_peptide_sequence_IUPred.temp", 'w')
                        outFile_peptide_IUPred.write(">"+neoskipping+"\n")
                        outFile_peptide_IUPred.write(str(peptide_neoskipping).replace("*","")+"\n")
                        outFile_peptide_IUPred.close()

                        # Run IUPred for obtaining the disordered regions
                        command4 = "module load "+python2+" ; python " + IUPred + "/iupred2a.py " + path1 + \
                                   "/Neoskippings_peptide_sequence_IUPred.temp long > " + path1 + \
                                   "/Neoskippings_peptide_sequence_IUPred.temp.out; " \
                                                                                "module unload "+python2+";"
                        os.system(command4)

                        # Process the output of IUPred
                        flag_new_interval = True
                        length = 0
                        with open(path1 + "/Neoskippings_peptide_sequence_IUPred.temp.out") as f:
                            for line in f:
                                if(re.search("#",line)):
                                    continue
                                else:
                                    tokens = line.rstrip().split("\t")
                                    prediction_value = float(tokens[2])
                                    AA_position = int(tokens[0])
                                    if(prediction_value>0.5 and flag_new_interval):
                                        flag_new_interval = False
                                        start = AA_position
                                        length = 1
                                    elif(not flag_new_interval):
                                        if(prediction_value>0.5):
                                            length += 1
                                            continue
                                        else:
                                            flag_new_interval = True
                                            end = AA_position-1
                                            #Close the interval and save it if the length of the interval is greater than 5
                                            if(length>=5):
                                                outFile_IUPred.write(neoskipping+"\tIDR\tDisordered_region\t"+str(start)+
                                                                     "\t"+str(end)+"\n")

                        # 5.5. If there is a peptide change, check if the exonized sequence will go to NMD
                        # if(gene_neoskipping not in peptide_change or gene_neoskipping not in NMD):
                        #     peptide_change[gene_neoskipping] = (not peptide_reference==peptide_neoskipping)
                        # else:
                        #     logger.info("Repeated neoskipping "+gene_neoskipping)
                        #     continue
                        peptide_change[gene_neoskipping] = (not peptide_reference == peptide_neoskipping)

                        if(peptide_reference==peptide_neoskipping):
                            NMD[gene_neoskipping] = False
                        else:
                            # Count the number of the exons in the file. Also check if the start codon is in on any
                            # of the exons
                            n_exons = 0
                            flag_start_codon = False
                            with open(path1 + "/aux_neoskipping_Exoniz.fa") as f:
                                for line in f:
                                    # Count the lines with a header
                                    if (re.search(">", line)):
                                        n_exons += 1
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
                                        start_coordinates = coordinates.split("-")[0]
                                        end_coordinates = coordinates.split("-")[1][:-4]
                                        if (int(start_coordinates) <= int(start_codon) <= int(end_coordinates)):
                                            flag_start_codon = True

                            # Read the DNA neoskipping file, checking where the stop codon is falling
                            cont3 = 0
                            flag_start, flag_end = False, False
                            ORF_counter_positions = len(ORF_EX)
                            with open(path1 + "/aux_neoskipping_Exoniz.fa") as f:
                                for line in f:
                                    # If its header, pass the line
                                    if (re.search(">", line)):
                                        cont3 += 1
                                        # coordinates = line.split(":")[2]
                                        #Take the last coordinates. These are the ones that MoSEA have used for obtaining the sequence
                                        coordinates = line.split(":")[-1]
                                        start_coordinates = coordinates.split("-")[0]
                                        end_coordinates = coordinates.split("-")[1][:-4]
                                        offset1 = -1
                                        offset2 = -1
                                        pass
                                    else:
                                        sequence = line.rstrip()
                                        # Find for the start codon
                                        if (not flag_start and not flag_end and int(start_coordinates) <= int(start_codon) <= int(
                                                end_coordinates)):
                                            # Get the number of bp of the ORF that are in this exon
                                            if (neoskipping_strand == "+"):
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
                                        #If ORF_counter_positions < 0, stop iterating
                                        # Check where the stop codon is falling
                                        if(ORF_counter_positions-1<=0):
                                            # If there are more than 50 bp remaining and is not the last exon, the transcript goes to NMD
                                            if (abs(ORF_counter_positions - 1) >= 50 and n_exons != cont3):
                                                NMD[gene_neoskipping] = True
                                            else:
                                                NMD[gene_neoskipping] = False
                                            break

                                if (ORF_counter_positions - 1 > 0):
                                    logger.info("ORF_counter_positions over 0!! for " + gene_neoskipping)
                                    # raise Exception("ORF_counter_positions over 0!! for "+neoskipping)
                                else:
                                    pass

                        # 5.6. Check if the sequence of the exon is not multiple of 3 (and could potentially produce a frame shift)
                        if (gene_neoskipping not in frame_shift):
                            if ((int(neoskipping_end) - (int(neoskipping_start)-1)) % 3 != 0):
                                frame_shift[gene_neoskipping] = True
                            else:
                                frame_shift[gene_neoskipping] = False
                        else:
                            pass

                else:
                    logger.info("Dont found exons matching the neoskipping")
                    continue

        outFile_peptide.close()
        outFile_sequence.close()
        outFile_peptide_Interpro.close()
        outFile_IUPred.close()

        # 6. Add the columns to the initial list
        logger.info("Processing output file...")
        logger.info("Size of dictionaries: Peptide change"+str(len(peptide_change))+" NMD: "+str(len(NMD))
                    + " frame_shift: " + str(len(frame_shift)) +" Stalling: "+str(len(Stalling)))
        outFile2 = open(output_path2,"w")
        outFile3 = open(output_path3,"w")
        index = 0
        with open(neoskipping_path) as f:
            header = next(f).rstrip()
            outFile2.write("Index\t"+header+"\tTranscript_id\tTranscript_TPM\tFrame_shift\tPeptide_change\tNMD\tStalling\n")
            outFile3.write("Index\tORF_reference\tORF_changed\tPeptide_reference\tPeptide_changed\n")
            Gene_id_pos = header.split("\t").index("Gene_id")
            Neoskipping_junction_pos = header.split("\t").index("Neoskipping_junction")
            Gene_name_pos = header.split("\t").index("external_gene_name")
            Sample_id_pos = header.split("\t").index("Sample_id")
            cont1 = 0
            for line in f:
                tokens = line.rstrip().split("\t")
                neoskipping = tokens[Neoskipping_junction_pos]
                gene = tokens[Gene_id_pos]
                gene_neoskipping = gene + ";" + neoskipping
                sample_id = tokens[Sample_id_pos]
                # cont1+=1
                # logger.info(str(cont1))
                # if(cont1==20):
                #     break
                try:
                    gene_name = tokens[Gene_name_pos]
                except IndexError:
                    gene_name = ""
                if(gene_neoskipping in neoskipping_transcript):
                    transcript = neoskipping_transcript[gene_neoskipping]
                    # tpm = get_expression(sample_id, transcript, CA46_transcript_expression, HL_60_transcript_expression, THP_1_transcript_expression)
                    tpm = get_expression(sample_id, transcript, transcript_expression)
                else:
                    transcript = ""
                    tpm = 0
                if(gene_neoskipping in index_DNA_ref and gene_neoskipping and index_DNA_ex and gene_neoskipping in index_AA_ref and gene_neoskipping in index_AA_ex):
                    outFile3.write(str(index)+"\t"+index_DNA_ref[gene_neoskipping]+"\t"+index_DNA_ex[gene_neoskipping]+"\t"+index_AA_ref[gene_neoskipping]+
                                   "\t"+index_AA_ex[gene_neoskipping]+"\n")
                else:
                    outFile3.write(str(index) + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\n")
                    # logger.info("Neoskipping " + gene_neoskipping + " not in dictionaries3")
                if(gene_neoskipping in peptide_change and gene_neoskipping in frame_shift and gene_neoskipping in NMD and gene_neoskipping in Stalling):
                    if(gene_name!=""):
                        outFile2.write(str(index)+"\t"+line.rstrip()+"\t"+transcript+"\t"+str(tpm)+"\t"+str(frame_shift[gene_neoskipping])
                                       +"\t"+str(peptide_change[gene_neoskipping])+"\t"+str(NMD[gene_neoskipping])+"\t"+str(Stalling[gene_neoskipping])+"\n")
                    else:
                        #If gene_name is empty add an extra space, for respecting the columns
                        outFile2.write(str(index) + "\t" + line.rstrip() + "\t" + " " + "\t" + transcript + "\t" + str(tpm) + "\t" +
                                       str(frame_shift[gene_neoskipping]) + "\t" + str(peptide_change[gene_neoskipping]) +
                                       "\t" + str(NMD[gene_neoskipping]) + "\t" + str(Stalling[gene_neoskipping]) + "\n")
                    index += 1
                else:
                    outFile2.write(str(index) + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\n")
                    logger.info("Neoskipping " + gene_neoskipping + " not in dictionaries2")

        outFile2.close()
        outFile3.close()

        # 7. Run Interpro
        logger.info("Run Interpro...")
        command3 = interpro + " -i " + path1 + "/Neoskippings_peptide_sequence_Interpro.temp -f tsv -o " + output_path4
        os.system(command3)

        logger.info("Saved "+output_sequence_path)
        logger.info("Saved "+output_peptide_path)
        logger.info("Saved "+output_path2)
        logger.info("Saved "+output_path3)
        logger.info("Saved "+output_path4)
        logger.info("Saved "+output_path5)
        logger.info("Done. Exiting program.")

        #Remove temporary files
        if(remove_temp_files):
            os.remove(path1 + "/aux_neoskipping_Exoniz.bed")
            os.remove(path1 + "/aux_reference_Exoniz.bed")
            os.remove(path1 + "/aux_neoskipping_Exoniz.fa")
            os.remove(path1 + "/aux_reference_Exoniz.fa")
            os.remove(path1 + "/aux_sequence_total_EX_Exoniz.fa")
            os.remove(path1 + "/aux_sequence_total_EX_ORF_Exoniz.fa")
            os.remove(path1 + "/Neoskippings_peptide_sequence_Interpro.temp")
            os.remove(path1 + "/Neoskippings_peptide_sequence_IUPred.temp")
            os.remove(path1 + "/Neoskippings_peptide_sequence_IUPred.temp.out")

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
