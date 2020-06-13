"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
exonizations_ISOTOPE.py: get significat exonizations
"""

import os
import pandas as pd

from lib.Exonization.extract_exonized_junctions import *
from lib.Exonization.get_reads_exonizations import *
from lib.Exonization.overlap_with_repeats import *
from lib.Exonization.get_significant_exonizations import *
from lib.Exonization.generate_random_intronic_positions import *
from lib.Exonization.get_coverageBed import *
from lib.Exonization.check_mutations_nearby import *

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


def main():
    try:

        logger.info("Starting execution")

        readcounts_path = "/projects_rg/SCLC_cohorts/Smart/STAR/readCounts.tab"
        bam_path = "/projects_rg/SCLC_cohorts/Smart/STAR"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        max_length = 500
        threshold = 5
        n_randomizations = 100
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/Smart/epydoor/exonizations"

        # 1. Identify the junctions that could generate an exonization
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_exonized_junctions.tab"
        extract_exonized_junctions(readcounts_path, gtf_path, max_length, output_path_aux)

        # 2. Given the list with the possible exonizations, get the reads associate to each of them
        logger.info("Part2...")
        output_path_aux2 = output_path+"/new_exonized_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)

        # 3. find the overlap between the nex exonizations and repeatitions (RepeatMasker)
        logger.info("Part3...")
        output_path_aux3 = output_path + "/new_exonized_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the exonizations with the reads counts,get those that are over a threshold
        logger.info("Part4...")
        output_path_aux4 = output_path + "/exonizations_by_sample.tab"
        get_significant_exonizations(output_path_aux3, threshold, output_path_aux4)

        # 5. generate a number of random position by exonization
        logger.info("Part5...")
        output_path_aux5 = output_path + "/random_exonizations.gtf"
        output_path_aux6 = output_path + "/random_exonizations.bed"
        generate_random_intronic_positions(output_path_aux4, gtf_path, n_randomizations, output_path_aux5, output_path_aux6)

        # 6. Run coverageBed on the samples in the cluster

        # 6.1. If there is any chr missing in the bed file, add an extra line with this info
        introns = pd.read_table(output_path_aux6,
                                names=["chr", "start", "end", "id", "strand", "zero"])
        chr_unique = introns.chr.unique().tolist()
        chr_set = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                   "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                   "chr22", "chrX", "chrY"]

        for element in chr_set:
            if (element not in chr_unique):
                with open(output_path_aux6, "a") as file:
                    file.write(element + "\t1\t1\tExonization_0_Random_0\t+\t0\n")

        # Sort the df by chr
        introns = pd.read_table(output_path_aux6,names=["chr", "start", "end", "id", "strand", "zero"])
        # Add a numeric columns associated with the chromosome
        introns["chr_num"] = introns["chr"].apply(lambda x: x[3:].rstrip())
        # X tranform it to 23 and Y to 24
        introns["chr_num"] = introns["chr_num"].replace('X', 23)
        introns["chr_num"] = introns["chr_num"].replace('Y', 24)
        introns.chr_num = pd.to_numeric(introns.chr_num, errors='coerce')
        introns.start = pd.to_numeric(introns.start, errors='coerce')
        introns.end = pd.to_numeric(introns.end, errors='coerce')
        introns_sorted = introns.sort_values(by=['chr_num', 'start', 'end'])
        # remove the last column and save
        del introns_sorted['chr_num']
        introns_sorted.to_csv(output_path_aux6, sep="\t", index=False)

        # 6.2. Run a job per sample
        logger.info("Part6...")
        command1="for sample in $(ls "+bam_path+"/*/*.bam);do " \
                "sample_id=$(echo $sample | awk -F '/' '{print $(NF-1)}');" \
                "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed "+dir_path+"/coverageBed.sh $(echo $sample) " \
                 +output_path_aux6+" "+output_path+"/$(echo $sample_id).coverage_sorted;done"
        os.system(command1)
        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()
