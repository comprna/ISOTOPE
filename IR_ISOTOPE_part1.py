"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
IR_ISOTOPE.py: get significat intron retention
"""

from lib.IR.extract_significant_IR import *
from lib.IR.IR_associate_gene_ids import *
from lib.IR.IR_kma_associate_gene_ids import *
from lib.IR.filter_IR import *
from lib.IR.filter_IR_CHESS import *
from lib.IR.generate_random_intronic_positions import *
import os
import csv

from argparse import ArgumentParser, RawTextHelpFormatter
import argparse

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


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

description = \
"Description: Get IR events\n\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-i", "--introns", required=True, help = "transcript expression to introns")
parser.add_argument("-b", "--bam", required=True, help = "path to STAR output")
parser.add_argument("-g", "--gtf", required=True, help = "gtf annotation")
parser.add_argument("-introns_normal", "--introns_normal", required=False, help = "transcript expression to introns on normal controls")
parser.add_argument("-introns_GTEX", "--introns_GTEX", required=False, help = "transcript expression to introns on GTEX samples")
parser.add_argument("-t", "--thres", required=False,  type=int, default=1, help="TPM threshold")
parser.add_argument("--tumor_specific", type=str2bool, nargs='?',const=True, default=False,help="Tumor specific mode")
parser.add_argument("--Rudin", type=str2bool, nargs='?',const=True, default=False,help="Rudin mode")
parser.add_argument("-o", "--output", required=True, help = "Output path")

def main(introns_path, bam_path, gtf_path, introns_Normal_path, introns_GTEX_path,
         TPM_threshold, tumor_specific, flag_Rudin, output_path):

    try:

        logger.info("Starting execution IR_ISOTOPE_part1")

        # introns_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_George_Peifer_Rudin_Yokota.txt"
        # bam_path = "/projects_rg/SCLC_cohorts/George/STAR/all_bams"
        # TPM_threshold = 1
        # tumor_specific = True
        # flag_Rudin = False
        # introns_Normal_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_Rudin_Normal.txt"
        # introns_GTEX_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.gtf"
        # gtf_protein_coding_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/SCLC/epydoor/IR"

        # 0.1. Create a gtf with only the exon information
        logger.info("Part0...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        gtf_path_exon = '{}.{}'.format(gtf_path, "exon")
        gtf = pd.read_table(gtf_path, delimiter="\t",header=None,comment="#")
        #Get only the information on the exons and on chromosomes from 1 to 22, X and Y
        gtf.columns = ['chr', 'type1', 'type2', 'start', 'end', 'dot', 'strand', 'dot2', 'rest_information']
        gtf = gtf[gtf['type2'].isin(["exon"])]
        gtf = gtf[gtf['chr'].isin(list(range(1,22)) + ["X","Y"])]
        #Add the chr suffix
        gtf['chr'] = 'chr' + gtf['chr'].astype(str)
        #Save the gtf in external file
        gtf.to_csv(gtf_path_exon,index=False,header=False,sep ='\t',quoting=csv.QUOTE_NONE)

        # 0.2. Format the intron file
        command0 = "Rscript " + dir_path + "/lib/IR/format_intron_file.R " + introns_path + " " + output_path + "/IR_formatted.tab"
        os.system(command0)

        # 1. Get the IR expressed
        logger.info("Part1...")
        extract_significant_IR(output_path + "/IR_formatted.tab", TPM_threshold, output_path + "/IR_expressed.tab")
        # extract_significant_IR(introns_path, TPM_threshold, output_path + "/IR_expressed.tab")

        # 2. Obtain the gene ids for the introns.
        logger.info("Part2...")
        # Separate between introns from kma (U2) and U12
        command1="head -n1 "+output_path + "/IR_expressed.tab > "+output_path + "/IR_kma_expressed.tab; grep kma_introns "\
                 +output_path + "/IR_expressed.tab >> "+output_path + "/IR_kma_expressed.tab"
        os.system(command1)
        command2 = "grep -v kma_introns "+output_path + "/IR_expressed.tab > "+output_path + "/IR_no_kma_expressed.tab"
        os.system(command2)
        IR_associate_gene_ids(output_path + "/IR_no_kma_expressed.tab", gtf_path, output_path + "/IR_no_kma_expressed_genes.tab")
        IR_kma_associate_gene_ids(output_path + "/IR_kma_expressed.tab", gtf_path, output_path + "/IR_kma_expressed_genes.tab")
        command3 = "cat "+output_path + "/IR_kma_expressed_genes.tab > "+output_path + "/IR_expressed_genes.tab; tail -n+2 "\
                   +output_path + "/IR_no_kma_expressed_genes.tab >> "+output_path + "/IR_expressed_genes.tab"
        os.system(command3)

        # 3. Get the IR tumor specific
        if(tumor_specific):

            if(flag_Rudin):
                #Get the significant introns for the set of normal
                extract_significant_IR(introns_Normal_path, TPM_threshold, output_path + "/IR_expressed_Normal.tab")

                #Filter by a set of Normal
                output_path_filtered = output_path + "/IR_expressed_genes_filtered.tab"
                filter_IR(output_path + "/IR_expressed_genes.tab", output_path + "/IR_expressed_Normal.tab", output_path_filtered)

                # Filter by a set of Normal (GTEX)
                output_path_filtered2 = output_path + "/IR_expressed_genes_filtered2.tab"
                filter_IR_CHESS(output_path_filtered, introns_GTEX_path, output_path_filtered2)

            else:
                # Filter by a set of Normal (GTEX)
                output_path_filtered2 = output_path + "/IR_expressed_genes_filtered2.tab"
                filter_IR_CHESS(output_path + "/IR_expressed_genes.tab", introns_GTEX_path, output_path_filtered2)


        else:
            output_path_filtered2 = output_path + "/IR_expressed_genes.tab"

        # 4. Generate random positions for each intron
        logger.info("Part4...")
        generate_random_intronic_positions(output_path_filtered2, gtf_path_exon, 100, output_path + "/random_introns.gtf",
                                           output_path + "/random_introns.bed")

        # 5. Run coverageBed on the samples in the cluster
        logger.info("Part5...")
        # dir_path = os.path.dirname(os.path.realpath(__file__))

        # 5.1. If there is any chr missing in the bed file, add an extra line with this info
        introns = pd.read_table(output_path + "/random_introns.bed",names=["chr", "start", "end", "id", "strand", "zero"])
        chr_unique = introns.chr.unique().tolist()
        chr_set = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                   "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

        for element in chr_set:
            if (element not in chr_unique):
                with open(output_path + "/random_introns.bed", "a") as file:
                    file.write(element+"\t1\t1\tExonization_0_Random_0\t+\t0\n")

        #Sort the df by chr
        introns = pd.read_table(output_path + "/random_introns.bed",names=["chr", "start", "end", "id", "strand", "zero"])
        #Add a numeric columns associated with the chromosome
        introns["chr_num"] = introns["chr"].apply(lambda x: x[3:].rstrip())
        #X tranform it to 23 and Y to 24
        introns["chr_num"] = introns["chr_num"].replace('X', 23)
        introns["chr_num"] = introns["chr_num"].replace('Y', 24)
        introns.chr_num = pd.to_numeric(introns.chr_num, errors='coerce')
        introns.start = pd.to_numeric(introns.start, errors='coerce')
        introns.end = pd.to_numeric(introns.end, errors='coerce')
        introns_sorted = introns.sort_values(by=['chr_num','start','end'])
        # remove the last column and save
        del introns_sorted['chr_num']
        introns_sorted.to_csv(output_path + "/random_introns.bed", sep="\t", index=False)

        # 5.2. Run a job per sample
        command3="for sample in $(ls "+bam_path+"/*/*.bam);do " \
                "sample_id=$(echo $sample | awk -F '/' '{print $(NF-1)}');" \
                "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed "+dir_path+"/coverageBed.sh $(echo $sample) " \
                 + output_path + "/random_introns.bed "+output_path+"/$(echo $sample_id).coverage_sorted;done"
        os.system(command3)
        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.introns,args.bam,args.gtf,args.introns_normal, args.introns_GTEX,
         args.thres,args.tumor_specific,args.Rudin, args.output)
