"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
IR_ISOTOPE.py: get significat intron retention
"""

import csv
from lib.IR.extract_significant_IR import *
from lib.IR.IR_associate_gene_ids import *
from lib.IR.filter_IR import *
from lib.IR.filter_IR_CHESS import *
from lib.IR.generate_random_intronic_positions import *
from lib.IR.get_coverageBed import *
from lib.IR.get_coverageBed_adapter import *
from lib.IR.get_peptide_sequence_RI import *
from lib.IR.select_fasta_candidates import *
from lib.IR.run_netMHC_classI_slurm_part1 import *
from lib.IR.run_netMHC_classI_slurm_part2 import *
from lib.IR.run_netMHCpan_classI_slurm_part1 import *
from lib.IR.run_netMHCpan_classI_slurm_part2 import *

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
parser.add_argument("-trans", "--transcript", required=True, help="transcript expression file")
parser.add_argument("-g", "--gtf", required=True, help="gtf annotation")
parser.add_argument("-genome", "--genome", required=True, help="Genome annotation")
parser.add_argument("-HLAclass", "--HLAclass", required=True, help="HLA genotype of the samples")
parser.add_argument("-HLAtypes", "--HLAtypes", required=True, help="HLA alelles recognized by NetMHC")
parser.add_argument("-HLAtypespan", "--HLAtypespan", required=True, help="HLA alelles recognized by NetMHCpan")
parser.add_argument("-netMHC", "--netMHC", required=True, help="netMHC path")
parser.add_argument("-netMHCpan", "--netMHCpan", required=True, help="netMHCpan path")
parser.add_argument("-t", "--thres", required=False, type=int, default=1,
                    help="Minimum expression to consider an intron")
parser.add_argument("-mosea", "--mosea", required=True, help="MoSEA path")
parser.add_argument("-mxfinder", "--mxfinder", required=True, help="MxFinder path")
parser.add_argument("-o", "--output", required=True, help="Output path")
parser.add_argument("--username", required=True, help="Cluster user name")
parser.add_argument("--tumor_specific", type=str2bool, nargs='?', const=True, default=False,
                    help="Tumor specific mode")
parser.add_argument("--temp", type=str2bool, nargs='?', const=True, default=False, help="Remove temp files")
parser.add_argument("-c", "--cluster", type=str2bool, nargs='?',const=True, default=False,help="Run in parallel on a cluster")

def main(transcript_expression_path, gtf_path, genome_path, HLAclass_path, HLAtypes_path,
         HLAtypes_pan_path, netMHC_path, netMHC_pan_path, threshold,
         mosea_path, mxfinder_path, output_path, tumor_specific, remove_temp_files, name_user, cluster):

    try:

        logger.info("Starting execution IR_ISOTOPE_part2")

        # transcript_expression_path = "/projects_rg/SCLC_cohorts/George/tables/iso_tpm_George_Peifer_Rudin_Yokota.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        # mosea = "/genomics/users/juanluis/Software/MoSEA-master/mosea.py"
        # fasta_genome = "/genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa"
        # orfs_scripts = "/genomics/users/juanluis/comprna/MxFinder/extract_orfs.py"
        # interpro = "/soft/EB_repo/bio/sequence/programs/noarch/interproscan/5.33-72.0/interproscan.sh"
        # IUPred = "/projects_rg/SCLC_cohorts/soft/IUPred2A"
        # HLAclass_path = "/projects_rg/SCLC_cohorts/tables/PHLAT_summary_ClassI_all_samples.out"
        # HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        # HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        # netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        # netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        # remove_temp_files = True
        # tumor_specific = True
        # name_user = "juanluis"
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/SCLC/epydoor/IR"
        # # ONLY FOR MARVIN
        # #python2 = "Python/2.7.14-foss-2017b"
        # # ONLY FOR HYDRA
        # python2 = "Python/2.7.11"

        # 0.1. Create a gtf with only the exon information
        dir_path = os.path.dirname(os.path.realpath(__file__))
        gtf_path_exon = '{}.{}'.format(gtf_path, "exon")
        gtf = pd.read_table(gtf_path, delimiter="\t",header=None,comment="#")
        #Get only the information on the exons and on chromosomes from 1 to 22, X and Y
        gtf.columns = ['chr', 'type1', 'type2', 'start', 'end', 'dot', 'strand', 'dot2', 'rest_information']
        gtf = gtf[gtf['type2'].isin(["exon"])]
        gtf = gtf[gtf['chr'].isin(list(range(1,23)) + ["X","Y"])]
        #Add the chr suffix
        gtf['chr'] = 'chr' + gtf['chr'].astype(str)
        #Save the gtf in external file
        gtf.to_csv(gtf_path_exon,index=False,header=False,sep ='\t',quoting=csv.QUOTE_NONE)

        # 6. Create the folder, if it doesn't exists
        logger.info("Part6...")
        if not os.path.exists(output_path + "/coverageBed"):
            os.makedirs(output_path + "/coverageBed")
        # Move all the coverage.sorted files to the created directory
        command1="mv "+output_path+"/*coverage_sorted "+output_path + "/coverageBed/"
        os.system(command1)

        # 7.1. Get the coverage for each exonization
        logger.info("Part7.1...")
        if(tumor_specific):
            output_path_filtered2 = output_path + "/IR_expressed_genes_filtered2.tab"
        else:
            output_path_filtered2 = output_path + "/IR_expressed_genes.tab"

        get_coverageBed_adapter(output_path_filtered2, output_path + "/random_introns.bed",output_path + "/coverageBed", output_path, name_user, cluster)

        # 7.2. Assemble all pieces into one single file
        logger.info("Part7.2...")
        command2="awk 'FNR==1 && NR!=1{next;}{print}' "+output_path+"/get_coverageBed_*.tab > "+output_path+"/IR_coverage.tab"
        os.system(command2)

        # 7.3. Get the introns with a significant p_value
        logger.info("Part7.3...")
        command3="head -n1 "+output_path+"/IR_coverage.tab > "+output_path+"/IR_significant_introns.tab; " \
                   "awk '{ if ($7 <= 0.05 && $6 > 0) print }' "+output_path+"/IR_coverage.tab >> "+output_path+"/IR_significant_introns.tab"
        os.system(command3)

        # 8. Get the peptide sequence associated
        logger.info("Part8...")
        get_peptide_sequence(output_path + "/IR_significant_introns.tab", transcript_expression_path, gtf_path,
                             output_path + "/IR_peptide_sequence.fa", output_path + "/IR_fasta_sequence.fa",
                             output_path + "/IR_ORF.tab", output_path + "/IR_ORF_sequences.tab", mosea_path,
                             genome_path, mxfinder_path, remove_temp_files)

        # 9. Filter the significant results
        logger.info("Part9...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        command4="Rscript "+dir_path+"/lib/IR/filter_results.R "+output_path + "/IR_ORF.tab"+" "+ \
                 output_path + "/IR_ORF_filtered.tab " + str(threshold) + " "+ output_path + "/IR_ORF_filtered_peptide_change.tab"
        os.system(command4)

        # 10. Select the fasta candidates for being run to the epitope analysis
        logger.info("Part10...")
        #Create the folder, if it doesn't exists
        if not os.path.exists(output_path + "/IR_fasta_files"):
            os.makedirs(output_path + "/IR_fasta_files")
        select_fasta_candidates(output_path + "/IR_ORF_filtered_peptide_change.tab", output_path + "/IR_peptide_sequence.fa", output_path + "/IR_peptide_sequence_filtered.fa", output_path + "/IR_fasta_files")

        #11. Run netMHC-4.0_part1
        logger.info("Part11...")
        if not os.path.exists(output_path + "/IR_NetMHC-4.0_files"):
            os.makedirs(output_path + "/IR_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHC-4.0_files", output_path + "/IR_NetMHC-4.0_neoantigens_type_gained.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_gained_all.tab", output_path + "/IR_NetMHC-4.0_neoantigens_type_lost.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_lost_all.tab", output_path + "/IR_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path, cluster)

        #12. Run netMHCpan-4.0_part1
        logger.info("Part12...")
        if not os.path.exists(output_path + "/IR_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/IR_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHCpan-4.0_files", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_gained.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_gained_all.tab", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_lost.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_lost_all.tab", output_path + "/IR_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_pan_path, cluster)

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.transcript,args.gtf,args.genome,args.HLAclass,args.HLAtypes,args.HLAtypespan,
         args.netMHC,args.netMHCpan,args.thres,args.mosea,args.mxfinder,args.output,args.tumor_specific,
         args.temp,args.username,args.cluster)
    # main("/home/trincadojl/Projects/SCLC/Smart/data/iso_tpm.txt",
    #      "/home/trincadojl/Projects/SCLC/Smart/annotation/Homo_sapiens.GRCh37.75.gtf",
    #      "/media/trincadojl/data/Projects/annotation/hg19.fa",
    #      "/home/trincadojl/Projects/SCLC/Smart/data/PHLAT_summary_ClassI.out",
    #      "/home/trincadojl/Projects/SCLC/Smart/data/NetMHC-4.0_HLA_types_accepted.tab",
    #      "/home/trincadojl/Projects/SCLC/Smart/data/NetMHCpan-4.0_HLA_types_accepted.tab",
    #      "/home/trincadojl/Software/netMHC-4.0/netMHC",
    #      "/home/trincadojl/Software/netMHCpan-4.0/netMHCpan",
    #      1,
    #      "/home/trincadojl/Software/MoSEA",
    #      "/home/trincadojl/Software/MxFinder",
    #      "/home/trincadojl/Projects/SCLC/Smart/test_ISOTOPE/IR",
    #      False,
    #      False,
    #      "juanluis",
    #      False)
