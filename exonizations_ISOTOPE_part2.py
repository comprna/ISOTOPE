"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
exonizations_ISOTOPE.py: get significat exonizations
"""

import os

from lib.Exonization.extract_exonized_junctions import *
from lib.Exonization.get_reads_exonizations import *
from lib.Exonization.overlap_with_repeats import *
from lib.Exonization.get_significant_exonizations import *
from lib.Exonization.generate_random_intronic_positions import *
from lib.Exonization.get_coverageBed import *
from lib.Exonization.get_coverageBed_adapter import *
from lib.Exonization.check_mutations_nearby import *
from lib.Exonization.select_fasta_candidates import *
from lib.Exonization.filter_exonizations import *
from lib.Exonization.filter_exonizations_CHESS import *
from lib.Exonization.get_peptide_sequence import *
from lib.Exonization.run_netMHC_classI_slurm_part1 import *
from lib.Exonization.run_netMHCpan_classI_slurm_part1 import *

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

description = \
"Description: Get exonization events\n\n"

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-r", "--reads", required=True, help = "reads mapped to junctions")
parser.add_argument("-g", "--gtf", required=True, help = "gtf annotation")
parser.add_argument("-genome", "--genome", required=True, help = "Genome annotation")
parser.add_argument("-trans", "--transcript", required=True, help = "transcript expression file")
parser.add_argument("-HLAclass", "--HLAclass", required=True, help = "HLA genotype of the samples")
parser.add_argument("-HLAtypes", "--HLAtypes", required=True, help = "HLA alelles recognized by NetMHC")
parser.add_argument("-HLAtypespan", "--HLAtypespan", required=True, help = "HLA alelles recognized by NetMHCpan")
parser.add_argument("-netMHC", "--netMHC", required=True, help = "netMHC path")
parser.add_argument("-netMHCpan", "--netMHCpan", required=True, help = "netMHCpan path")
parser.add_argument("-mosea", "--mosea", required=True, help = "MoSEA path")
parser.add_argument("-orfs", "--orfs", required=True, help = "MxFinder path")
parser.add_argument("-o", "--output", required=True, help = "Output path")
parser.add_argument("--username", required=True, help = "Cluster user name")
parser.add_argument("-rep", "--repeats", required=True, help = "Regions of the genome with repeats from maskerDB",default=None)
parser.add_argument("-t", "--thres", required=False,  type=int, default=5, help="Minimum number of reads mapping the event")
parser.add_argument("-m", "--max", required=False,  type=int, default=500)
parser.add_argument("-rand", "--rand", required=False,  type=int, default=100, help="Number of rounds for calculating significance of each event")
parser.add_argument("--tumor_specific", type=str2bool, nargs='?',const=True, default=False,help="Tumor specific mode")
parser.add_argument("--Rudin", type=str2bool, nargs='?',const=True, default=False,help="Rudin mode")
parser.add_argument("--remove_temp_files", type=str2bool, nargs='?',const=True, default=True,help="Remove temp files")
parser.add_argument("-mut","--mutations", required=False, help = "Mutations path")
parser.add_argument("--chess", required=False, help = "CHESS SE path")

def main(readcounts_path, gtf_path, genome_path, transcript_expression_path, HLAclass_path, HLAtypes_path, HLAtypes_pan_path, netMHC_path, netMHC_pan_path,
         mosea_path, orfs_scripts, output_path, repeats_path, threshold, tumor_specific, mutations_path, CHESS_SE_path, flag_Rudin, remove_temp_files, name_user):
    try:

        logger.info("Starting execution exonizations_ISOTOPE_part2")

        # tumor_specific = True
        # readcounts_path = "/projects_rg/SCLC_cohorts/Smart/STAR/readCounts.tab"
        # transcript_expression_path = "/projects_rg/SCLC_cohorts/Smart/Salmon/iso_tpm.txt"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        # mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        # repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        # CHESS_SE_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe"
        # mosea = "/genomics/users/juanluis/Software/MoSEA-master/mosea.py"
        # genome_path = "/genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa"
        # orfs_scripts = "/genomics/users/juanluis/comprna/MxFinder/extract_orfs.py"
        interpro = None
        IUPred = None
        # HLAclass_path = "/projects_rg/SCLC_cohorts/Smart/PHLAT/PHLAT_summary_ClassI.out"
        # HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        # HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        # netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        # netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        # remove_temp_files = True
        # flag_Rudin = False
        # threshold = 10
        # name_user = "juanluis"
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/Smart/epydoor/exonizations"

        # ONLY FOR MARVIN
        #python2 = "Python/2.7.14-foss-2017b"
        # ONLY FOR HYDRA
        python2 = "Python/2.7.11"

        # 6. Create the folder, if it doesn't exists
        logger.info("Part6...")
        if not os.path.exists(output_path + "/coverageBed"):
            os.makedirs(output_path + "/coverageBed")
        # Move all the coverage.sorted files to the created directory
        command1 = "mv " + output_path + "/*coverage_sorted " + output_path + "/coverageBed/"
        os.system(command1)

        # 7.1. Get the coverage for each exonization
        logger.info("Part7...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        get_coverageBed_adapter(output_path + "/exonizations_by_sample.tab", output_path + "/random_exonizations.bed",
                        output_path + "/coverageBed", output_path, name_user)

        # 7.2. Assemble all pieces into one single file
        command2 = "awk 'FNR==1 && NR!=1{next;}{print}' " + output_path + "/get_coverageBed_*.tab > " + output_path + "/exonizations_by_sample_coverage.tab"
        os.system(command2)

        # 8. Check if in the exonizations there are mutations nearby
        logger.info("Part8...")
        check_mutations_nearby(output_path + "/exonizations_by_sample_coverage.tab", mutations_path, 200, output_path + "/exonizations_by_sample_coverage_mut.tab")

        # 9. Separate between mutated and non-mutated cases
        logger.info("Part9...")
        command2="Rscript "+dir_path+"/lib/Exonization/separate_mutated_cases.R "+output_path + \
                 "/exonizations_by_sample_coverage_mut.tab"+" "+output_path + "/mutated_exonizations.tab"+" "+output_path + "/non_mutated_exonizations.tab"
        # print(command2)
        os.system(command2)

        # 10. Get the tumor specific events
        if(tumor_specific):

            # Get also the significant exonizations from Rudin and Intropolis
            output_Rudin_path_aux2 = output_path + "/new_exonized_junctions_Rudin_normal_reads.tab"
            readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
            get_reads_exonizations(output_path+"/new_exonized_junctions.tab", readCounts_Rudin_path, output_Rudin_path_aux2)
            output_Rudin_path_aux3 = output_path + "/new_exonized_junctions_Rudin_normal_reads_repeatitions.tab"
            overlap_with_repeats(output_Rudin_path_aux2, repeats_path, output_Rudin_path_aux3)
            output_Rudin_path_aux4 = output_path + "/exonizations_by_sample_Rudin_normal.tab"
            get_significant_exonizations(output_Rudin_path_aux3, threshold, output_Rudin_path_aux4)

            output_Intropolis_path_aux2 = output_path + "/new_exonized_junctions_Intropolis_reads.tab"
            get_reads_exonizations(output_path+"/new_exonized_junctions.tab", readcounts_path, output_Intropolis_path_aux2)
            output_Intropolis_path_aux3 = output_path + "/new_exonized_junctions_Intropolis_reads_repeatitions.tab"
            overlap_with_repeats(output_Intropolis_path_aux2, repeats_path, output_Intropolis_path_aux3)
            output_Intropolis_path_aux4 = output_path + "/exonizations_by_sample_Intropolis.tab"
            get_significant_exonizations(output_Intropolis_path_aux3, threshold, output_Intropolis_path_aux4)

            output_Rudin_path_aux4 = output_path + "/exonizations_by_sample_Rudin_normal.tab"
            output_Intropolis_path_aux4 = output_path + "/exonizations_by_sample_Intropolis.tab"
            output_path_aux11 = output_path + "/non_mutated_exonizations_filtered.tab"
            filter_exonizations(output_path + "/non_mutated_exonizations.tab", output_Rudin_path_aux4, output_Intropolis_path_aux4, output_path_aux11, flag_Rudin)
            output_path_aux12 = output_path + "/non_mutated_exonizations_filtered2.tab"
            filter_exonizations_CHESS(output_path_aux11, CHESS_SE_path, output_path_aux12)

            # 11. Join the mutated and non_mutated cases
            logger.info("Part10...")
            output_path_aux13 = output_path + "/all_exonizations.tab"
            command3 = "cat " + output_path + "/mutated_exonizations.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path_aux12 + " >> " + output_path_aux13
            os.system(command3)

        else:

            # 11. Join the mutated and non_mutated cases
            logger.info("Part10...")
            output_path_aux13 = output_path + "/all_exonizations.tab"
            command3 = "cat " + output_path + "/mutated_exonizations.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path + "/non_mutated_exonizations.tab" + " >> " + output_path_aux13
            os.system(command3)

        # 12. Get the peptide sequence associated
        logger.info("Part11...")
        output_path_aux13 = output_path + "/all_exonizations.tab"
        output_path_peptide = output_path + "/exonizations_peptide_sequence.fa"
        output_path_dna = output_path + "/exonizations_fasta_sequence.fa"
        output_path_aux14 = output_path + "/all_exonizations_ORF.tab"
        output_path_aux15 = output_path + "/all_exonizations_ORF_sequences.tab"
        output_path_aux16 = output_path + "/all_exonizations_Interpro.tab"
        output_path_aux17 = output_path + "/all_exonizations_IUPred.tab"
        get_peptide_sequence(output_path_aux13, transcript_expression_path, gtf_path,
                             output_path_peptide, output_path_dna, output_path_aux14,
                             output_path_aux15, output_path_aux16, output_path_aux17, mosea_path, genome_path, orfs_scripts,
                             interpro,IUPred, remove_temp_files, python2)

        # 13. Filter the significant results
        logger.info("Part12...")
        output_path_aux18 = output_path + "/all_exonizations_filtered.tab"
        output_path_aux19 = output_path + "/all_exonizations_filtered_peptide_change.tab"
        command4="module load R; Rscript "+dir_path+"/lib/Exonization/filter_results.R "+output_path_aux14+" "+output_path_aux18+" "+output_path_aux19
        os.system(command4)

        # 14. Select the fasta candidates for being run to the epitope analysis
        logger.info("Part13...")
        output_path_aux20 = output_path + "/exonizations_peptide_sequence.fa"
        output_path_aux21 = output_path + "/exonizations_peptide_sequence_filtered.fa"
        #Create the folder, if it doesn't exists
        if not os.path.exists(output_path + "/exonization_fasta_files"):
            os.makedirs(output_path + "/exonization_fasta_files")
        select_fasta_candidates(output_path_aux19, output_path_aux20, output_path_aux21, output_path + "/exonization_fasta_files")

        # 15. Run netMHC-4.0_part1
        logger.info("Part14...")
        if not os.path.exists(output_path + "/exonizations_NetMHC-4.0_files"):
            os.makedirs(output_path + "/exonizations_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path_aux19, HLAclass_path, HLAtypes_path,
                                      output_path + "/exonization_fasta_files",output_path + "/exonizations_NetMHC-4.0_files", output_path + "/exonizations_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/exonizations_NetMHC-4.0_neoantigens_type_3_all.tab", output_path + "/exonizations_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/exonizations_NetMHC-4.0_neoantigens_type_2_all.tab", output_path + "/exonizations_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        # 16. Run netMHCpan-4.0_part1
        logger.info("Part15...")
        if not os.path.exists(output_path + "/exonizations_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/exonizations_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path_aux19, HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/exonization_fasta_files",output_path + "/exonizations_NetMHCpan-4.0_files", output_path + "/exonizations_NetMHCpan-4.0_neoantigens_type_3.tab",
                                      output_path + "/exonizations_NetMHCpan-4.0_neoantigens_type_3_all.tab", output_path + "/exonizations_NetMHCpan-4.0_neoantigens_type_2.tab",
                                      output_path + "/exonizations_NetMHCpan-4.0_neoantigens_type_2_all.tab", output_path + "/exonizations_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_pan_path)
        logger.info("Wait until all jobs have finished. Then, go on with part3")

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.reads,args.gtf,args.genome,args.transcript,args.HLAclass,args.HLAtypes,args.HLAtypespan,args.netMHC,
         args.netMHCpan,args.mosea,args.orfs,args.output,args.repeats,args.thres,args.tumor_specific,args.mutations,
         args.chess,args.Rudin,args.temp,args.username)

    # readcounts_path = "/media/trincadojl/WINDOWS 10/Work/SCLC/ISOTOPE_TEST/data/readCounts_TEST.tab"
    # bam_path = "/media/trincadojl/WINDOWS 10/Work/SCLC/ISOTOPE_TEST/data/STAR"
    # gtf_path =  "/media/trincadojl/WINDOWS 10/Work/SCLC/ISOTOPE_TEST/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
    # genome_path =  "/media/trincadojl/data/Projects/annotation/hg19.fa"
    # mosea_path =  "/home/trincadojl/Software/MoSEA"
    # #mosea_path =  "/home/shinoda/Software/MoSEA-py3"
    # max_length = 500
    # threshold = 5
    # n_randomizations = 100
    # repeats_path = "/media/trincadojl/WINDOWS 10/Work/SCLC/ISOTOPE_TEST/annotation/hg19_repeats_TEST.bed"
    # output_path = "/home/trincadojl/Desktop/ISOTOPE_test/exonizations"
    #
    # main(readcounts_path, bam_path, gtf_path, genome_path, mosea_path, output_path, repeats_path, 500, 5, 100)
