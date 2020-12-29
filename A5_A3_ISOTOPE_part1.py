"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
A5_A3_ISOTOPE.py: get significant alternative splice site events
"""


from lib.A5_A3.extract_exonized_junctions import *
from lib.A5_A3.get_reads_exonizations import *
from lib.A5_A3.overlap_with_repeats import *
from lib.A5_A3.get_significant_exonizations import *
from lib.A5_A3.compare_reads_random_junctions import *
from lib.A5_A3.check_mutations_nearby import *
from lib.A5_A3.filter_exonizations import *
from lib.A5_A3.filter_exonizations_CHESS import *
from lib.A5_A3.get_peptide_sequence import *
from lib.A5_A3.select_fasta_candidates import *
from lib.A5_A3.run_netMHC_classI_slurm_part1 import *
from lib.A5_A3.run_netMHCpan_classI_slurm_part1 import *
from itertools import chain, islice
import os

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
"Description: Get alternative splice site events\n\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-r", "--reads", required=True, help = "reads mapped to junctions")
parser.add_argument("-trans", "--transcript", required=True, help = "transcript expression file")
parser.add_argument("-g", "--gtf", required=True, help = "gtf annotation")
parser.add_argument("-c", "--conversion", required=True, help = "gene name conversion")
parser.add_argument("-m", "--max", required=False, type=int, default=500)
parser.add_argument("-t", "--thres", required=False, type=int, default=5, help="Minimum number of reads mapping the event")
parser.add_argument("-chunks", "--size_chunks", required=False, type=int, default=1, help="For paralellization, this values indicates number of jobs to run")
parser.add_argument("-rep", "--repeats", required=True, help = "Regions of the genome with repeats from maskerDB",default=None)
parser.add_argument("-mut","--mutations", required=False, help = "Mutations path")
parser.add_argument("--chessA5", required=False, help = "CHESS A5 path")
parser.add_argument("--chessA3", required=False, help = "CHESS A3 path")
parser.add_argument("--tumor_specific", type=str2bool, nargs='?',const=True, default=False,help="Tumor specific mode")
parser.add_argument("-mosea", "--mosea", required=True, help = "MoSEA path")
parser.add_argument("-orfs", "--orfs", required=True, help = "MxFinder path")
parser.add_argument("-genome", "--genome", required=True, help = "Genome annotation")
parser.add_argument("-HLAclass", "--HLAclass", required=True, help = "HLA genotype of the samples")
parser.add_argument("-HLAtypes", "--HLAtypes", required=True, help = "HLA alelles recognized by NetMHC")
parser.add_argument("-HLAtypespan", "--HLAtypespan", required=True, help = "HLA alelles recognized by NetMHCpan")
parser.add_argument("-netMHC", "--netMHC", required=True, help = "netMHC path")
parser.add_argument("-netMHCpan", "--netMHCpan", required=True, help = "netMHCpan path")
parser.add_argument("--temp", type=str2bool, nargs='?',const=True, default=False,help="Remove temp files")
parser.add_argument("--Rudin", type=str2bool, nargs='?',const=True, default=False,help="Rudin mode")
parser.add_argument("--username", required=True, help = "Cluster user name")

def chunks(iterable, n):
   "chunks(ABCDE,2) => AB CD E"
   iterable = iter(iterable)
   while True:
       yield chain([next(iterable)], islice(iterable, n-1))


def main(readcounts_path, transcript_expression_path, gtf_path, conversion_names, max_length,
         threshold, size_chunks, repeats_path, mutations_path, CHESS_A5_path, CHESS_A3_path,
         tumor_specific, mosea, orfs_scripts, fasta_genome, HLAclass_path, HLAtypes_path,
         HLAtypes_pan_path, netMHC_path, netMHC_pan_path, remove_temp_files, flag_Rudin,
         output_path, name_user):

    try:

        logger.info("Starting execution")

        # readcounts_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/readCounts_George_Peifer_Rudin_Yokota.tab"
        # transcript_expression_path = "/projects_rg/SCLC_cohorts/George/tables/iso_tpm_George_Peifer_Rudin_Yokota.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        # conversion_names = "/projects_rg/SCLC_cohorts/tables/Ensembl_gene_conversion.txt"
        # max_length = 500
        # threshold = 5
        # size_chunks = 100
        # repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        # mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        # CHESS_A5_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe"
        # CHESS_A3_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe"
        # tumor_specific = False
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
        # flag_Rudin = False
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/SCLC/epydoor/A5_A3"
        # name_user = "juanluis"

        # 1. Identify the junctions that could generate an alternative splice site
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_A5_A3_junctions.tab"
        extract_exonized_junctions(readcounts_path, gtf_path, fasta_genome, max_length, output_path_aux, mosea)

        # 2. Given the list with the possible A5_A3, get the reads associate to each of them
        logger.info("Part2...")
        output_path_aux2 = output_path+"/new_A5_A3_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)

        # 3. find the overlap between the nex A5_A3 and repeatitions (RepeatMasker)
        logger.info("Part3...")
        output_path_aux3 = output_path + "/new_A5_A3_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the A5_A3 with the reads counts,get those that are over a threshold
        logger.info("Part4...")
        get_significant_exonizations(output_path_aux3, threshold, output_path + "/A5_A3_by_sample.tab")

        # 5. for applying some filtering on the list of A5_A3 junctions, we are gonna compare the readcounts for each
        # junction against other new junctions associated to the same gene
        logger.info("Part5...")
        compare_reads_random_junctions(output_path + "/A5_A3_by_sample.tab", readcounts_path, gtf_path, output_path + "/A5_A3_by_sample_coverage.tab")

        # 6. Check if in the A5_A3 there are mutations nearby
        logger.info("Part6...")
        check_mutations_nearby(output_path + "/A5_A3_by_sample_coverage.tab", mutations_path, 200, output_path + "/A5_A3_by_sample_coverage_mut.tab")

        # 7. Separate the mutated from the non-mutated cases
        logger.info("Part7...")
        command1="Rscript "+dir_path+"/lib/A5_A3/separate_mutated_cases.R "+ output_path + "/A5_A3_by_sample_coverage_mut.tab " \
                 + output_path + "/A5_A3_mutated.tab " + output_path + "/A5_A3_non_mutated.tab "
        os.system(command1)

        # 8. Get the tumor specific events
        if(tumor_specific):

            # Get also the significant A5_A3 from Rudin and Intropolis
            logger.info("Part8.1...")
            output_Rudin_path_aux2 = output_path + "/new_A5_A3_junctions_Rudin_normal_reads.tab"
            readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
            get_reads_exonizations(output_path+"/new_A5_A3_junctions.tab", readCounts_Rudin_path, output_Rudin_path_aux2)
            output_Rudin_path_aux3 = output_path + "/new_A5_A3_junctions_Rudin_normal_reads_repeatitions.tab"
            overlap_with_repeats(output_Rudin_path_aux2, repeats_path, output_Rudin_path_aux3)
            output_Rudin_path_aux4 = output_path + "/A5_A3_by_sample_Rudin_normal.tab"
            get_significant_exonizations(output_Rudin_path_aux3, threshold, output_Rudin_path_aux4)

            logger.info("Part8.2...")
            output_Intropolis_path_aux2 = output_path + "/new_A5_A3_junctions_Intropolis_reads.tab"
            get_reads_exonizations(output_path+"/new_A5_A3_junctions.tab", readcounts_path, output_Intropolis_path_aux2)
            output_Intropolis_path_aux3 = output_path + "/new_A5_A3_junctions_Intropolis_reads_repeatitions.tab"
            overlap_with_repeats(output_Intropolis_path_aux2, repeats_path, output_Intropolis_path_aux3)
            output_Intropolis_path_aux4 = output_path + "/A5_A3_by_sample_Intropolis.tab"
            get_significant_exonizations(output_Intropolis_path_aux3, threshold, output_Intropolis_path_aux4)

            logger.info("Part8.3...")
            output_Rudin_path_aux4 = output_path + "/A5_A3_by_sample_Rudin_normal.tab"
            output_Intropolis_path_aux4 = output_path + "/A5_A3_by_sample_Intropolis.tab"
            output_path_aux11 = output_path + "/A5_A3_non_mutated_filtered.tab"
            filter_exonizations(output_path + "/A5_A3_non_mutated.tab", output_Rudin_path_aux4, output_Intropolis_path_aux4, output_path_aux11, flag_Rudin)
            output_path_aux12 = output_path + "/A5_A3_non_mutated_filtered2.tab"
            filter_exonizations_CHESS(output_path_aux11, CHESS_A5_path, CHESS_A3_path, output_path_aux12)

            # 9. Join the mutated and non_mutated cases
            logger.info("Part8.4...")
            output_path_aux13 = output_path + "/all_A5_A3.tab"
            command3 = "cat " + output_path + "/A5_A3_mutated.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path_aux12 + " >> " + output_path_aux13
            os.system(command3)

        else:

            # 9. Join the mutated and non_mutated cases
            logger.info("Part8...")
            output_path_aux13 = output_path + "/all_A5_A3.tab"
            command3 = "cat " + output_path + "/A5_A3_mutated.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path + "/A5_A3_non_mutated.tab" + " >> " + output_path_aux13
            os.system(command3)

        # 10. Get the peptide sequence associated

        # 10.1. Split the input file into n pieces. Run a job per piece. When all jobs have finished, we will assemble all the pieces
        logger.info("Part9...")
        logger.info("get_peptide_sequence: Split the file into pieces and run get_peptide_sequence by chunk")
        dir_path = os.path.dirname(os.path.realpath(__file__))

        with open(output_path_aux13) as f_aux:
            header = f_aux.readline().strip()

        dict_jobs = {}
        with open(output_path_aux13) as bigfile:
            for i, lines in enumerate(chunks(bigfile, size_chunks)):
                file_split = '{}.{}'.format(output_path_aux13, i)
                f = open(file_split, 'w')
                #Output the header, if it's not the first chunk
                if(i!=0):
                    f.writelines(header)
                with f:
                    f.writelines(lines)
                #Run a job per file
                logger.info("Processing " + "chunk_" + str(i) + "...")
                command1 = "module load Python; python " + dir_path + "/get_peptide_sequence.py " + output_path_aux13 + " " + \
                transcript_expression_path + " " + gtf_path + " " + output_path + "/A5_A3_peptide_sequence.fa " + \
                output_path + "/A5_A3_fasta_sequence.fa " + output_path + "/A5_A3_ORF.tab " + output_path + "/A5_A3_ORF_sequences.tab " + \
                mosea + " " + fasta_genome + " " + orfs_scripts + " " + remove_temp_files
                open_peptides_file = open(output_path + "/aux.sh", "w")
                open_peptides_file.write("#!/bin/sh\n")
                # open_peptides_file.write("#SBATCH --partition=normal\n")
                open_peptides_file.write("#SBATCH --mem 3000\n")
                open_peptides_file.write(
                    "#SBATCH -e " + output_path + "/" + "get_peptide_sequence" + "_chunk_" + str(i) + ".err" + "\n")
                open_peptides_file.write(
                    "#SBATCH -o " + output_path + "/" + "get_peptide_sequence" + "_chunk_" + str(i) + ".out" + "\n")
                open_peptides_file.write(command1 + ";\n")
                open_peptides_file.close()
                command2 = "sbatch -J " + "get_peptide_sequence" + "_chunk_" + str(i) + " " + output_path + "/aux.sh; sleep 0.5;"
                # os.system(command2)
                job_message = subprocess.check_output(command2, shell=True)
                # Get the job id and store it
                job_id = (str(job_message).rstrip().split(" ")[-1])[:-3]
                dict_jobs[job_id] = 1

        logger.info("get_peptide_sequence: Waiting for all the jobs to finished...")
        flag_exit = False
        while (not flag_exit):
            # Initialize the dictionary with the pending jobs in the cluster
            pending_jobs = {}
            os.system("sleep 10")
            p = subprocess.Popen(["squeue", "-u", name_user], stdout=subprocess.PIPE)
            # Skip the first line (the header)
            line = p.stdout.readline()
            for line in p.stdout:
                flag_exit = True
                # Get the id of the job
                job_id_aux = str(line).rstrip().split()[1]
                # Save the id of the jobs
                pending_jobs[job_id_aux] = 1
                # If there is any job on the cluster on dict_jobs, break the loop and wait for another 10 seconds
                # to check the status of the jobs in the cluster
                if (job_id_aux in dict_jobs):
                    flag_exit = False
                    break

        logger.info("get_peptide_sequence:All jobs finished.\n\n")

        # 11. Filter the relevant results
        command4 = "Rscript " + dir_path + "/lib/A5_A3/filter_results.R " + output_path + "/A5_A3_ORF.tab " \
                   + conversion_names + " " + output_path + "/A5_A3_ORF_filtered.tab " + output_path + "/A5_A3_ORF_filtered_peptide_change.tab"
        os.system(command4)

        # 12. Select the fasta candidates for being run to the epitope analysis
        logger.info("Part10...")
        # Create the folder, if it doesn't exists
        if not os.path.exists(output_path + "/A5_A3_fasta_files"):
            os.makedirs(output_path + "/A5_A3_fasta_files")
        select_fasta_candidates(output_path + "/A5_A3_ORF_filtered_peptide_change.tab",
                                output_path + "/A5_A3_peptide_sequence.fa", output_path + "/A5_A3_peptide_sequence_filtered.fa",
                                output_path + "/A5_A3_fasta_files")

        # 13. Run netMHC-4.0_part1
        logger.info("Part11...")
        if not os.path.exists(output_path + "/A5_A3_NetMHC-4.0_files"):
            os.makedirs(output_path + "/A5_A3_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/A5_A3_fasta_files",
                                      output_path + "/A5_A3_NetMHC-4.0_files",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        # 14. Run netMHCpan-4.0_part1
        logger.info("Part12...")
        if not os.path.exists(output_path + "/A5_A3_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/A5_A3_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                         output_path + "/A5_A3_fasta_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                         netMHC_pan_path)

        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Abor  ting execution")
        sys.exit(1)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.reads,args.transcript,args.gtf,args.conversion,args.conversion,args.max,args.thres,args.size_chunks,
         args.repeats,args.mutations,args.chessA5,args.chessA3,args.tumor_specific,args.mosea,args.orfs,
         args.genome,args.HLAclass,args.HLAtypes,args.HLAtypespan,args.netMHC,args.netMHCpan,args.temp,
         args.Rudin,args.username)