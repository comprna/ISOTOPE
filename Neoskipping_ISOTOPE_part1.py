"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
Neoskipping_ISOTOPE.py: get significant neoskipping events
"""

import os

from lib.Neoskipping.extract_neoskipping_junctions import *
from lib.Neoskipping.extract_neoskipping_junctions_Intropolis import *
from lib.Neoskipping.check_mutations_nearby import *
from lib.Neoskipping.filter_neoskipping import *
from lib.Neoskipping.filter_neoskipping_CHESS import *
from lib.Neoskipping.get_peptide_sequence import *
from lib.Neoskipping.select_fasta_candidates import *
from lib.Neoskipping.run_netMHC_classI_slurm_part1 import *
from lib.Neoskipping.run_netMHCpan_classI_slurm_part1 import *

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
        transcript_expression_path = "/projects_rg/SCLC_cohorts/Smart/Salmon/iso_tpm.txt"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        conversion_names = "/projects_rg/SCLC_cohorts/tables/Ensembl_gene_conversion.txt"
        max_length = 500
        threshold = 5
        threshold2 = 10
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        CHESS_SE_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe"
        tumor_specific = True
        mosea = "/genomics/users/juanluis/Software/MoSEA-master/mosea.py"
        fasta_genome = "/genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa"
        orfs_scripts = "/genomics/users/juanluis/comprna/MxFinder/extract_orfs.py"
        interpro = "/soft/EB_repo/bio/sequence/programs/noarch/interproscan/5.33-72.0/interproscan.sh"
        IUPred = "/projects_rg/SCLC_cohorts/soft/IUPred2A"
        HLAclass_path = "/projects_rg/SCLC_cohorts/Smart/PHLAT/PHLAT_summary_ClassI.out"
        HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        remove_temp_files = True
        flag_Rudin = False
        output_path = "/users/genomics/juanluis/SCLC_cohorts/Smart/epydoor/neoskipping"
        # ONLY FOR MARVIN
        #python2 = "Python/2.7.14-foss-2017b"
        # ONLY FOR HYDRA
        python2 = "Python/2.7.11"

    # 1. Identify the junctions that could generate an alternative splice site
    logger.info("Part1...")
    dir_path = os.path.dirname(os.path.realpath(__file__))
    output_path_aux = output_path + "/new_Neoskipping_junctions.tab"
    extract_neoskipping_junctions(readcounts_path, gtf_path, threshold, output_path_aux)

    # 2. Get the tumor specific neoskipping events
    if (tumor_specific):

        logger.info("Part2...")
        # Get also the significant neoskipping from Rudin and Intropolis
        output_Rudin_path_aux = output_path + "/new_Neoskipping_junctions_Rudin_normal_reads.tab"
        readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
        extract_neoskipping_junctions(readCounts_Rudin_path, gtf_path, threshold, output_Rudin_path_aux)

        output_Intropolis_path_aux = output_path + "/new_Neoskipping_junctions_Rudin_normal_reads.tab"
        readCounts_Intropolis_path = "/projects_rg/Annotation/Intropolis/intropolis.v1.hg19.filtered.tsv"
        extract_neoskipping_junctions_Intropolis(readcounts_path, readCounts_Intropolis_path, gtf_path, threshold,
                                                 output_Intropolis_path_aux)

        filter_neoskipping(output_path_aux, output_Rudin_path_aux, output_Intropolis_path_aux,
                           output_path + "/new_Neoskipping_junctions_filtered.tab", flag_Rudin)
        filter_neoskipping_CHESS(output_path + "/new_Neoskipping_junctions_filtered.tab", CHESS_SE_path,
                                 output_path + "/new_Neoskipping_junctions_filtered2.tab")
        output_path2 = output_path + "/new_Neoskipping_junctions_filtered2.tab"

    else:
        output_path2 = output_path + "/new_Neoskipping_junctions.tab"

    # 3. Get the mutations nearby
    logger.info("Part3...")
    check_mutations_nearby(output_path2, mutations_path, 200, output_path + "/new_Neoskipping_junctions_mut.tab")

    # 4. Get the gene ids
    logger.info("Part4...")
    command1 = "module load R; Rscript " + dir_path + "/lib/Neoskipping/get_Gene_ids_BiomaRt.R " + output_path + "/new_Neoskipping_junctions_mut.tab " + output_path + "/new_Neoskipping_junctions_mut2.tab"
    os.system(command1)

    # 5. Get the peptide sequences
    logger.info("Part5...")
    output_path_peptide = output_path + "/neoskipping_peptide_sequence.fa"
    output_path_dna = output_path + "/neoskipping_fasta_sequence.fa"
    output_path_aux14 = output_path + "/all_neoskipping_ORF.tab"
    output_path_aux15 = output_path + "/all_neoskipping_ORF_sequences.tab"
    output_path_aux16 = output_path + "/all_neoskipping_Interpro.tab"
    output_path_aux17 = output_path + "/all_neoskipping_IUPred.tab"
    get_peptide_sequence(output_path + "/new_Neoskipping_junctions_mut2.tab", transcript_expression_path, gtf_path,
                         codons_gtf_path,
                         output_path_peptide, output_path_dna, output_path_aux14,
                         output_path_aux15, output_path_aux16, output_path_aux17, mosea, fasta_genome, orfs_scripts,
                         interpro, IUPred, remove_temp_files, python2)

    # 6. Filter the cases for running netMHC
    logger.info("Part6...")
    output_path_aux18 = output_path + "/all_neoskipping_filtered.tab"
    command2 = "module load R; Rscript " + dir_path + "/lib/Neoskipping/filter_results.R " + output_path_aux14 + " " + output_path_aux18 + " " + output_path + "/all_neoskipping_filtered_peptide_change.tab"
    os.system(command2)

    # 7. Select the fasta candidates for being run to the epitope analysis
    logger.info("Part7...")
    output_path_aux20 = output_path + "/neoskipping_peptide_sequence.fa"
    output_path_aux21 = output_path + "/neoskipping_peptide_sequence_filtered.fa"
    # Create the folder, if it doesn't exists
    if not os.path.exists(output_path + "/neoskipping_fasta_files"):
        os.makedirs(output_path + "/neoskipping_fasta_files")
    select_fasta_candidates(output_path + "/all_neoskipping_filtered_peptide_change.tab", output_path_aux20,
                            output_path_aux21, output_path + "/neoskipping_fasta_files")

        # 8. Run netMHC-4.0_part1
        logger.info("Part8...")
        if not os.path.exists(output_path + "/neoskipping_NetMHC-4.0_files"):
            os.makedirs(output_path + "/neoskipping_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path + "/all_neoskipping_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/neoskipping_fasta_files",output_path + "/neoskipping_NetMHC-4.0_files", output_path + "/neoskipping_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/neoskipping_NetMHC-4.0_neoantigens_type_3_all.tab", output_path + "/neoskipping_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/neoskipping_NetMHC-4.0_neoantigens_type_2_all.tab", output_path + "/neoskipping_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        # 9. Run netMHCpan-4.0_part1
        logger.info("Part9...")
        if not os.path.exists(output_path + "/neoskipping_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/neoskipping_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path + "/all_neoskipping_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/neoskipping_fasta_files",output_path + "/neoskipping_NetMHCpan-4.0_files", output_path + "/neoskipping_NetMHCpan-4.0_neoantigens_type_3.tab",
                                      output_path + "/neoskipping_NetMHCpan-4.0_neoantigens_type_3_all.tab", output_path + "/neoskipping_NetMHCpan-4.0_neoantigens_type_2.tab",
                                      output_path + "/neoskipping_NetMHCpan-4.0_neoantigens_type_2_all.tab", output_path + "/neoskipping_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_pan_path)

        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()

