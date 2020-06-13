"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
A5_A3_ISOTOPE.py: get significat A5_A3
"""

from lib.A5_A3.run_netMHC_classI_slurm_part1 import *
from lib.A5_A3.run_netMHCpan_classI_slurm_part1 import *
from lib.A5_A3.run_netMHC_classI_slurm_part2 import *
from lib.A5_A3.run_netMHCpan_classI_slurm_part2 import *
from lib.A5_A3.format_to_SPADA import *

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

        HLAclass_path = "/projects_rg/SCLC_cohorts/tables/PHLAT_summary_ClassI_all_samples.out"
        HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/SCLC/epydoor/A5_A3"


        #17. Run netMHC-4.0_part2
        logger.info("Part13...")
        run_netMHC_classI_slurm_part2(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/A5_A3_fasta_files",
                                      output_path + "/A5_A3_NetMHC-4.0_files",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        #18. Run netMHCpan-4.0_part2
        logger.info("Part14...")
        run_netMHCpan_classI_slurm_part2(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path,
                                         HLAtypes_pan_path,
                                         output_path + "/A5_A3_fasta_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                         netMHC_pan_path)

        # 19. Run format_to_SPADA
        logger.info("Part18...")
        format_to_SPADA(output_path + "/A5_A3_ORF.tab", output_path + "/A5_A3_ORF_sequences.tab",
                        output_path + "/A5_A3_Interpro.tab",
                        output_path + "/A5_A3_IUPred.tab", output_path + "/A5_A3_SPADA.tab",
                        output_path + "/A5_A3_SPADA.fasta", output_path + "/A5_A3_SPADA_features.tab")
        logger.info("Done.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()
