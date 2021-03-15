"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
IR_ISOTOPE.py: get significat intron retention
"""

import os

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
from lib.IR.format_to_SPADA import *

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
parser.add_argument("-HLAclass", "--HLAclass", required=True, help="HLA genotype of the samples")
parser.add_argument("-HLAtypes", "--HLAtypes", required=True, help="HLA alelles recognized by NetMHC")
parser.add_argument("-HLAtypespan", "--HLAtypespan", required=True, help="HLA alelles recognized by NetMHCpan")
parser.add_argument("-netMHC", "--netMHC", required=True, help="netMHC path")
parser.add_argument("-netMHCpan", "--netMHCpan", required=True, help="netMHCpan path")
parser.add_argument("-o", "--output", required=True, help="Output path")

def main(HLAclass_path, HLAtypes_path, HLAtypes_pan_path, netMHC_path, netMHC_pan_path, output_path):

    try:

        logger.info("Starting execution IR_ISOTOPE_part3")

        # HLAclass_path = "/projects_rg/SCLC_cohorts/tables/PHLAT_summary_ClassI_all_samples.out"
        # HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        # HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        # netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        # netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/SCLC/epydoor/IR"

        #13. Run netMHC-4.0_part2
        logger.info("Part13...")
        run_netMHC_classI_slurm_part2(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHC-4.0_files", output_path + "/IR_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_2_all.tab", output_path + "/IR_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        #14. Run netMHCpan-4.0_part2
        logger.info("Part14...")
        run_netMHCpan_classI_slurm_part2(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHCpan-4.0_files", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2_all.tab", output_path + "/IR_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_pan_path)

        # 15. Run format_to_SPADA
        # logger.info("Part15...")
        # format_to_SPADA(output_path + "/IR_ORF.tab", output_path + "/IR_ORF_sequences.tab", output_path + "/IR_ORF_Interpro.tab",
        #                 output_path + "/IR_ORF_IUPred.tab", output_path + "/IR_SPADA.tab", output_path + "/IR_SPADA.fasta",
        #                 output_path + "/IR_SPADA_features.tab")

        logger.info("Done.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.HLAclass,args.HLAtypes,args.HLAtypespan,args.netMHC,args.netMHCpan,args.output)
