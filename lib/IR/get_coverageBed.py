"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_coverageBed: get the read counts per randomization. Get a significance for each intron

"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
from statsmodels.distributions.empirical_distribution import ECDF


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

# description = \
#     "Description:\n\n"
#
# parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
#                         add_help=True)
# parser.add_argument("-i", "--input", required=True,
#                     help="Input file")
# parser.add_argument("-g", "--gtf", required=True,
#                     help="Gtf file")
# parser.add_argument("-l", "--length", required=False, type=int, default=0,
#                     help="Max length for the exonizations")
# parser.add_argument("-o", "--output", required=True,
#                     help="Output file")

def extract_number(id):
    '''
    Extract the number id. If there is not, return a 0:
    '''
    try:
        return int(id.split("_")[3])
    except:
        return 0


def get_coverageBed(input_path, gtf_path, coverage_path, output_path, sclc_flag):

    # args = parser.parse_args()

    try:
        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # coverage_path = sys.argv[3]
        # output_path = sys.argv[4]
        # sclc_flag = sys.argv[5]

        # input_path = "/projects_rg/SCLC_cohorts/Smart/IR/IR_significant_genes_filtered2.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/Smart/IR/random_introns.bed"
        # coverage_path = "/projects_rg/SCLC_cohorts/Smart/IR/coverageBed"
        # output_path = "/projects_rg/SCLC_cohorts/Smart/IR/IR_significant_genes_filtered3.tab"
        # sclc_flag = False

        # input_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered2.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/random_introns.bed"
        # coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/intron_retention/IR_v5"
        # output_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/IR_significant_genes_filtered3.tab"

        # Load the gtf file
        logger.info("Loading gtf with the randomizations...")
        gtf = pd.read_table(gtf_path, delimiter="\t")
        gtf.columns = ['chr','start','end','id','strand','score']

        # Load the input file
        logger.info("Loading introns...")
        introns = pd.read_table(input_path, delimiter="\t")
        # Sort by sample id
        introns = introns.sort_values("Sample_id")
        aux_sample = ""
        pvalue_list, reads_coverage_intron_list = [], []
        for i in range(0,len(introns.index)):
            logger.info(str(i)+"...")
            new_sample = introns["Sample_id"].iloc[i].rstrip()
            #Get the id of this intron
            intron = introns["Event_id"].iloc[i].rstrip()
            # intron = introns["IR"].iloc[i].rstrip()
            chr_intron = intron.split(":")[0]
            start_intron = int(intron.split(":")[1].split("-")[0])
            end_intron = int(intron.split(":")[1].split("-")[1][:-3])
            id_intron = gtf.loc[(gtf['chr'] == chr_intron) & (gtf['start'] == start_intron) &
                                     (gtf['end'] == end_intron),"id"]
            if(id_intron.size!=0):
                main_n_intron = int(id_intron.apply(lambda x: x.split("_")[1]).iloc[0])
                #If it's a new sample, load the new file
                if(aux_sample!=new_sample):
                    aux_sample = new_sample
                    #If SCLC samples analysis is executed
                    if(sclc_flag=="True"):
                        sample_formatted = aux_sample.replace("T", "").replace("X", "").replace(".", "-")
                    else:
                        sample_formatted = aux_sample.replace(".", "-")
                    # sample_formatted = aux_sample
                    logger.info("Processing sample "+sample_formatted+"...")
                    #Load the corresponding coverageBed file
                    file = coverage_path + "/" + sample_formatted + ".coverage_sorted"
                    # coverage_file = pd.read_table(file, delimiter="\t", header=None, skiprows=1)
                    coverage_file = pd.read_table(file, delimiter="\t", header=None)
                    coverage_file.columns = ["chr","start","end","id","strand","score","n_reads","n_bases","length","fraction_bases"]
                    #Create an id like the intron
                    if(sclc_flag=="True"):
                        id_exons = coverage_file["chr"] + ";" + list(map(str,coverage_file["start"])) + ";" + \
                                   list(map(str,coverage_file["end"])) + ";" + coverage_file["strand"]
                    else:
                        id_exons = "chr" + coverage_file["chr"].apply(str) + ";" + coverage_file["start"].apply(str) + ";" + \
                                   coverage_file["end"].apply(str) + ";" + coverage_file["strand"]
                    coverage_file["id_exons"] = id_exons
                    # Extract the number of the random introns
                    n_intron = coverage_file["id"].apply(lambda x: int(x.split("_")[1]))
                    n_randomization = coverage_file["id"].apply(extract_number)
                    coverage_file["n_intron"] = n_intron
                    coverage_file["n_randomization"] = n_randomization
                    #Sort the df by the two previous numbers
                    coverage_file_sorted = coverage_file.sort_values(["n_intron","n_randomization"])

                # Get the number of reads of the intron and their randomizations and do the test
                reads_intron = int(coverage_file_sorted.loc[(coverage_file_sorted['n_intron'] == main_n_intron) &
                                                             (coverage_file_sorted['n_randomization'] == 0)]["n_reads"].iloc[0])
                reads_randomizations = coverage_file_sorted.loc[(coverage_file_sorted['n_intron'] == main_n_intron) &
                                                                (coverage_file_sorted['n_randomization'] != 0)]["n_reads"].tolist()
                ecdf = ECDF(reads_randomizations)
                pvalue = (1.0 - ecdf(reads_intron))
                pvalue_list.append(pvalue)
                reads_coverage_intron_list.append(reads_intron)

            else:   #No generated introns for this event
                pvalue_list.append(1)
                reads_coverage_intron_list.append(0)

        introns["reads_intron"] = reads_coverage_intron_list
        introns["pvalue_coverage"] = pvalue_list
        # Save the file
        introns.to_csv(output_path, sep="\t", index=False)
        logger.info("Saved "+output_path)
        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    input_path = sys.argv[1]
    gtf_path = sys.argv[2]
    coverage_path = sys.argv[3]
    output_path = sys.argv[4]
    sclc_flag = sys.argv[5]
    get_coverageBed(input_path, gtf_path, coverage_path, output_path, sclc_flag)