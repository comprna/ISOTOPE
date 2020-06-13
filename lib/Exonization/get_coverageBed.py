"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

get_coverageBed: get the read counts per randomization. Get a significance for each exonization

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
# parser.add_argument("-c", "--coverage", required=True,
#                     help="Coverage file")
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

        # input_path = args.input
        # gtf_path = args.gtf
        # coverage_path = args.coverage
        # output_path = args.output

        # input_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample.tab"
        # gtf_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/random_exonizations.bed"
        # coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/"
        # output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample_coverage.tab"

        # Load the gtf file
        logger.info("Loading gtf with the exonizations...")
        gtf = pd.read_table(gtf_path, delimiter="\t")
        gtf.columns = ['chr','start','end','id','strand','score']

        # Load the input file
        logger.info("Loading exonizations...")
        exonizations = pd.read_table(input_path, delimiter="\t")
        # Sort by sample id
        exonizations = exonizations.sort_values("Sample_id")
        aux_sample = ""
        pvalue_list, reads_coverage_exonization_list = [], []
        for i in range(0,len(exonizations.index)):
            new_sample = exonizations["Sample_id"].iloc[i].rstrip()
            #Get the id of this exonization
            exonization = exonizations["New_exon"].iloc[i].rstrip()
            chr_exonization = exonization.split(";")[0]
            start_exonization = int(exonization.split(";")[1])
            end_exonization = int(exonization.split(";")[2])
            id_exonization = gtf.loc[(gtf['chr'] == chr_exonization) & (gtf['start'] == start_exonization) &
                                     (gtf['end'] == end_exonization),"id"]
            main_n_exonization = int(id_exonization.apply(lambda x: x.split("_")[1]).iloc[0])
            #If it's a new sample, load the new file
            if(aux_sample!=new_sample):
                aux_sample = new_sample
                # If SCLC samples analysis is executed
                if (sclc_flag == "True"):
                    sample_formatted = aux_sample.replace("T", "").replace("X", "").replace(".", "-")
                else:
                    sample_formatted = aux_sample.replace(".", "-")
                logger.info("Processing sample "+aux_sample+"...")
                #Load the corresponding stringtie file
                file = coverage_path + "/" + aux_sample + ".coverage_sorted"
                coverage_file = pd.read_table(file, delimiter="\t", header=None, skiprows=1)
                coverage_file.columns = ["chr","start","end","id","strand","score","n_reads","n_bases","length","fraction_bases"]
                #Create an id like the exonization
                id_exons = coverage_file["chr"] + ";" + list(map(str,coverage_file["start"])) + ";" + \
                           list(map(str,coverage_file["end"])) + ";" + coverage_file["strand"]
                coverage_file["id_exons"] = id_exons
                # Extract the number of the random exonizations
                n_exonization = coverage_file["id"].apply(lambda x: int(x.split("_")[1]))
                n_randomization = coverage_file["id"].apply(extract_number)
                coverage_file["n_exonization"] = n_exonization
                coverage_file["n_randomization"] = n_randomization
                #Sort the df by the two previous numbers
                coverage_file_sorted = coverage_file.sort_values(["n_exonization","n_randomization"])

            # Get the number of reads of the exonization and their randomizations and do the test
            reads_exonization = int(coverage_file_sorted.loc[(coverage_file_sorted['n_exonization'] == main_n_exonization) &
                                                         (coverage_file_sorted['n_randomization'] == 0)]["n_reads"].iloc[0])
            reads_randomizations = coverage_file_sorted.loc[(coverage_file_sorted['n_exonization'] == main_n_exonization) &
                                                            (coverage_file_sorted['n_randomization'] != 0)]["n_reads"].tolist()
            ecdf = ECDF(reads_randomizations)
            pvalue = (1.0 - ecdf(reads_exonization))
            pvalue_list.append(pvalue)
            reads_coverage_exonization_list.append(reads_exonization)

        exonizations["reads_exonization"] = reads_coverage_exonization_list
        exonizations["pvalue"] = pvalue_list
        # Save the file
        exonizations.to_csv(output_path, sep="\t", index=False)
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