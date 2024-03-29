# ISOTOPE (ISOform-guided prediction of epiTOPEs in cancer)

![ISOTOPE_pipeline.jpg](https://user-images.githubusercontent.com/23315833/117435640-5be4a680-af2e-11eb-823a-01865a490134.png)

ISOTOPE identifies cancer-specific splicing-derived epitopes from short-read RNA-seq data. The pipeline
operates on individual tumor samples, without the requirement of additional controls or multiple tumor samples. 

The full description and testing of ISOTOPE and its applications to melanoma and small cell lung cancer samples can be found in our publication:

* Trincado JL, Reixachs-Solé M, Pérez-Granado J, Fugmann T, Sanz F, Yokota J, Eyras E. ISOTOPE: ISOform-guided prediction of epiTOPEs in cancer. PLoS Comput Biol. 2021 Sep 16;17(9):e1009411. https://doi.org/10.1371/journal.pcbi.1009411.


The ISOTOPE pipeline is divided into 4 parts, depending of the type of alternative splicing event of interest:

   * Pseudoexons (**Exonizations**)
   * New exons skipping events (**Neoskipping**)
   * Alternative splice site (**A5_A3**)
   * Intron retention (**IR**)
   
To obtain exonizations, neoskipping, and A5_A3 events, the input is a file with the counts of supporting reads for all posible splicing junctions in the genome. This file (readCounts.tab) is created with Junckey (https://github.com/comprna/Junckey#1-format-star-output). From these junctions, ISOTOPE will obtain all splicing events that show significant expression based on the sequence coverage of the gene locus.

For the identification of IR events, a normalized expression like TPMs for all possible intronic regions is needed. ISOTOPE first creates a transcriptome with all possible intronic regions using kma (https://github.com/pachterlab/kma) and then quantifies this transcriptome with a pseudoalligner (e.g. Salmon https://combine-lab.github.io/salmon/ or Kallisto https://pachterlab.github.io/kallisto/about). With these values, ISOTOPE filters out intron retention events that are lowly expressed.


The workflow is quite similar across event types, but there are some specificities that are important to take into account 
   
The user must run each of the 3 parts sequentially. The pipeline could be run on a single computer or in a slurm cluster with multithreading. If the user has several samples to analyze we strongly recommend the use of a cluster since some parts are computationally intensive.

We provide an extensive tutorial detailling how to run each of these steps. If you have a query related to the tool or any problem regarding running the pipeline, please create an issue in the repository. 

# How to run ISOTOPE

ISOTOPE allows you to obtain exonizations, neoskipping, and alternative 5'/3' splice site (A5_A3) events from a file of splicing junction from sequencing reads, and evaluate their immunogenic potential. In addition, it also evaluates the immunogenic potential of introns from retained-intron abundance files. In this tutorial, we provide all the necessary steps to recover each type of event. We provide some toy datasets for testing. If there is anything unclear or not working feel free to post an issue on GitHub and we will respond as soon as possible.

To obtain exonizations, neoskipping, and A5_A3 events, the first required input is a file with the read counts for all possible junctions in the genome. This file (readCounts.tab) is created with Junckey (https://github.com/comprna/Junckey#1-format-star-output). From these junctions, ISOTOPE will obtain splicing events.

## Exonizations

We have created a toy dataset ([readCounts_TEST.tab](https://github.com/comprna/ISOTOPE_supplementary/blob/main/readCounts_TEST.tab)) with some of the genes for which we have detected neoepitopes from splicing events (source data from Smart et al, 2018). We also provide some toy bam files obtained from STAR to run with the tutorial ([https://figshare.com/s/e42fb6db093bdb1ef046](https://figshare.com/s/e42fb6db093bdb1ef046)). Changes in the putative coding sequence created by the splicing alterations detected could give rise to new regions potentially bound by the MHC complex, increasing the likelihood to be recognized as neoantigens. This file has been generated with Junckey (https://github.com/comprna/Junckey#1-format-star-output). Visit Junckey GitHub page to see how to extract junction read counts.

### Part1:

In this first part, ISOTOPE will extract all the significant exonizations. Any exon whose junction boundaries are supported by a number of reads above the selected threshold (by default 5) will be recovered. This threshold could be changed by --thres. To assess the significance of each predicted exonization, ISOTOPE generates random intronic positions of the gene locus in which the exonization is located and compares the coverage of the exonizations of interest with that from the random introns. This step is computationally intensive and has been prepared to be parallelized (see the option **--cluster** below). When all jobs have finished, we can proceed with part2.

The pipeline uses MoSEA software (https://github.com/comprna/MoSEA) to test for splice site motifs. A genome fasta file needs to be provided (see option **--genome** below).

Additionally, if a file of repetitive elements is provided, ISOTOPE will look if the new exonizations fall on repetitive regions or regions with low complexity. Here we provide a [test bed file](https://github.com/comprna/ISOTOPE_supplementary/blob/main/hg19_repeats_TEST_f.bed) obtained from RepeatMasker (http://www.repeatmasker.org/).

**ISOTOPE_path** is the path where the whole GitHub repository is located. **data_path** is the root path where the data will be stored. **output_path** is the path where the output files will be saved. Users should change these paths according to their system folders. We will refer to these paths in the rest of the tutorial.

```
ISOTOPE_path="/users/genomics/juanluis/Software/ISOTOPE"
data_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE"
output_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE/exonizations"

python "${ISOTOPE_path}"/exonizations_ISOTOPE_part1.py -r "${data_path}"/data/readCounts_TEST.tab -b "${data_path}"/data/STAR --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa -mosea /users/genomics/juanluis/Software/MoSEA --rep "${data_path}"/annotation/hg19_repeats_TEST_f.bed -o "${output_path}"
```

We detail here the description for each of the parameters:

- **-r**  | **--reads**: Required. Reads mapped to junctions (file obtained through [Junckey](https://github.com/comprna/Junckey#1-format-star-output))

- **-o**  | **--output**: Required. output path

- **-b**  | **--bam**: Required. Path to STAR output. ISOTOPE will make use of the mapped reads to do the bootstrapping tests

- **-g**  | **--gtf**: Required. GTF annotation file. Only the exonic regions will be used.

- **-genome**  | **--genome**: Required. Genome annotation file.

- **-mosea**  | **--mosea**: Required. [MoSEA](https://github.com/comprna/MoSEA) path

- **-m**  | **--max**: Maximum length for the exonizations. Default: 500nt.

- **-t**  | **--thres**: Minimum number of reads mapping the event. Default: 5 reads.

- **-rand**  | **--rand**: Number of rounds for bootstrap testing. Default: 100.

- **-rep**  | **--repeats**: Regions of the genome with repeats from maskerDB

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is run on a slurm cluster. Default: False (it will run locally).

- **-h**  | **--help**: display the help message describing the different paramenters



### Part2:

In part 2, ISOTOPE tests the coverage of the exonizations with respect to the intronic regions using an empirical distribution. If a file with genomic mutations is provided, ISOTOPE will look for mutations falling nearby the exonizations (-mut | --mutations option (TODO: provide example file). For all the significant exonizations, ISOTOPE will evaluate if they would produce new open reading frames (ORFs) by comparing against the canonical associated transcript (all_exonizations_ORF.tab). 

If selected by the user (--tumor_specific), ISOTOPE will extract the tumor-specific events by comparing them against the GTEX database (https://www.gtexportal.org/home/) and Intropolis database (https://github.com/nellore/intropolis). We have extracted the SE values from GTEX using SUPPA (https://github.com/comprna/SUPPA). We provide this file ([chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe](https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe)). We also have extracted from Intropolis all junctions related with normal samples. This file is available via figshare, cause size issues (https://figshare.com/s/ae2283696db044488c32). Additionally, the user could provide some extra control files to compare with (--control_path). They should be in the same format as readCounts_TEST.tab. Any exonizations over the threshold count on these controls will be censored from our exonization candidates.

Quantification of isoform abundances for the same samples being tested should be provided (-trans | --transcript). This is used to ensure the observed splicing changes occur on minimally expressed isoforms. This information is also used to identify the most expressed transcript, which is defined as the "wild-type" and is used as the reference to check if the splicing change gives rise to changes in the open reading frame. This file could be computed with pseudoaligners like Salmon (https://combine-lab.github.io/salmon/) or Kallisto (https://pachterlab.github.io/kallisto/about). We provide this [file](https://github.com/comprna/ISOTOPE_supplementary/blob/main/iso_tpm.txt) for the samples tested in the tutorial. 

We make use of another software developed in the lab, MxFinder (https://github.com/JLTrincado/MxFinder). This package allows you to obtain the resulting ORF sequences produced as a consequence of the splicing changes.

For all those events giving rise to a new open reading frame and not likely to be degraded by the non-sense mediated decay (NMD) pathway (filtered by checking the 50nt rule), using NetMHC (https://services.healthtech.dtu.dk/service.php?NetMHC-4.0) and NetMHCpan (https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) we evaluate the likelihood the translated peptides could generate epitopes bound by the MHC complex. We used both NetMHC and NetMHCpan. Both methods are very similar but there are some HLA types that are only covered by one of the methods. Users should feel free to choose either one or both outputs to analyze the predicted epitopes. HLA types accepted by each method are provided in this tutorial ([HLAtypes](https://github.com/comprna/ISOTOPE_supplementary/blob/main/NetMHC-4.0_HLA_types_accepted.tab), [HLAtypespan](https://github.com/comprna/ISOTOPE_supplementary/blob/main/NetMHCpan-4.0_HLA_types_accepted.tab)). In addition, the HLA genotyping for your samples needs to be provided. This could be computed through PHLAT (https://sites.google.com/site/phlatfortype/). We include here [HLAclass](https://github.com/comprna/ISOTOPE_supplementary/blob/main/PHLAT_summary_ClassI.out) for the tested samples 

Important note: this step is computationally intensive and has been parallelized (see option **--cluster**). When all jobs have finished, we can proceed with the last part3:

```
python "${ISOTOPE_path}"/exonizations_ISOTOPE_part2.py -r "${data_path}"/data/readCounts_TEST.tab --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa --transcript /projects_rg/SCLC_cohorts/Smart/Salmon/iso_tpm.txt --HLAclass /projects_rg/SCLC_cohorts/Smart/PHLAT/PHLAT_summary_ClassI.out --HLAtypes /projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan /projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC --netMHCpan /projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan -mosea /users/genomics/juanluis/Software/MoSEA --mxfinder /genomics/users/juanluis/comprna/MxFinder --rep "${data_path}"/annotation/hg19_repeats.bed -o "${output_path}" --username juanluis --Intropolis "${data_path}"/annotation/intropolis.v1.hg19.filtered.tsv --tumor_specific --chess "${data_path}"/annotation/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe
```

We detail here the description for each of the parameters:

- **-r**  | **--reads**: Required. Reads mapped to junctions (file obtained through [Junckey](https://github.com/comprna/Junckey#1-format-star-output))

- **-o**  | **--output**: Required. output path.

- **-b**  | **--bam**: Required. Path to STAR output. ISOTOPE will make use of the mapped reads to do the empirical test of significance.

- **-trans**  | **--transcript**: Required. Transcript expression file.

- **-g**  | **--gtf**: Required. GTF annotation file. Only the exonic regions will be used.

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-genome**  | **--genome**: Required. Genome annotation file.

- **-mosea**  | **--mosea**: Required. [MoSEA](https://github.com/comprna/MoSEA) path

- **-mxfinder**  | **--mxfinder**: Required. [MxFinder](https://github.com/JLTrincado/MxFinder) path

- **-mut**  | **--mutations**: Mutations path file.

- **-m**  | **--max**: Maximum length for the exonizations. Default: 500nt.

- **-t**  | **--thres**: Minimum number of reads supporting the splicing event. Default: 5 reads.

- **-rand**  | **--rand**: Number of rounds for the empirical test. Default: 100.

- **--tumor_specific**: Obtains tumor-specific changes. Default: False

- **--control_path**: If --tumor_specific is activated, user can provide their own control samples in the same format as --reads file.

- **--Intropolis**: If --tumor_specific is activated, user can provide path to Intropolis file (https://figshare.com/s/ae2283696db044488c32)

- **--chess**: If --tumor_specific is activated, user can provide path to the CHESS file (https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe)

- **-rep**  | **--repeats**: Regions of the genome with repeats from maskerDB.

- **--temp**: Remove temporary files

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is executed on a slurm cluster. Default: False (it will run locally).

- **-user**  | **--user**: Only necessary if --cluster option is enabled. Cluster user name. 

- **-h**  | **--help**: display the help message describing the different parameters.

### Part3:

All the epitopes are evaluated to detect if these are actually neoepitopes derived from the splicing alteration (exonizations_NetMHC-4.0_neoantigens_type_gained.tab, exonizations_NetMHCpan-4.0_neoantigens_type_gained.tab) or potentially depleted as they correspond to self-antigens originating from the parts of the ORF in the reference transcripts that would be not present in the modified ORF after the splicing change (exonizations_NetMHC-4.0_neoantigens_type_lost.tab, exonizations_NetMHCpan-4.0_neoantigens_type_lost.tab).

```
python "${ISOTOPE_path}"/exonizations_ISOTOPE_part3.py --HLAclass /projects_rg/SCLC_cohorts/Smart/PHLAT/PHLAT_summary_ClassI.out --HLAtypes /projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan /projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC --netMHCpan /projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan -o "${output_path}"
```

We detail here the description for each of the parameters:

- **-o**  | **--output**: Required. output path

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-h**  | **--help**: display the help message describing the different parameters.

## Alternative splice site

Here we describe the steps to obtain alternative splice site events. The initial file is the same toy dataset ([readCounts_TEST.tab](https://github.com/comprna/ISOTOPE_supplementary/blob/main/readCounts_TEST.tab)) as before and most of the steps are pretty similar as with exonizations. It's composed of 2 parts.

### Part1:

In a similar way as with the exonizations, ISOTOPE will extract all the associated significant alternative splice site events. To assess the significance of each event recovered, ISOTOPE compares the reads mapping to each event against random junctions obtained on the same gene by empirical distribution. For all the events passing the previous filters, ISOTOPE will evaluate the likelihood to generate new open reading frames with respect to the canonical associated transcript (A5_A3_ORF.tab). Finally, the translated peptides are evaluated for their affinity for the MHC complex. Again this step is computationally intensive and has been parallelized (see option **--cluster**). When all jobs have finished, we can proceed with part2.

This time we deactivate --tumor_specific flag (we just do not include it on the command), because the A5_A3 events on our toy dataset are not tumor-specific. We include the Intropolis and CHESS paths if the user is interested to apply it on their data. CHESS files for A5_A3 events are provided ([chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe](https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe), [chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe](https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe)) as well as the Intropolis dataset (https://figshare.com/s/ae2283696db044488c32)

```
ISOTOPE_path="/users/genomics/juanluis/Software/ISOTOPE"
data_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE"
output_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE/A5_A3"

python "${ISOTOPE_path}"/A5_A3_ISOTOPE_part1.py -r "${data_path}"/data/readCounts_TEST.tab  --transcript "${data_path}"/data/iso_tpm.txt --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -conv "${data_path}"/annotation/Ensembl_gene_conversion.txt -genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa -mosea /users/genomics/juanluis/Software/MoSEA --mxfinder /users/genomics/juanluis/Software/MxFinder --rep "${data_path}"/annotation/hg19_repeats.bed -o "${output_path}" --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan --username juanluis --Intropolis "${data_path}"/annotation/intropolis.v1.hg19.filtered.tsv --chessA5 "${data_path}"/annotation/chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe --chessA3 "${data_path}"/annotation/chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe
```
We detail here the description for each of the parameters:

- **-r**  | **--reads**: Required. Reads mapped to junctions (file obtained through [Junckey](https://github.com/comprna/Junckey#1-format-star-output))

- **-o**  | **--output**: Required. output path.

- **-b**  | **--bam**: Required. Path to STAR output. ISOTOPE will make use of the mapped reads to do the empirical tests.

- **-trans**  | **--transcript**: Required. Transcript expression file.

- **-g**  | **--gtf**: Required. GTF annotation file. Only the exonic regions will be used.

- **-conv**  | **--conversion**: Required. File with the correspondence gene ID-gene symbol. We provide a [file](https://github.com/comprna/ISOTOPE_supplementary/blob/main/Ensembl_gene_conversion.txt) from biomart.

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-genome**  | **--genome**: Required. Genome annotation file.

- **-mosea**  | **--mosea**: Required. [MoSEA](https://github.com/comprna/MoSEA) path.

- **-mxfinder**  | **--mxfinder**: Required. [MxFinder](https://github.com/JLTrincado/MxFinder) path.

- **-mut**  | **--mutations**: Mutations path file.

- **-m**  | **--max**: Maximum length for the exonizations. Default: 500nt.

- **-t**  | **--thres**: Minimum number of reads mapping the event. Default: 5 reads.

- **-rand**  | **--rand**: Number of rounds for bootstrap testing. Default: 100.

- **--tumor_specific**: Obtains tumor-specific changes. Default: False.

- **--control_path**: If --tumor_specific is activated, user can provide their own control samples in the same format as --reads file.

- **--Intropolis**: If --tumor_specific is activated, user can provide path to Intropolis file (https://figshare.com/s/ae2283696db044488c32)

- **--chessA5**: If --tumor_specific is activated, the user must provide the path to the CHESS file with A5' events (https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe)

- **--chessA3**: If --tumor_specific is activated, the user must provide the path to the CHESS file with A3' events  (https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe)

- **-rep**  | **--repeats**: Regions of the genome with repeats from maskerDB

- **--temp**: Remove temporary files

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is executed on a slurm cluster. Default: False (it will run locally).

- **-h**  | **--help**: display the help message describing the different parameters

### Part2:

All the epitopes are evaluated to detect if these are actually neoepitopes gained because of the splicing alteration (A5_A3_NetMHC-4.0_neoantigens_type_gained.tab, A5_A3_NetMHCpan-4.0_neoantigens_type_gained.tab) or lost comparing against the self-antigens created by reference transcripts (A5_A3_NetMHC-4.0_neoantigens_type_lost.tab, A5_A3_NetMHCpan-4.0_neoantigens_type_lost.tab).

```
python "${ISOTOPE_path}"/A5_A3_ISOTOPE_part2.py --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan -o "${output_path}"
```

We detail here the description for each of the parameters:

- **-o**  | **--output**: Required. output path

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-h**  | **--help**: display the help message describing the different parameters

## Neoskipping

### Part1:

In a similar way as with the exonizations, ISOTOPE will extract all the neoskipping events. To assess the significance of each event recovered by --fold option only events with a greater fold of reads mapping the neoskipping with respect to the spanned junctions will be retained (Default: 1). For all the events passing the filters, ISOTOPE will evaluate the likelihood to generate new open reading frames with respect to the canonical associated transcript (all_neoskipping_filtered_peptide_change.tab). Finally, the translated peptides are evaluated for their affinity for the MHC complex. This step is computationally intensive and has been parallelized (see option **--cluster**). When all jobs have finished, we can proceed with part2:

```
ISOTOPE_path="/users/genomics/juanluis/Software/ISOTOPE"
data_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE"
output_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE/neoskipping"

python "${ISOTOPE_path}"/Neoskipping_ISOTOPE_part1.py -r "${data_path}"/data/readCounts_TEST.tab --transcript "${data_path}"/data/iso_tpm.txt --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa -mosea /users/genomics/juanluis/Software/MoSEA --mxfinder /users/genomics/juanluis/Software/MxFinder -o "${output_path}" --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan --temp -o "${output_path}" --Intropolis "${data_path}"/annotation/intropolis.v1.hg19.filtered.tsv --chess "${data_path}"/annotation/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe  
```
We detail here the description for each of the parameters:

- **-r**  | **--reads**: Required. Reads mapped to junctions (file obtained through [Junckey](https://github.com/comprna/Junckey#1-format-star-output))

- **-o**  | **--output**: Required. Output path.

- **-trans**  | **--transcript**: Required. Transcript expression file.

- **-g**  | **--gtf**: Required. GTF annotation file. Only the exonic regions will be used.

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-genome**  | **--genome**: Required. Genome annotation file.

- **-mosea**  | **--mosea**: Required. [MoSEA](https://github.com/comprna/MoSEA) path

- **-mxfinder**  | **--mxfinder**: Required. [MxFinder](https://github.com/JLTrincado/MxFinder) path

- **-mut**  | **--mutations**: Mutations path file.

- **-t**  | **--thres**: Minimum number of reads mapping the event. Default: 5 reads.

- **-f**  | **--fold**: Minimum fold of reads mapping the neoskipping with respect to the spanned junctions. Default: 0.

- **--tumor_specific**: Obtains tumor-specific changes. Default: False

- **--control_path**: If --tumor_specific is activated, user can provide their own control samples in the same format as --reads file.

- **--Intropolis**: If --tumor_specific is activated, user can provide path to Intropolis file (https://figshare.com/s/ae2283696db044488c32)

- **--chess**: If --tumor_specific is activated, the user must provide the path to the CHESS file (https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe)

- **--temp**: Remove temporary files

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is executed on a slurm cluster. Default: False (it will run locally).

- **-user**  | **--user**: Only necessary if --cluster option is enabled. Cluster user name.

- **-h**  | **--help**: display the help message describing the different parameters

### Part2:

All the epitopes are evaluated to detect if these are actually neoepitopes gained because of the splicing alteration (neoskipping_NetMHC-4.0_neoantigens_type_gained.tab, neoskipping_NetMHCpan-4.0_neoantigens_type_gained.tab) or lost comparing against the self-antigens created by reference transcripts (neoskipping_NetMHC-4.0_neoantigens_type_lost.tab, neoskipping_NetMHCpan-4.0_neoantigens_type_lost.tab).

```
python "${ISOTOPE_path}"/Neoskipping_ISOTOPE_part2.py --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan -o "${output_path}"
```

We detail here the description for each of the parameters:

- **-o**  | **--output**: Required. output path

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-h**  | **--help**: display the help message describing the different parameters

## Intron retention (IR)

For the analysis of IR changes the quantification of intron abundances needs to be computed. We provide a toy dataset for the tutorial [iso_tpm_introns_TEST.txt](https://github.com/comprna/ISOTOPE_supplementary/blob/main/iso_tpm_introns_TEST.txt). This file has been obtained using kma (https://github.com/pachterlab/kma). Each line is an intronic region and the columns are all the normalized abundances values for each of the samples (in this case TPMs, but any normalized expression could be valid).

In this part1, ISOTOPE will extract the tumor-specific introns by comparing them against the GTEX database (https://www.gtexportal.org/home/) if flag --tumor_specific is set. We have extracted the IR events from GTEX using SUPPA (https://github.com/comprna/SUPPA). We provide this file ([chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe](https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe)). Additionally, the user could provide some extra control files to compare with (--introns_normal). They should be in the same format as iso_tpm_introns_TEST.txt. Any IR event that is minimally expressed (by default 1, this value can be changed with the parameter --thres) in these control samples will be censored from our IR candidates.

For each of the surviving IR events from the previous step, we will compare its expression level with n random intronic positions to test significance. This step is computationally intensive and has been parallelized (see option **--cluster**). When all jobs have finished, we can proceed with part2: 

### Part1:

```
ISOTOPE_path="/users/genomics/juanluis/Software/ISOTOPE"
data_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE"
output_path="/users/genomics/juanluis/SCLC_cohorts/Smart/test_ISOTOPE/IR"

python "${ISOTOPE_path}"/IR_ISOTOPE_part1.py -i "${data_path}"/annotation/iso_tpm_introns_TEST.txt -b "${data_path}"/data/STAR --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -o "${output_path}" --tumor_specific --chess "${data_path}"/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe
```
- **-i**  | **--introns**: Required. Intron quantification file (obtained through [kma]((https://github.com/pachterlab/kma)))

- **-o**  | **--output**: Required. output path

- **-b**  | **--bam**: Required. Path to STAR output. ISOTOPE will make use of the mapped reads to do the empirical tests

- **-g**  | **--gtf**: Required. GTF annotation file. Only the exonic regions will be used.

- **-t**  | **--thres**: Minimum expression for the intron candidates. Default: 1.

- **-rand**  | **--rand**: Number of rounds for bootstrap testing. Default: 100.

- **--tumor_specific**: Obtains tumor-specific changes. Default: False

- **--control_path**: If --tumor_specific is activated, user can provide their own control samples in the same format as --introns file.

- **--chess**: If --tumor_specific is activated, the user must provide a path to the CHESS file (https://github.com/comprna/ISOTOPE_supplementary/blob/main/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe)

- **--temp**: Remove temporary files

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is executed on a slurm cluster. Default: False (it will run locally).

- **-h**  | **--help**: display the help message describing the different parameters

### Part2:

From all significant introns according to the previous test, ISOTOPE will evaluate if they would drive new open reading frames by comparing against the canonical associated transcript (IR_ORF_filtered_peptide_change.tab). For all those introns giving rise to a new open reading frame and not likely to be degraded by nonsense-mediated decay pathway, using NetMHC (https://services.healthtech.dtu.dk/service.php?NetMHC-4.0) and NetMHCpan (https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) we evaluate the likelihood that the translated peptides could generate epitopes bound by MHC complex. Again, this step is computationally intensive and has been parallelized (see option **--cluster**). When all jobs have finished, we can proceed with the last part3:

```
python "${ISOTOPE_path}"/IR_ISOTOPE_part2.py --transcript /projects_rg/SCLC_cohorts/Smart/Salmon/iso_tpm.txt --gtf "${data_path}"/annotation/Homo_sapiens.GRCh37.75.gtf -genome /genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa -mosea /users/genomics/juanluis/Software/MoSEA --mxfinder /users/genomics/juanluis/Software/MxFinder -o "${output_path}" --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan --username juanluis --tumor_specific
```

We detail here the description of the parameters:

- **-trans**  | **--transcript**: Required. Transcript expression file

- **-o**  | **--output**: Required. output path

- **-genome**  | **--genome**: Required. Genome annotation file.

- **-g**  | **--gtf**: Required. Gtf annotation file. Only the exonic regions will be used.

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-mosea**  | **--mosea**: Required. [MoSEA](https://github.com/comprna/MoSEA) path

- **-mxfinder**  | **--mxfinder**: Required. [MxFinder](https://github.com/JLTrincado/MxFinder) path

- **-t**  | **--thres**: Minimum expression to consider an intron valid. Default: 1.

- **--tumor_specific**: Obtains tumor specific changes. Default: False

- **--temp**: Remove temporary files

- **-c**  | **--cluster**: If activated, ISOTOPE will parallelize the execution if the software is executed on a slurm cluster. Default: False (it will run locally).

- **-user**  | **--user**: Only necessary if --cluster option is enabled. Cluster user name.

- **-h**  | **--help**: display the help message describing the different parameters


### Part3:

Part3 just summarises all calculations from part2. Epitopes are evaluated to detect if these are actually neoepitopes derived from a splicing alteration (IR_NetMHC-4.0_neoantigens_type_gained.tab, IR_NetMHCpan-4.0_neoantigens_type_gained.tab) or potentially depleted by comparing against the candidate self-antigens in the reference transcripts (IR_NetMHC-4.0_neoantigens_type_lost.tab, IR_NetMHCpan-4.0_neoantigens_type_lost.tab)

```
python "${ISOTOPE_path}"/IR_ISOTOPE_part3.py --HLAclass "${data_path}"/annotation/PHLAT_summary_ClassI.out --HLAtypes "${data_path}"/annotation/NetMHC-4.0_HLA_types_accepted.tab --HLAtypespan "${data_path}"/annotation/NetMHCpan-4.0_HLA_types_accepted.tab --netMHC /genomics/users/juanluis/Software/netMHC-4.0/netMHC --netMHCpan /genomics/users/juanluis/Software/netMHCpan-4.0/netMHCpan -o "${output_path}" 
```
We detail here the description of the parameters:

- **-o**  | **--output**: Required. output path

- **-HLAclass**  | **--HLAclass**: Required. HLA genotype of the samples.  

- **-HLAtypes**  | **--HLAtypes**: Required. HLA alelles recognized by NetMHC.  

- **-HLAtypespan**  | **--HLAtypespan**: Required. HLA alelles recognized by NetMHCpan.  

- **-netMHC**  | **--netMHC**: Required. NetMHC path.  

- **-netMHCpan**  | **--netMHCpan**: Required. NetMHCpan path.  

- **-h**  | **--help**: display the help message describing the different parameters

## Required packages:

- python3
- pandas
- biopython
- R (>=3.0.0)
- biomaRt
- bedtools
- MoSEA (https://github.com/comprna/MoSEA)
- MxFinder (https://github.com/JLTrincado/MxFinder)
