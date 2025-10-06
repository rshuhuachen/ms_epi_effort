# Temporal epigenetic changes

This is the workflow used for the manuscript titled "Epigenetic changes associated with reproductive investment and life-history trade-offs in lekking male black grouse (Lyrurus tetrix)" (in review) which can be read on EcoEvoRxiv with DOI [10.32942/X2506F](https://ecoevorxiv.org/repository/view/10410/).

## Overview

This directory consists of four main sub-directories: data, scripts, results and plots. This project uses 118 samples that are part of a larger project where we sequenced 450 epigenomes of black grouse. We used epigbs, a reduced representation bisulphite sequencing method, to quantify DNA methylation at CpG sites. The epigbs method includes both a lab protocol executed at the Netherlands Institute of Ecology, as well as a bioinformatic pipeline. The forked pipeline can be found on [github](https://github.com/rshuhuachen/ms_epigbs_grouse/tree/main), which has been forked to be modified for the black grouse data. The current directory only contains the analyses from the point where the epigbs pipeline has been completed, meaning that the pooled libraries that contain 22 samples each were demultiplexed and CpG methylation was called. The samples used in the current study are divided over all 21 sequenced libraries. The metadata for this can also be found in the current project.

You can clone this directory by executing the following command:

`git clone https://github.com/rshuhuachen/ms_epi_effort.git`

All analyses were executed in R using RStudio on a HPC as many models are run simultaneously, and the epigenetic data are large in size. All R package versions can be found in the main manuscript.

## Data

The raw data used in this manuscript are not included in this github directory, but need to be downloaded from NCBI. The birds included in this dataset have also been genotyped using WGS to quantify inbreeding and mutation load, see the publication on [mutation load](https://www.nature.com/articles/s41559-025-02802-8) and the corresponding [github directory](https://github.com/rshuhuachen/ms_load_grouse/).

### Reference genome

The reference genome can be found on NCBI BioProject [PRJNA1085187](https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJNA1085187). This file is not required for analysis, but is required for the epigbs pipeline.

### Annotation

The genome annotation can be found in [this repository](https://github.com/rshuhuachen/ms_load_grouse/) in the folder `data/genomic/annotation/` under name PO2979_Lyrurus_tetrix_black_grouse.annotation.gff. The RNA seq data generated to annotate the genome can be found under NCBI Bio Accession [SRX24353984](https://www.ncbi.nlm.nih.gov/sra/SRX24353984%5Baccn%5D%5D). There are files and scripts that were produced to extract e.g. promoter regions, which can be found [here](https://github.com/rshuhuachen/grouse-annotation.git).

### EpiGBS data

The raw sequencing data can be found at NCBI BioProject PRJNA1085187 with Bio Accessions [SAMN51757821](https://www.ncbi.nlm.nih.gov/biosample/51757821) - [SAMN51757844](https://www.ncbi.nlm.nih.gov/biosample/51757844). Note that these raw sequences are pooled libraries that contain 22 barcoded samples each. Demultiplexing is done within the epigbs bioinformatic pipeline (see above). I also copied the barcodes used per library in the `data/barcodes_epigbs` subdirectory which were used for demultiplexing.

### Phenotypic and metadata

Phenotypic and metadata can be found in the subdirectory ./`data.` This folder contains one subdirectory: and two other files.

The `metadata_samples.csv`file includes the following column names: id (bird ID), site (lek site), age_cat (yearling or adult age category), prepost (whether the sample was taken pre or post lekking), fulldate (the full date when the sample was taken), attend (lek attendance), dist (lek centrality measured as distance to the lek centre), MS (mating success), surv (survival to the next autumn), blue_nextyear (blue chroma expressed in the next lekking season), lyre_nextyear (lyre size expressed in the next lekking season). This is the main data file.

The `scaffold_names_dovetail.RData` file contains the names of a scaffold, the size of each scaffold in bases and bytes, how many lines of bases and bites, and the assigned scaffold number based on its length (scaffold 1 is the largest).

The `phenotypes` folder contains all raw phenotypic data for all epigbs samples (all 450, not all used for the current analysis). Many of these contain additional phenotypic data not used for the current analysis, and were used to pre-process the final phenotypic file.

The raw methylation dataset (2GB in size) which has been pre-processed with methylKit to merge strands (i.e. a critical intermediate file that can be used for analysis without executing the entire epigbs pipeline) can be downloaded from Figshare with DOI [10.6084/m9.figshare.30186571](https://figshare.com/articles/dataset/Raw_methylkit_object_for_pre-lekking_and_post-lekking_samples_from_black_grouse_males/30186571) and should be deposited in `./data/processed`.

### Whole genome bisulphite sequencing

The pipeline executed for the single WGBS sample can be found here <https://github.com/rshuhuachen/wgbs_grouse>. The data used can be found under the same NCBI BioProject with SRA BioAccession Number [SAMN51785344](https://www.ncbi.nlm.nih.gov/biosample/51785344). The processed intermediate datafile to produce the figure (Supplementary Figure 1) can be found in `results/wgbs/`

## Scripts

The R scripts can be found in the `scripts` directory. They are numbered according to the order in which they should be executed.

1_prepare_data. Here, we load in some raw phenotypic data files and process it to isolate the data that were used for the current paper. Only phenotypic/metadata is involved in this file.

2_filter_methdata. Here, we load in the methylation data produced by the epigbs bioinformatic pipeline (bismark, to be exact), do some exploration, and filter it to produce our high-quality CpG site dataset. The processing is mostly done with the R package methylkit. Within the file, we also call a custom script stored in `scripts/function_convert_methfile`which reorganises the datafile to make it easier to analyse the data.

3_explore. Here, we get some summary statistics from the high quality datafile, compute a principal component analysis (PCA) on the raw data. We also test here with linear models if library had an effect on methylation %

4_calculate_delta. Here, we use the processed methylation file to calculate delta methylation.

5_changing_meth. Here, we identify CpG sites that significantly change in methylation across the lekking season, i.e. dynamic CpG sites. We make a function which can be used iteratively to compute a model per CpG site and run these models in parallel. We explore overdispersion and exclude CpG sites where the model cannot run successfully. We also explore some of the raw data and model output, and annotate the CpG sites according to the gene annotation file. Lastly, we execute several binomial tests to see if dynamic CpG sites are located in certain genes more often than expected.

6_effort_ms_meth. Here, we model how variation in reproductive effort is associated with methylation differences. Within this script, we call another script `./scripts/function_models` where we have specified a few functions that make iteration of models across different traits easier. We randomly select one pair of samples for those individuals who have multiple samples, as we do not have enough repeated samples to include ID as a random effect without causing model fititng issues.

7_nextyear_meth and 8_survival_meth. In these scripts, we do something similar to script 6 but instead, we model how delta methylation is associated with future traits: survival to the next year and plumage expression in the next year.

9_all_almost_sig_cpg. In this last script, we annotate all sites that were significant (before multiple testing correction) and we execute a few binomial tests to see if certain genes contain more than expeced significant CpG sites

plots. Within the `./scripts/plots` you will find the 4 scripts used to create the main and supplementary figures, and a script used to make qqplots.

## Plots

Within the ./plots subdirectory, you will find a few folders. All created plots in scripts 1 to 9 can be found here, and the final plots used for the manuscript are stored in the 'final' folder.

## Results

In this directory, most results are outputted but not all can be found on github due to filesize restrictions. The file that are synced are the outputs from GO analyses and gene lists of significant CpG sites.

## 
