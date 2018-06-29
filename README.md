# Reproduce results for cancer eQTL mapping study

Note: There are additional annotation files required before you run the scripts below, which can be obtained here: https://osf.io/e467a/?action=download There are files for SNP array annotation, tumor purity estimates, TCGA clinical annotations, or gene annotation. These files should all be copied to "theRootDir", before the scripts below are run.

Note: These scripts need to be run in the sequence below: the early scripts create / download data that is used by the subsequent scripts!

Note: Substantial computing resources are required to run some of these scripts. The scripts were tested on the Bionimbus Protected Data Cloud, using Ubuntu Linux 14.04.3 LTS, under a configuration with 32 cores and 128Gb or RAM. R version 3.2.2 was used. 


# Instructions to install required packages and setup environment

Please install the following packges. From the R prompt:

```
> source("http://bioconductor.org/biocLite.R")
> biocLite(c("parallel", "grid", "scales", "ggthemes", "pd.genomewidesnp.6", "GO.db", "RColorBrewer", 
"reshape2", "goseq", "PEER", "car", "GenomicRanges", "GenomicFeatures", "ggplot2", "org.Hs.eg.db", 
"TxDb.Hsapiens.UCSC.hg19.knownGene", "glmnet", "gdata", "knitr", "preprocessCore", "Hmisc", "qqman", "MASS"))
```

# Scripts to download data and reproduce analysis

**theRootDir.R**
This file defines the variable "theRootDir", which is used as the root directory for all other scripts, which build various directory structures on top of this. This is assumed to be located at "/mnt/data_scratch/prediXcanProj/". This should either be updated to match the location you have used on your machine, or you must create this directory.

**simulatedData.R**
This script creates the simulated data and all of the result that pertain to the simulated data, i.e. Figure 1 and all associated results and tables.

**read_raw_tcga_breast_cancer_snp_data.R**
This script processes the raw Affymetrix SNP 6.0 array data that was downloaded from GDC (instructions on downloading the raw (protected) SNP data from GDC are contained in the script), however authorization must be obtained from TCGA to access these germline genetic data and we are NOT permitted to re-post these data. The script outputs the genotypes in a matrix of with row names rsIDs and column names TCGA sample name, which is saved as an RData file and used by subsequent scripts. This script also relies on a number of Affymetrix annotation files that we provide here. Details of where these files came from (and how to obtain them yourself) are contained in the script.

**download_tcga_expression_data.R**
A script to download the TCGA gene expression data from FireBrowse.

**batch_correct_tcga_data.R**
Script to process TCGA expression data and batch correct it using the RUV approach. We described this approach in a previous publication (see Methods section of Geeleher et al., Genome Research 2017, PMID:28847918))

**tcga_expression_and_snp_data_matrixEqtl_format.R**
Script to read the TCGA SNP and gene expression data (from the file formats output in the previous scripts) and output them in Matrix eQTL (PMID: 22492648) format . Note, these Matrix eQTL format files are used in subsequent scripts, but Matrix eQTL was not used for the analysis in the paper. Note this file also includes tumor purity estimates from the following paper (Onuchic et al., Cell Reports 2016, PMID:27851969); these tumor purity estimates were not used in the final analysis, but they are used in some analysis as a comparison of methods for estimating tumor purity. These tumor purity estimates are available with the supplementary data of Onuchic et al., and we have also made them available here.

**filter_tcga_eqtl_data.R**
This script filters the TCGA expression and SNP data based on the guidelines of GTEx (see script for details). It processes and outputs the CPE tumor purity estimates used in the main analysis. CPE purity estimates are available in the supplementary data of Aran et al., Nature Communications 2017 (PMID:26634437); we also provide them. This script also uses some TCGA clinical data (to filter males from the breast cancer dataset), provided here. The script saves an RData file (containing subsetted and ordered SNP, expression, purity data etc.) that is used in subsequent eQTL analysis by the next scripts.

**run_cis_eqtl_analysis_in_r.R**
Runs an eQTL analysis in R. Written in a way to be reasonably efficient (by grouping SNPs by gene only one matrix inversion per gene, by far the most computationally intensive step in model fitting). Script takes approx 40hrs to run on a standard desktop (no multi-threading). This script performs various eQTL mappings, for example, with an without PEER factors and encoding tumor purity in different ways. It also calculates the PEER factors. The output is processed by the next script (prepare_results_allAnalyses.R).

**prepare_results_allAnalyses.R**
This script takes the output of the script above (run_cis_eqtl_analysis_in_r.R) and processes the results, outputting them in a usable format. I.e. outputted model coefficients and creating matrices of P-values, effect sizes and standard errors, which this script takes care of in an efficient way.

**compareVariousEqtls.R**
This script inputs the data prepared above (i.e. the eQTL data) and outputs several of the results and figures (Fig. 2(d-f)) included in the main text. For example the numbers of significant eQTLs with the various approaches.

**cpeExpressionCorsFig2.R**
Creates Fig. 2(a-c) and associated results. i.e. the results pertaining to the estimated tumor purities and how those are related to gene expression.

**gtex_and_tcga_integrative_analysis.R**
This script does the integrative analysis of GTEx and TCGA data, everything in Figure 3 and associated results and tables.
Note: This script relies on summary statistics that were downloaded from GTEx, specifically, the release GTEx v6 data. The following files must be downloaded from GTEx and place in your "data/" directory (i.e. theRootDir/data, which should have already been created by previous scripts):

Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt
Cells_EBV-transformed_lymphocytes_Analysis.v6p.all_snpgene_pairs.txt
Cells_Transformed_fibroblasts_Analysis.v6p.all_snpgene_pairs.txt
Breast_Mammary_Tissue_Analysis.v6p.all_snpgene_pairs.txt

Furthermore, some annotation files are required to convert the GTEx SNP IDs to RS IDs etc. We have included these with the other annotation data here. These are:

gencode.v19.genes.v6p.patched_contigs.parsed.txt
GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt

One other key step in this script is flipping the direction-of-effect for some SNPs; this is because TCGA and GTEx have not always used the same reference allele. Furthermore, the annotation files provided by Affymetrix, and the annotations used by GTEx, do not always report the SNP nucleotide info (G,C,T,A) on the same strand, i.e. GTEx have always reported this information from the forward strand, but the Affymetrix annotation files used to process the TCGA data are mixed between the forward and reverse strands; this is defined in this file (originally part of the Birdsuite Metadata and included with the other annotation files we provide here):

affy6arrayStrandInfo.csv
Finally, the GTEx summary data in particular is very large and requires a lot of memory to load. Tested on a machine with 128Gb or RAM.

**imputed_snps_with_gtex_analysis.R**
This script does the analysis with the imputed SNPs from the large breast cancer GWAS meta-analysis and integrative analysis of GTEx, TCGA and GWAS data, i.e. everything in Figure 4 and associated results and supplementary figures/tables. Note, this script requires the imputed TCGA germline genotype data, which is protected data, that we are not permitted to repost (contained in brcaGwasSigSnpsImputed.RData). NB: Authorization must be obtained from TCGA to access all germline genotype data.
Additional Script:

**processImputedVcfs_TcgaBrca.R**
Additionally, we provide a script we used to process the imputed TCGA SNP data, outputted as VCF files by the Michigan Imputation Server. Details of how these imputations were performed are complex, and are included in the Methods section of our paper. These files are enormous and we are not permitted to re-post them, as these are protected data, which can be used to de-identify patients, and thus authorization must be obtained from TCGA in order to access them. We provide the output of this script as an RData file, as it is necessary to run the "imputed_snps_with_gtex_analysis.R" script. For completeness, we also provide this script, which was used to process these files. 

### TCGA eQTL analysis by PAM50 subtype
**run_cis_eqtl_analysis_bySubtype.R** #' This script will run conventaional and interaction models, stratified by PAM50 defined subtypes. We have enough samples for  Basal-like, HER2-enriched, Luminal A, Luminal B. The output of this script is processed by prepare_results_allAnalyses.R.

**compare_results_between_subtypes.R** Compares the results between different subtypes and against the pooled analysis. This script also creates the relevant supplementary tables (4 & 5)

### TCGA analysis controlling for copy number
**controlForCnv.R** Runs the eQTL analysis in TCGA controlling for copy number. The output of this script is processed by prepare_results_allAnalyses.R. Note, gene level copy number was calculated using the script "map_cnvs_to_genes_for_eQTL_proj.R"

### TCGA analysis additionally controlling for methylation
**controlForMethylation.R** Runs the eQTL analysis in TCGA controlling for methylation.
**compare_results_with_without_cnvsAndMethy.R** Compares the output when we have controlled for CNVs, methylation, or only PEER. Creates the supplementary figure for this analysis.

### TCGA analysis for only 10% of samples with the highest tumor purity
**analysis_on_top_10pct.R** Runs the analysis, the output of which is processed by prepare_results_allAnalyses.R
**interpret_results_for_top_10pct.R** Compares these results to the analyses on the entire TCGA dataset.

### METABRIC analysis and comparison to TCGA
**processMetabric_genotypes.R** Reads the METABRIC genotype files output by Affymetrix Genotyping Console and converts them into a usable format, i.e. changes the genotype calls to a 0,1,2 numeric matrix.
**processMetabric.R** Reads the genotype and expression data, processes and filters them, and outputs the covaraiates for expression heterogeneity, purity and anscestry.
**doMetaBricEqtlAnalysis.R** Runs the actual METABRIC eQTL association analysis.
**compare_results_against_metabric.R** Compares the METABRIC results to the TCGA results and outputs the results.




