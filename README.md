**MORFEE** (**Mutation on Open Reading FramE annotation**) is a tool
(R package) that, from a VCF file, detects and annotates single nucleotide
variants creating premature ATG codons.

MORFEE algorithm is written in R language and can run on all operating
systems which have an R interpreter (including Linux, macOS and Windows).
MORFEE starts with a minimal VCF file (i.e. with at least the *chr*, *position*,
*reference allele* and *alternate allele* fields) as an input, but can also
work with already [ANNOVAR](http://annovar.openbioinformatics.org/)-annotated
([Wang et al., 2010](https://doi.org/10.1093/nar/gkq603)) VCF files.
MORFEE has some R packages dependencies which are available on
[CRAN](https://cran.r-project.org/) or
[Bioconductor](https://www.bioconductor.org/) repositories and uses the
[GENCODE](https://www.gencodegenes.org/) database.

When using MORFEE, please cite our [publication](https://doi.org/10.1101/2020.03.29.012054).

## Installation
### From CRAN/Bioconductor
Not yet available!

### From GitHub
```
library("devtools")
install_github("daissi/MORFEE")
```

## Usage

In a first step, MORFEE reads the input VCF file and use ANNOVAR (that has then
to be beforehand installed) through the wrapper function
*vcf.annovar.annotation* to annotate all variants.
This step is skipped if the input file has already been annotated.
The minimal ANNOVAR annotations required by MORFEE are:
- **Func.refGene** to extract 5'UTR variants.
- **GeneDetail.refGene** to extract transcript ID.

### A. Annotation with ANNOVAR directly
```
# Annotation with ANNOVAR
table_annovar.pl \
        inst/extdata/VCF_example.vcf \
        ${PATH_ANNOVAR_DB} \
        -buildver hg19 \
        -remove -protocol refGene,gwasCatalog,avsnp150,clinvar_20190305,gnomad211_genome,dbnsfp35a \
        -operation gx,r,f,f,f,f -nastring . -vcfinput \
        -xref ${PATH_ANNOVAR_DB}/morfee.lof_metrics.by_gene.txt
```

### B. Annotation with ANNOVAR through MORFEE
```
library(VariantAnnotation)
library(MORFEE)

my_raw_vcf <- "inst/extdata/VCF_example.vcf"

# Annotation with ANNOVAR through the wrapper function
vcf.annovar.annotation(my_raw_vcf)
```
### Annotation with MORFEE
```
library("VariantAnnotation")
library("MORFEE")

# Download database and create an R object used by MORFEE
MORFEE_DATA <- get.morfee.data()

# To avoid too much Internet traffic, save locally this object
# save(MORFEE_DATA, file="MORFEE_DATA.RData")
# load("MORFEE_DATA.RData")

MY_VCF <- "VCF_example.hg19_multianno.vcf"

# Read/Load the VCF into R
my_vcf <- readVcf(MY_VCF)

# Perform the MORFEE annotation
my_vcf_morfee <- morfee.annotation(my_vcf, MORFEE_DATA)

# Keep only variants which create a new upstream ATG sequence
my_vcf_morfee_NEW <- my_vcf_morfee[!is.na(info(my_vcf_morfee)$MORFEE)]

# Write a VCF file with the new MORFEE field/annotation
writeVcf(my_vcf_morfee_NEW, filename="VCF_example.hg19_multianno.morfee_anno_NEW.vcf")

# Write an XLSX file with the new MORFEE field/annotation
vcf_2_xlsx(my_vcf_morfee_NEW, file="VCF_example.hg19_multianno.morfee_anno_NEW.xlsx")
```

Please see below the description of the different columns of the xlsx file:
- **GeneDetail.refGene**: Sequence Variant Nomenclature (ANNOVAR output).
       Please note, this nomenclature does not respect the [HGVS Recommendations](https://varnomen.hgvs.org/) ([doi:10.1002/humu.22981](https://doi.org/10.1002/humu.22981)).
       The position of the variant used with a prefix *c.* should be the position of the variant regarding the ATG position based one the coding sequence.
       In other words, the position on a sequence WITHOUT introns. ANNOVAR reports the position on the sequence WITH introns.
- **orfSNVs**: *in_frame* or *out_of_frame* to determine whether there is a frame shift between the new ATG and the reference one.
- **Ratio_length_pred_obs**: ratio between the predicted vs the observed length of the protein (<1 for a shorter protein and >1 for a longer protein). In case of several transcripts, only one value is given. That corresponding to the predicted length the closest to the observed length, i.e the ratio closest to 1.
- **NewAALength**: predicted length of the new protein in amino acids.
- **MORFEE**: MORFEE VCF-style annotation

*All fields below come from ANNOVAR annotation*:

- **pLI.refGene**: *probability of Loss-of-function Intolerance* score from [gnomAD](https://gnomad.broadinstitute.org/)
- **exp_lof.refGene**: expected *predicted Loss-of-Function* variants from gnomAD
- **oe_lof.refGene**: *observed / expected* score
- **oe_lof_lower.refGene**: lower limit of 90% confidence interval of *oe* metric
- **oe_lof_upper.refGene**: upper limit of 90% confidence interval of *oe* metric
- **gwasCatalog**: gene-associated phenotype from [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- **CLNDN**: disease name from [ClinVar](https://www.ncbi.nlm.nih.gov/variation/)
- **CLNDISDB**: database name and identifier for the disease name from ClinVar
- **CLNSIG**: clinical significance from ClinVar
- **AF**: allele frequency from [gnomAD](https://gnomad.broadinstitute.org/)
- **AF_popmax**: maximum population allele frequency from gnomAD
- **SIFT_pred**: *Sorting Tolerant From Intolerant* function prediction from [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)
- **Polyphen2_HDIV_pred**: Polyphen2 function prediction based on HumDiv from dbNSFP
- **Polyphen2_HVAR_pred**: Polyphen2 function prediction based on HumVar from dbNSFP
- **LRT_pred**: *Likelihood Ratio Test* function prediction from dbNSFP
- **MutationTaster_pred**: MutationTaster function prediction from dbNSFP
- **MutationAssessor_pred**: MutationAssessor function prediction from dbNSFP
- **FATHMM_pred**: *Functional Analysis through Hidden Markov Models* function prediction from dbNSFP
- **PROVEAN_pred**: *Protein Variation Effect Analyzer* function prediction from dbNSFP
- **MetaSVM_pred**: *Support Vector Machine* combination of function prediction score from dbNSFP
- **MetaLR_pred**: *Logistic Regression* combination of function prediction score from dbNSFP
- **M.CAP_pred**: *Mendelian Clinically Applicable Pathogenicity* likelihood score from dbNSFP
- **REVEL_rankscore**: *Rare Exome Variant Ensemble Learner* pathogenicity score from dbNSFP
- **MutPred_rankscore**: pathogenicity score from dbNSFP
- **CADD_phred**: *Combined Annotation Dependent Depletion* pathogenicity score from dbNSFP
- **fathmm.MKL_coding_pred**: fathmm-MKL function prediction based on HumDiv from dbNSFP
- **integrated_fitCons_score**: fitCons function prediction based on HumDiv from dbNSFP
- **GERP.._RS**: GERP++ RS conservation score from dbNSFP
- **phyloP100way_vertebrate_rankscore**: phyloP conservation score from dbNSFP
- **phyloP20way_mammalian_rankscore**: phyloP conservation score from dbNSFP
- **phastCons100way_vertebrate_rankscore**: phastCons conservation score from dbNSFP
- **phastCons20way_mammalian_rankscore**: conservation score from dbNSFP
- **SiPhy_29way_logOdds**: SiPhy conservation score from dbNSFP
