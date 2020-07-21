### -------------------------------------------------------------------------
### MorfeeVCF 
###

setClass("MorfeeVCF", 
    contains="CollapsedVCF",
    prototype=prototype(
      fixed=DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList(),
                      QUAL=numeric(), FILTER=character())) )

### Automatically generated "coerce<-" method is broken so we fix it.
### See S4Vectors/R/S4-utils.R in the S4Vectors package for more information.
S4Vectors:::setReplaceAs("MorfeeVCF", "RangedSummarizedExperiment",
    S4Vectors:::canonical_replace_as_2
)
