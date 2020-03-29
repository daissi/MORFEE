#' parse mutation
#'
#' @param x a string mutation: "-94G>A"
#'
#' @return a vector of different element of the mutation: c("-","94","G","A")
#'
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_mutation <- function(x){

  if(!is.character(x)){
    stop("x must be a string like '-94G>A'")
  }

  sig <- str_extract(x, "^-")
  if(is.na(sig)){
    sig <- "+"
  }
  pos <- str_extract(x, "^-*[0-9]+")
  mut <- gsub(pos, "", x)
  mut <- str_split_fixed(mut, ">", 2)
  pos <- gsub("-","",pos)

  c(sig,pos,mut)
}


#' parse refGene
#'
#' @param x a string refGene mutation: "NM_001010939:c.-94G>A"
#'
#' @return a vector of different element of the mutation: c("NM_001010939","c.-94G>A","-","94","G","A")
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_GeneDetail.refGene <- function(x){

  if(!is.character(x)){
    stop("x must be a string like 'NM_001010939:c.-94G>A'")
  }

  # Fix wrong character encoding!?!?!
  x <- gsub("\\\\x3b",";",x)

  x <- str_split(x, ";")[[1]]

  x <- str_split_fixed(x, ":", 2)

  y <- str_split_fixed(x[,2], "\\.", 2)[,2]
  z <- t(sapply(y, parse_mutation))

  a <- cbind(x,z)
  dimnames(a) <- NULL

  a
}

#' get data used by MORFEE
#'
#' @param path is a path to the temporary folder to store downloaded databases
#' @param force a boolean wheter the function should check if the database is already stored before to download them
#'
#' @importFrom utils download.file
#' @importFrom rtracklayer readGFF
#' @importFrom data.table fread
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
get.morfee.data <- function(path=tempdir(), force=FALSE){

  MORFEE.DATA <- list()

  MORFEE.DATA[["DB_PATH"]] <- paste0(path,"/MORFEE_DB/") # Path to the DB
  MORFEE.DATA[["DB_PATH_ANNOVAR"]] <- paste0(MORFEE.DATA[["DB_PATH"]],"ANNOVAR","/") # Path to the ANNOVAR DB
  MORFEE.DATA[["DB_PATH_GENCODE"]] <- paste0(MORFEE.DATA[["DB_PATH"]],"GENCODE","/") # Path to the GENCODE DB

  MORFEE.DATA[["GENCODE_URL"]] <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/"
  MORFEE.DATA[["GENCODE_FILE_ANNOT"]] <- "gencode.v33lift37.annotation.gff3.gz" # Position codon start + strand
  MORFEE.DATA[["GENCODE_FILE_METAD"]] <- "gencode.v33lift37.metadata.RefSeq.gz" # NM code
  MORFEE.DATA[["GENCODE_FILE_SEQUE"]] <- "GRCh37.primary_assembly.genome.fa.gz" # Whole Human Sequence

  MORFEE.DATA[["SEQ_INIT"]] <- Biostrings::DNAString("ATG") # START codon
  MORFEE.DATA[["SEQ_STOP"]] <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA")) # STOP codon

  if(!dir.exists(MORFEE.DATA[["DB_PATH_GENCODE"]])){
    dir.create(MORFEE.DATA[["DB_PATH_GENCODE"]] , recursive = TRUE)
  }

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]]))
  }

  MORFEE.DATA[["GENCODE_ANNOT"]] <- readGFF(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                   MORFEE.DATA[["GENCODE_FILE_ANNOT"]]) , version=3)

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_METAD"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_METAD"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_METAD"]]))
  }

  MORFEE.DATA[["GENCODE_METAD"]] <- fread(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                 MORFEE.DATA[["GENCODE_FILE_METAD"]]))

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]]))
  }

  MORFEE.DATA[["GENCODE_SEQ"]] <- readDNAStringSet(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                          MORFEE.DATA[["GENCODE_FILE_SEQUE"]]))

  MORFEE.DATA[["GENCODE_SEQ_ORDER"]] <- str_split_fixed(MORFEE.DATA[["GENCODE_SEQ"]]@ranges@NAMES, " ", n=2)[,1]

  return(MORFEE.DATA)
}

#' write a VCF object as an .ods OpenDocument file
#'
#' @param myvcf a VCF object
#' @param file name of the file to be created
#'
#' @importFrom readODS write_ods
#'
#' @export
#'
vcf_2_ods <- function(myvcf, file){

  df.myvcf <- combine.vcf.slot(myvcf)

  write_ods(x=df.myvcf, path=file)

}


#' write a VCF object as an .xlsx Office Open XML file
#'
#' @param myvcf a VCF object
#' @param file name of the file to be created
#'
#' @importFrom writexl write_xlsx
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
vcf_2_xlsx <- function(myvcf, file){

  df.myvcf <- combine.vcf.slot(myvcf)

  df.myvcf$orfSNVs <- unlist(get.orfSNVs(df.myvcf$MORFEE))
  df.myvcf$NewAALength <- unlist(get.NewAALength(df.myvcf$MORFEE))
  df.myvcf$Ratio_length_pred_obs <- get.Ratio_length_pred_obs(df.myvcf$NewAALength)

  col2keep <- c("seqnames","start","REF","ALT","Gene.refGene","avsnp150",
                "GeneDetail.refGene","orfSNVs","Ratio_length_pred_obs","NewAALength","MORFEE",
                "pLI.refGene","exp_lof.refGene","oe_lof.refGene","oe_lof_lower.refGene","oe_lof_upper.refGene",
                "gwasCatalog","CLNDN","CLNDISDB","CLNSIG","AF","AF_popmax",
                "SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred",
                "MutationAssessor_pred","FATHMM_pred","PROVEAN_pred","MetaSVM_pred","MetaLR_pred",
                "M.CAP_pred","REVEL_rankscore","MutPred_rankscore","CADD_phred","fathmm.MKL_coding_pred",
                "integrated_fitCons_score","GERP.._RS","phyloP100way_vertebrate_rankscore","phyloP20way_mammalian_rankscore",
                "phastCons100way_vertebrate_rankscore","phastCons20way_mammalian_rankscore",
                "SiPhy_29way_logOdds")# M-CAP_pred, fathmm-MKL_coding_pred, GERP++_RS

  df.ok <- as.data.frame(df.myvcf[,col2keep])

  for(i in 1:ncol(df.ok)){
    temp_aschar <- lapply(df.ok[,i], as.character)
    temp_paste <- lapply(temp_aschar, paste, sep="", collapse=",")
    temp_unlist <- unlist(temp_paste)
    df.ok[,i] <- gsub("\\\\x3b"," ; ",temp_unlist) # Fix encoding problem
    df.ok[,i] <- gsub("Name\\\\x3d","",temp_unlist) # Clean gwasCatalog
  }

  write_xlsx(x=df.ok, path=file)

}

get.Ratio_length_pred_obs <- function(x){

  z <- sapply(x, get.Ratio_length_pred_obs.meth)

  return(z)
}

get.Ratio_length_pred_obs.meth <- function(x){

  x <- str_split_fixed(x,";", n=Inf)
  y <- sapply(x, function(x) eval(parse(text=x)))
  y <- as.numeric(y)
  y.min.id <- which.min(abs(y-1))[1]
  z <- y[y.min.id]

  return(z)
}

get.NewAALength <- function(x){

  morfee_temp <-  lapply(x, function(x) paste0(x, sep="", collapse="|"))
  morfee_temp <-  lapply(morfee_temp, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])

  orf <- lapply(morfee_temp, get.NewAALength.meth)

  return(orf)
}

get.NewAALength.meth <- function(y){

  length_id <- seq(from = 4, to = length(y), by = 4)
  length_temp <- paste(y[length_id], collapse="; ", sep="")
  length_temp <- gsub("\\(aa\\)","",length_temp)
  length_temp <- gsub("\\[","",length_temp)
  length_temp <- gsub("\\]","",length_temp)
  return(length_temp)
}

get.orfSNVs <- function(x){

  morfee_temp <-  lapply(x, function(x) paste0(x, sep="", collapse="|"))
  morfee_temp <-  lapply(morfee_temp, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])

  orf <- lapply(morfee_temp, get.orfSNVs.meth)

  return(orf)
}

get.orfSNVs.meth <- function(y){

  orf_id <- seq(from = 3, to = length(y), by = 4)
  paste(y[orf_id], collapse="; ", sep="")
}


#' combine different slots of VCF object
#'
#' @param myvcf a VCF object
#'
#' @return a data.frame
#'
#' @importFrom DelayedArray rowRanges
#'
#' @export
#'
combine.vcf.slot <- function(myvcf){

  if(class(myvcf)!="CollapsedVCF"){
    stop("Class of myvcf must be 'CollapsedVCF'!!!")
  }

  df <- as.data.frame(rowRanges(myvcf))
  df.info <- as.data.frame(info(myvcf))
  df <- cbind(df,df.info)

  df.geno <- as.data.frame(geno(myvcf))
  if(nrow(df.geno)!=0){
    df <- cbind(df,df.geno)
  }

 return(df)
}

