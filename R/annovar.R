#' check if ANNOVAR is installed
#' otherwise open the website of ANNOVAR
#'
#' @param path_to_annovar define path to ANNOVAR binaries
#' @param open_browser booleen if TRUE, when ANNOVAR is not detected a web-browser is open to the ANNOVAR's website
#'
#' @importFrom utils browseURL
#'
#' @export
#'
is.annovar.installed <- function(path_to_annovar=NULL, open_browser=FALSE){
  message("Detection of ANNOVAR...")

  test.annovar <- ((Sys.which("annotate_variation.pl") == "") | (Sys.which("table_annovar.pl") == ""))

  if(test.annovar){
    message("Please install ANNOVAR (annotate_variation and table_annovar) in your PATH!")

    if(open_browser){
      browseURL("http://annovar.openbioinformatics.org/")
    }

    return(FALSE)
  }else{
    message("ANNOVAR is detected in your PATH!")
    return(TRUE)
  }
}

#' download ANNOVAR DB for annotation
#' https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/#additional-databases
#'
#' @param force a boolean wheter the function should check if the database is already stored before to download them
#' @param path_annovar_db a path where the database will be stored
#'
#' @importFrom R.utils gunzip
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
#'
download.annovar.db <- function(force=FALSE, path_annovar_db="DB_annovar"){

  if(!is.annovar.installed()){
    stop("ANNOVAR is not installed!")
  }


  if(!dir.exists(path_annovar_db)){
    dir.create(path_annovar_db)
  }

  ## RefGene
  if(!file.exists(paste(path_annovar_db,"hg19_refGene.txt", sep="/"))){

    my_annovar_cmd_A <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb -webfrom annovar refGene",
                                path_annovar_db)

    system(my_annovar_cmd_A)
  }

  ## dbSNP 150
  if(!file.exists(paste(path_annovar_db,"hg19_avsnp150.txt", sep="/"))){

    my_annovar_cmd_B <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb -webfrom annovar avsnp150",
                                path_annovar_db)

    system(my_annovar_cmd_B)
  }

  ## gnomAD 2.1.1
  if(!file.exists(paste(path_annovar_db,"hg19_gnomad211_genome.txt", sep="/"))){

    my_annovar_cmd_C <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb -webfrom annovar gnomad211_genome",
                                path_annovar_db)

    system(my_annovar_cmd_C)
  }

  ## clinvar 20190305
  if(!file.exists(paste(path_annovar_db,"hg19_clinvar_20190305.txt", sep="/"))){

    my_annovar_cmd_D <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb -webfrom annovar clinvar_20190305",
                                path_annovar_db)

    system(my_annovar_cmd_D)
  }

  ## dbNSFP 3.5a (for PhyloP score)
  if(!file.exists(paste(path_annovar_db,"hg19_dbnsfp35a.txt", sep="/"))){

    my_annovar_cmd_D <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb -webfrom annovar dbnsfp35a",
                                path_annovar_db)

    system(my_annovar_cmd_D)
  }

  ## gwasCatalog
  if(!file.exists(paste(path_annovar_db,"hg19_gwasCatalog.txt", sep="/"))){

    my_annovar_cmd_E <- paste("annotate_variation.pl",
                                "-buildver hg19",
                                "-downdb gwasCatalog",
                                path_annovar_db)

    system(my_annovar_cmd_E)
  }

  ## Constraint pLoF Metrics (gnomAD)
  if(!file.exists(paste(path_annovar_db,"morfee.lof_metrics.by_gene.txt", sep="/"))){

    download.file("https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
                  paste(path_annovar_db,"gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", sep="/"))

    gunzip(paste(path_annovar_db,"gnomad.v2.1.1.lof_metrics.by_gene.txt.gz", sep="/"))

    lof_gene <- fread(paste(path_annovar_db,"gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="/"))

    colnames(lof_gene)[colnames(lof_gene)=="gene"] <- "#Gene_name" # Rename first column to allow ANNOVAR to detect the header

    col2keep <- c("#Gene_name","exp_lof","pLI","oe_lof","oe_lof_lower","oe_lof_upper")

    fwrite(lof_gene[, col2keep, with = FALSE], sep="\t", na = "NA", quote = FALSE,
           file=paste(path_annovar_db,"morfee.lof_metrics.by_gene.txt", sep="/"))

  }

}


#' vcf annotation with ANNOVAR
#'
#' @param my_raw_vcf a string to define name of raw VCF to annotate
#' @param path_results a path where the annoated VCF will be saved
#' @param path_annovar_db a path where the database will be stored
#'
#' @export
#'
vcf.annovar.annotation <- function(my_raw_vcf, path_results="MORFEE_results", path_annovar_db="DB_annovar"){

  if(!is.annovar.installed()){
    stop()
  }

  if(!dir.exists(path_annovar_db)){
    dir.create(path_annovar_db)
  }

  if(!dir.exists(path_annovar_db)){
    dir.create(path_annovar_db)
  }

  my_annovar_cmd_Z <- paste("table_annovar.pl", my_raw_vcf, paste0(path_annovar_db,"/"),
                              "-buildver hg19",
                              "-out", paste(path_results,gsub(".vcf", "", basename(my_raw_vcf)), sep="/"),
                              "-remove -protocol refGene,gwasCatalog,avsnp150,clinvar_20190305,gnomad211_genome,dbnsfp35a -operation gx,r,f,f,f,f -nastring . -vcfinput",
                              "-xref",paste(path_annovar_db,"morfee.lof_metrics.by_gene.txt", sep="/"))

  system(my_annovar_cmd_Z)

}

