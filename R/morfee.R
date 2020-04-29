#' vcf annotation with MORFEE
#'
#' @param myvcf_annot an ANNOVAR annotated VCF object to annotate with MORFEE
#' @param morfee_data data object obtained from get.morfee.data()
#'
#' @importFrom stats na.omit
#' @import foreach
#' @import VariantAnnotation
#' @import Biostrings
#' @import GenomicRanges
#'
#' @export
#'
morfee.annotation <- function(myvcf_annot, morfee_data){

  myvcf_annot_info <- info(myvcf_annot)
  myvcf_annot_info$MORFEE <- NA
  myvcf_annot_header <- rbind(info(header(myvcf_annot)),
                              data.frame(Number = ".", Type = "String", Description = "New ATG annotation provided by MORFEE", stringsAsFactors = FALSE) )

  rownames(myvcf_annot_header) <- c(rownames(info(header(myvcf_annot))),"MORFEE")
  info(header(myvcf_annot)) <- myvcf_annot_header
  info(myvcf_annot) <- myvcf_annot_info

  i <- NULL
  for(i in 1:nrow(myvcf_annot_info)){

    my_func <- as.character(myvcf_annot_info[i,"Func.refGene"])
    if(my_func!="UTR5"){
      message("Skip a variant which is not in UTR5 region")
      next
    }

    my_refgene <- as.character(myvcf_annot_info[i,"GeneDetail.refGene"])

    if(my_refgene=="."){
      message("Skip a variant which has no GeneDetail.refGene annotation")
      next
    }

    my_gene <- as.character(myvcf_annot_info[i,"Gene.refGene"])
    my_nm_list <- parse_GeneDetail.refGene(my_refgene)

    my_snp_pos_geno <- start(ranges(rowRanges(myvcf_annot)))[i]
    my_chr <- paste0("chr",as.character(seqnames(rowRanges(myvcf_annot))[i]))
    my_gencode_seque <- morfee_data[["GENCODE_SEQ"]][which(morfee_data[["GENCODE_SEQ_ORDER"]]==my_chr)]

    # Loop for each transcript (row)
    for(nm in 1:nrow(my_nm_list)){
      my_nm <- my_nm_list[nm,1]
      my_snp <- my_nm_list[nm,2]
      my_upordown <- my_nm_list[nm,3]
      my_snp_pos_rel <- as.numeric(my_nm_list[nm,4]) # Position from ATG, ex: 94
      my_nm_id <- grep(paste0(my_nm,"\\."), morfee_data[["GENCODE_METAD"]]$V2)

      if(nchar(my_nm_list[nm,5])!=1){
        message(" indel detected! Not YET supported!")
        next
      }
      if(my_nm_list[nm,5]=="-"){
        message(" indel detected! Not YET supported!")
        next
      }
      if(nchar(my_nm_list[nm,6])!=1){
        message(" indel detected! Not YET supported!")
        next
      }
      if(my_nm_list[nm,6]=="-"){
        message(" indel detected! Not YET supported!")
        next
      }

      if(length(my_nm_id)==0){
        message("Annotation for",my_nm_id," was not found in GENCODE !?!")
        message(" skipping this variant for now... NOT YET IMPLEMENTED")
        next
      }

      my_enst <- morfee_data[["GENCODE_METAD"]][my_nm_id,1]
      my_transcript_id <- grep(paste0(my_enst,"_"), morfee_data[["GENCODE_ANNOT"]]$transcript_id)
      gencode_annot_sub <- morfee_data[["GENCODE_ANNOT"]][my_transcript_id,]
      gencode_annot_cds <- gencode_annot_sub[gencode_annot_sub$type=="CDS",]
      gencode_annot_exon <- gencode_annot_sub[gencode_annot_sub$type=="exon",]

      my_init_codon_r <- gencode_annot_sub[gencode_annot_sub$type=="start_codon",]

      if(nrow(my_init_codon_r)==0){
        message("Variant without start_codon! Skip")
        next
      }

      my_init_codon_start <- as.numeric(my_init_codon_r[,"start"])[1]  # Several start codon?
      my_init_codon_end   <- as.numeric(my_init_codon_r[,"end"])[1]  # Several start codon?

      my_stop_codon_r <- gencode_annot_sub[gencode_annot_sub$type=="stop_codon",]

      if(nrow(my_stop_codon_r)==0){
        message("Variant without stop_codon! Skip")
        next
      }

      my_stop_codon_start <- as.numeric(my_stop_codon_r[,"start"])[1]  # Several stop codon?
      my_stop_codon_end   <- as.numeric(my_stop_codon_r[,"end"])[1]  # Several stop codon?

      if(my_init_codon_end < my_stop_codon_end){
        my_init_codon_5 <- my_init_codon_start
        my_init_codon_3 <- my_init_codon_end
        my_stop_codon_5 <- my_stop_codon_start
        my_stop_codon_3 <- my_stop_codon_end
      }else{
        my_init_codon_5 <- my_init_codon_end
        my_init_codon_3 <- my_init_codon_start
        my_stop_codon_5 <- my_stop_codon_end
        my_stop_codon_3 <- my_stop_codon_start
      }

      ########################################
      # Test gene orientation
      if(my_init_codon_end < my_stop_codon_end){

        gencode_annot_exon <- gencode_annot_exon[order(as.numeric(gencode_annot_exon$exon_number), decreasing = FALSE),]
        exons_length <- (gencode_annot_exon[,"end"]+1)-gencode_annot_exon[,"start"]
        my_snp_exon <- which(gencode_annot_exon[,"start"] <= my_snp_pos_geno & gencode_annot_exon[,"end"] >= my_snp_pos_geno)
        my_start_exon <- which(gencode_annot_exon[,"start"] <= my_init_codon_5 & gencode_annot_exon[,"end"] >= my_init_codon_5)
        my_cdna_length_A <- my_cdna_length_B <- 0

        if(length(my_snp_exon)==0){
          message("Variant not in an exonic part. Skip")
          next
        }

        # Reference sequence stats
        for(exon_i in 1:nrow(gencode_annot_exon)){
          if(exon_i==1){
            my_cdna   <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
          }else{
            my_cdna_i <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
            my_cdna <- c(my_cdna, my_cdna_i)
          }

          # Calcul position of variant in cDNA
          if(exon_i < my_snp_exon){
            my_cdna_length_A <- my_cdna_length_A + exons_length[exon_i]
          }else if(exon_i == my_snp_exon){
            my_snp_pos_cdna <- my_cdna_length_A + ((my_snp_pos_geno+1) - gencode_annot_exon[exon_i,"start"])
          }

          # Calcul position of reference A(TG) in cDNA
          if(exon_i < my_start_exon){
            my_cdna_length_B <- my_cdna_length_B + exons_length[exon_i]
          }else if(exon_i == my_start_exon){
            my_init_codon_5_cdna <- my_cdna_length_B + ((my_init_codon_5 + 1) - gencode_annot_exon[exon_i,"start"])
          }
        }

        # Find codon start in reference sequence
        stats_orig <- matchPattern(morfee_data[["SEQ_INIT"]], my_cdna)


          # Found all STOP in reference sequence
          for(j in 1:length(morfee_data[["SEQ_STOP"]])){

            stats_stop_orig_j <- matchPattern(morfee_data[["SEQ_STOP"]][[j]], my_cdna)

            if(j==1){
              stats_stop_orig <- stats_stop_orig_j
            }else{
              stats_stop_orig <- c(stats_stop_orig, stats_stop_orig_j)
            }
          }


        my_ref_allele <- as.character(subseq(my_cdna, start=my_snp_pos_cdna, end=my_snp_pos_cdna))
        if(my_ref_allele!=my_nm_list[nm,5]){
          message("Mismatch between alleles")
          next
        }

        my_ref_ATG <- as.character(subseq(my_cdna, start=my_init_codon_5_cdna, end=(my_init_codon_5_cdna+2)))
        if(my_ref_ATG!="ATG"){
          message("Reference ATG no detected!")
          next
        }

        # Replace my sequence with SNP
        my_cdna_updated <- replaceLetterAt(my_cdna, my_snp_pos_cdna, my_nm_list[nm,6])

        ########################################
        }else{
        # message("Opposite orientation!")

          gencode_annot_exon <- gencode_annot_exon[order(as.numeric(gencode_annot_exon$exon_number), decreasing = TRUE),]
          exons_length <- (gencode_annot_exon[,"end"]+1)-gencode_annot_exon[,"start"]
          my_snp_exon <- which(gencode_annot_exon[,"start"] <= my_snp_pos_geno & gencode_annot_exon[,"end"] >= my_snp_pos_geno)
          my_start_exon <- which(gencode_annot_exon[,"start"] <= my_init_codon_5 & gencode_annot_exon[,"end"] >= my_init_codon_5)
          my_cdna_length_A <- my_cdna_length_B <- 0

          if(length(my_snp_exon)==0){
            message("Variant not in an exonic part. Skip")
            next
          }

          # Reference sequence stats
          for(exon_i in 1:nrow(gencode_annot_exon)){
            if(exon_i==1){
              my_cdna   <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
              my_cdna   <- reverse(complement(my_cdna))
            }else{
              my_cdna_i <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
              my_cdna_i <- reverse(complement(my_cdna_i))

              my_cdna <- c(my_cdna, my_cdna_i)
            }

            # Calcul position of variant in cDNA
            if(exon_i < my_snp_exon){
              my_cdna_length_A <- my_cdna_length_A + exons_length[exon_i]
            }else if(exon_i == my_snp_exon){
              my_snp_pos_cdna <- my_cdna_length_A + (gencode_annot_exon[exon_i,"end"] - (my_snp_pos_geno-1))
            }

            # Calcul position of reference A(TG) in cDNA
            if(exon_i < my_start_exon){
              my_cdna_length_B <- my_cdna_length_B + exons_length[exon_i]
            }else if(exon_i == my_start_exon){
              my_init_codon_5_cdna <- my_cdna_length_B + (gencode_annot_exon[exon_i,"end"] - (my_init_codon_3+1))
            }
          }

          # Find codon start in reference sequence
          stats_orig <- matchPattern(morfee_data[["SEQ_INIT"]], my_cdna)

          # Found all STOP in reference sequence
          for(j in 1:length(morfee_data[["SEQ_STOP"]])){

            stats_stop_orig_j <- matchPattern(morfee_data[["SEQ_STOP"]][[j]], my_cdna)

            if(j==1){
              stats_stop_orig <- stats_stop_orig_j
            }else{
              stats_stop_orig <- c(stats_stop_orig, stats_stop_orig_j)
            }
          }

          my_ref_allele <- as.character(subseq(my_cdna, start=my_snp_pos_cdna, end=my_snp_pos_cdna))
          if(my_ref_allele!=my_nm_list[nm,5]){
            message("Mismatch between alleles")
            next
          }

          my_ref_ATG <- as.character(subseq(my_cdna, start=my_init_codon_5_cdna, end=(my_init_codon_5_cdna+2)))
          if(my_ref_ATG!="ATG"){
            message("Reference ATG no detected!")
            next
          }

          # Replace my sequence with SNP
          my_cdna_updated <- replaceLetterAt(my_cdna, my_snp_pos_cdna, my_nm_list[nm,6])

        }
        ########################################

        # Find codon start in mutated sequence
        stats_mut <- matchPattern(morfee_data[["SEQ_INIT"]], my_cdna_updated)

        # Found all STOP in mutated sequence
        for(j in 1:length(morfee_data[["SEQ_STOP"]])){

          stats_stop_mut_j <- matchPattern(morfee_data[["SEQ_STOP"]][[j]], my_cdna_updated)

          if(j==1){
            stats_stop_mut <- stats_stop_mut_j
          }else{
            stats_stop_mut <- c(stats_stop_mut, stats_stop_mut_j)
          }
        }

        # Comparer codon start in reference and mutated sequences
        new.atg <- ranges(stats_mut)[!c(ranges(stats_mut) %in% ranges(stats_orig)) ,]

        # Compare stop codons in reference and mutated sequences
        del.stop <- ranges(stats_stop_orig)[!c(ranges(stats_stop_orig) %in% ranges(stats_stop_mut)) ,]

        if(length(new.atg)>0){
          message("New ATG detected!")

          if(my_init_codon_end < my_stop_codon_end){
            my.strand <- "forward"
          }else{
            my.strand <- "reverse"
          }

          new.atg.distance <- my_init_codon_5_cdna - (start(new.atg)[1])

          test.frame <- (new.atg.distance%%3)

          if(test.frame==0){
            in.frame <- "in_frame"
          }else{
            in.frame <- paste0("out_of_frame_(",test.frame,")")
          }

          # Found all STOP
          for(j in 1:length(morfee_data[["SEQ_STOP"]])){

            my_stop_j <- matchPattern(morfee_data[["SEQ_STOP"]][[j]], my_cdna_updated[ start(range(new.atg)) : length(my_cdna_updated) ])

            if(j==1){
              my_stops <- my_stop_j
            }else{
              my_stops <- c(my_stops, my_stop_j)
            }
          }

          # First stop in phase to the new codon start
          my_stops_sort <- sort(start(ranges(my_stops)))
          my_first_stop <- my_stops_sort[my_stops_sort%%3==1][1]

          # Lenght of proteins
          generated.prot.length <- (my_first_stop-1)/3
          ref.prot.length <- (sum(gencode_annot_cds[,"end"]+1 - gencode_annot_cds[,"start"] ) -3)/3

          print( paste("For",my_gene,"-",my_nm,"and",my_snp))
          print(paste0(" - New ATG detected at: ",new.atg.distance," from the main ATG!"))
          print( paste(" - new ATG is",in.frame,"to the main ATG!"))
          print( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
          cat("\n\n")


          # Update myvcf_annot_info
          new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE"],
                                 paste0(my_nm,":",my.strand,",",new.atg.distance,",",in.frame,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)")) )
                             , collapse="|")

          myvcf_annot_info[i,"MORFEE"] <- new_field

        }# END new ATG

        if(length(del.stop)>0){
          message("STOP deletion detected!")

          # Test whether an ATG is upstream and in-frame with STOP
          # if true:  orf.stop <- TRUE
          # if false: orf.stop <- FALSE

          # Use stats_orig, but could use stats_mut
#         uatg <- ranges(stats_orig)[c(ranges(stats_orig) < ranges(del.stop)) ,]
          uatg <- start(stats_orig)[ c(start(del.stop) - start(stats_orig)) > 0]

          uatg <- c(21, uatg)# DEBUG
#         uatg <- c(20, 21)# DEBUG

          uatg_in.frame <- uatg[((start(del.stop) - uatg) %% 3)==0]

          if(length(uatg_in.frame)>0){

            # TODO: compute distance and length



            message(" -  uSTOP deletion in ORF detected!")
            print( paste("For",my_gene,"-",my_nm,"and",my_snp))
#            print(paste0(" - Deletion of a STOP codon detected at: ",del.stop.distance," from the main ATG!"))
#            print( paste(" - new ATG is",in.frame,"to the main ATG!"))
#            print( paste(" - new generated protein has a length of",stop.generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
            cat("\n\n")
          }
        }
      }
  }
  info(myvcf_annot) <- myvcf_annot_info
  return(myvcf_annot)
}
