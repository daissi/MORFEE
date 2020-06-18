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
  myvcf_annot_info$MORFEE_uATG <- NA
  myvcf_annot_info$MORFEE_uSTOP <- NA
  myvcf_annot_header <- rbind(info(header(myvcf_annot)),
                              data.frame(Number = ".", Type = "String", Description = "New ATG annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "Deletion STOP annotation provided by MORFEE", stringsAsFactors = FALSE) )

  rownames(myvcf_annot_header) <- c(rownames(info(header(myvcf_annot))),"MORFEE_uATG", "MORFEE_uSTOP")
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

      my_enst <- morfee_data[["GENCODE_METAD"]][my_nm_id,1][1]
      my_transcript_id <- grep(paste0(my_enst,"_"), morfee_data[["GENCODE_ANNOT"]]$transcript_id)
      gencode_annot_sub <- morfee_data[["GENCODE_ANNOT"]][my_transcript_id,]
      gencode_annot_cds <- gencode_annot_sub[gencode_annot_sub$type=="CDS",]
      gencode_annot_exon <- gencode_annot_sub[gencode_annot_sub$type=="exon",]
      gencode_annot_transcript_type <- unique(gencode_annot_exon$transcript_type)

      if(length(gencode_annot_transcript_type)<1){
        message("Transcript type is not protein coding!")
        next
      }else if(!(gencode_annot_transcript_type %in% "protein_coding")){
        message("Transcript type is not protein coding!")
        next
      }

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

          message( paste("For",my_gene,"-",my_nm,"and",my_snp))
          message(paste0(" - New uATG detected at: ",new.atg.distance," from the main ATG!"))
          message( paste(" - new uATG is",in.frame,"to the main ATG!"))
          message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
          message("\n\n")


          # Update myvcf_annot_info
          new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE_uATG"],
                                 paste0(my_nm,":",my.strand,",",new.atg.distance,",",in.frame,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)")) )
                             , collapse="|")

          myvcf_annot_info[i,"MORFEE_uATG"] <- new_field

        }# END new ATG

        if(length(del.stop)>0){

          # Use stats_orig, but could use stats_mut
          uatg <- start(stats_orig)[ c(start(del.stop) - start(stats_orig)) > 0]

          uatg_in_frame <- uatg[((start(del.stop) - uatg) %% 3)==0]

          del.stop.distance <- my_init_codon_5_cdna - (start(del.stop)[1])

          if(length(uatg_in_frame)>0){

            if(my_init_codon_end < my_stop_codon_end){
              my.strand <- "forward"
            }else{
              my.strand <- "reverse"
            }

              print(       "STOP deletion detected!")
              print(       " -  uSTOP deletion in ORF detected!")
              print( paste("For",my_gene,"-",my_nm,"and",my_snp))
              print(paste0(" - Deletion of a uSTOP codon detected at: ",-del.stop.distance," from the main ATG!"))
              print( paste(" --- "   ,as.character(        my_cdna[start(del.stop)[1]:end(del.stop)[1]] ),
                           " becomes ",as.character(my_cdna_updated[start(del.stop)[1]:end(del.stop)[1]] ) ))
              print( paste(" --- Gene direction:",my.strand))


            # several uATG could be present, so the protein length will be different
            for(uatg_i in uatg_in_frame){
            # uatg_i = uatg_in_frame[1]

              # Find next stop in frame with uatg_i
              uatg_i_in_frame <- start(stats_stop_mut)[ ((uatg_i - start(stats_stop_mut)) %%3)==0 ]
              first_new_stop <- min( uatg_i_in_frame[uatg_i_in_frame > uatg_i] )

              # Compute distance and length
              stop.generated.prot.length <- (first_new_stop-uatg_i)/3
              ref.prot.length <- (sum(gencode_annot_cds[,"end"]+1 - gencode_annot_cds[,"start"] ) -3)/3

              uatg_used <- -(my_init_codon_5_cdna - uatg_i)
              stop_used <- (my_init_codon_5_cdna - 1 - first_new_stop)

              if(uatg_used>=0){
                message("Position of uATG is positive! Probably an error in the used reference database")
                next
              }

              if(stop_used<0){

                overlapping.perc <- (-stop_used/(ref.prot.length*3))*100

                overlapping.prot <- paste0("overlapping_",round(overlapping.perc, digits = 2),"%")
              }else{
                overlapping.prot <- "not_overlapping"
              }

              # FIXME: add case: overlapping.prot <- "elongated"

              stop.codon <- as.character(stats_stop_mut[start(stats_stop_mut)==first_new_stop])

              print(       " --")
              print( paste(" --- using uATG at",uatg_used,"to the main ATG!"))
              print(paste0(" --- using STOP (",stop.codon,") at ",-stop_used," to the main ATG!"))
              print( paste(" --- new predicted ORF has a length of",stop.generated.prot.length,"(aa) vs",ref.prot.length,"(aa) for the main protein"))
              print( paste(" --- new predicted ORF is",overlapping.prot,"with the main protein"))
              print(paste0(" - DEBUG: i=",i," ; nm=",nm))

              # Update myvcf_annot_info
              new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE_uSTOP"],
                                     paste0(my_nm,":",my.strand,",",-del.stop.distance,",",overlapping.prot,",",stop.generated.prot.length,"[/",ref.prot.length,"]","(aa)")) )
                                 , collapse="|")

              myvcf_annot_info[i,"MORFEE_uSTOP"] <- new_field
            }
            cat("\n\n")

          }else{

              message(       "STOP deletion detected!")
              message(       " -  uSTOP deletion detected BUT without an upstream ATG (not in an ORF region)!")
              message( paste("For",my_gene,"-",my_nm,"and",my_snp))
              message(paste0(" - Deletion of a uSTOP codon detected at: ",-del.stop.distance," from the main ATG!"))
              message( paste(" --- "   ,as.character(        my_cdna[start(del.stop)[1]:end(del.stop)[1]] ),
                             " becomes ",as.character(my_cdna_updated[start(del.stop)[1]:end(del.stop)[1]] ) ))
              message("\n\n")

          }
        }# END del STOP
      }
  }
  info(myvcf_annot) <- myvcf_annot_info
  return(myvcf_annot)
}
