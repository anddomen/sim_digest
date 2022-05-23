########################
#Sequence simulation digest program
#Function: condense exercise 5 of genotyping 
#by sequencing to a single function
#Author: Andrea Domen
########################

#Load needed libraries
library(BiocManager)
library(SimRAD)

#############################################
##HIGHLY ADVISED: assign test parameters to objects
#it'll make your life easier I promise

#Do you want your results to be reproducible?
#Put TRUE if yes, FALSE if no/ you don't care
seed <- TRUE

#What size (in nucleotides) do you want your simulated sequence to be?
nt_size <- 156040895
#What do you want the GC frequency to be?
gc_freq <- 0.40
#What enzyme(s) are you using?
#Supported enzymes: ApeKI, EcoRI, MseI, TaqI, SbfI, PstI, AciI, AgeI, NcoI, NotI
#if you are using a single enzyme (single digest), put NA for enzyme 2, not in quotes!
enzyme1 <- "ApeKI"
enzyme2 <- NA
#if the enzyme you want isn't supported, enter the cutsites below (in quotes)
#and put NA for enzyme 1 and enzyme 2 above
#basically just put in NA for anything you're not using
cs_5p1 <- NA
cs_3p1 <- NA
cs_5p2 <- NA
cs_3p2 <- NA
#what are the fragment sizes you're interested in? (size selection)
#Note: if you want to see the resulting graph, make sure to click 'plots'
min.frag.size <- 80
max.frag.size <- 999
#what read length do you want?
read_length_in <- 150
#is your sample paired? TRUE/FALSE
paired <- TRUE
#how many samples are you planning on having?
num.samp <- 96
#what is the max output you are expecting?
#eg for a Illumina HiSeq 3000 it is 90GB
max.output <- 90000000000

#Now that all the parameters have been specified, run the whole script
#FYI it might take a minute or two

#############################################

X_digest <- function(seed, nt_size, gc_freq, enzyme1, enzyme2, 
                     min.frag.size, max.frag.size, read_length_in, paired, 
                     num.samp, max.output) {
  #set seed to achieve reproducible results
  if(seed == TRUE){
    set.seed(818) 
  }
  #enzyme selection
  if (enzyme1 %in% c("ApeKI", "EcoRI", "MseI", "TaqI", "SbfI", "PstI", "AciI", "AgeI", "NcoI", "NotI")){
    if(enzyme1 == "ApeKI"){
      cs_5p1 <- "G"
      cs_3p1 <- "CAGC"
      cs_5p2 <- "G"
      cs_3p2 <- "CTGC"
    } else if(enzyme1 == "EcoRI")
    {
      cs_5p1 <- "G"
      cs_3p1 <- "AATTC"
    } else if(enzyme1 == "TaqI")
    {
      cs_5p1 <- "T"
      cs_3p1 <- "CGA"
    } else if(enzyme1 == "SbfI")
    {
      cs_5p1 <- "CCTGCA"
      cs_3p1 <- "GG"
    } else if(enzyme1 == "PstI")
    {
      cs_5p1 <- "CTGCA"
      cs_3p1 <- "G"
    } else if(enzyme1 == "AciI")
    {
      cs_5p1 <- "C"
      cs_3p1 <- "CGC"
    } else if(enzyme1 == "AgeI")
    {
      cs_5p1 <- "A"
      cs_3p1 <- "TGGCCA"
    } else if(enzyme1 == "NcoI")
    {
      cs_5p1 <- "C"
      cs_3p1 <- "CATGG"
    } else if(enzyme1 == "MseI")
    {
      cs_5p1 <- "T"
      cs_3p1 <- "TAA"
    } else if (enzyme1 == "NotI")
    {
      cs_5p1 <- "CG"
      cs_3p1 <- "GGCCGC"
    }}
  #second enzyme if using a double digest
  if(enzyme2 %in% c("EcoRI", "MseI", "TaqI", "SbFI", "PstI", "AciI", "AgeI", "NcoI", "NotI")){
    if(enzyme2 == "TaqI") {
      cs_5p2 <- "T"
      cs_3p2 <- "CGA"
    } else if(enzyme2 == "EcoRI")
    {
      cs_5p2 <- "G"
      cs_3p2 <- "AATTC"
    } else if (enzyme2 == "SbFI")
    {
      cs_5p2 <- "CCTGCA"
      cs_3p2 <- "GG"
    } else if(enzyme2 == "PstI")
    {
      cs_5p2 <- "CTGCA"
      cs_3p2 <- "G"
    } else if(enzyme2 == "AciI")
    {
      cs_5p2 <- "C"
      cs_3p2 <- "CGC"
    } else if(enzyme2 == "AgeI")
    {
      cs_5p2 <- "A"
      cs_3p2 <- "CCGGT"
    } else if(enzyme2 == "NcoI")
    {
      cs_5p2 <- "C"
      cs_3p2 <- "CATGG"
    } else if(enzyme2 == "MseI")
    {
      cs_5p2 <- "T"
      cs_3p2 <- "TAA"
    } else if(enzyme2 == "NotI")
    {
      cs_5p2 <- "GC"
      cs_3p2 <- "GGCCGC"
    }}
  #simulate the squence with the desired size and gc freq
  x.sim <- sim.DNAseq(size = nt_size, GCfreq = gc_freq)
  #digest the simulated chromosome with the desired cut sites
  x.sim.dig <- insilico.digest(x.sim, cs_5p1, cs_3p1, cs_5p2, cs_3p2)
  #select fragments of desired size
  x.size.sel <- size.select(x.sim.dig, min.frag.size, max.frag.size, graph = TRUE)
  ##The following 4 functions are from the cover_cals.R code file:
  computeSeq_reads <-function(frag_length, read_length, paired = FALSE) {
    if (paired == TRUE) {
      read_length <- 2 * read_length
    }
    return(min(read_length, frag_length))
  }
  # Returns the estimated coverage the digest will give, for a simulated digest
  get_genome_coverage <- function(fragments, read_length_in, paired_in, genome_size) {
    total_bases_sequencable <- get_bases_sequencable(fragments, read_length_in, paired_in)
    genome_coverage <- (total_bases_sequencable / genome_size) * 100
    return(genome_coverage)
  }
  # given a list of fragments (character vector), read length (e.g. 150), and paired_in (e.g. TRUE)
  # returns the total amount of the fragments sequencable by reads (as a single number, e.g. 1451414)
  get_bases_sequencable <- function(fragments, read_length_in, paired_in = FALSE) {
    frag_lengths <- lapply(fragments, nchar)
    
    frag_length_vector <- unlist(frag_lengths)
    
    sequencable_bases <-lapply(frag_lengths,
                               computeSeq_reads,
                               read_length = read_length_in,
                               paired = paired_in)
    
    sequencable_bases_vector <- unlist(sequencable_bases)
    
    total_bases_sequencable <- sum(sequencable_bases_vector)
    return(total_bases_sequencable)
  }
  # Calculates the depth per base on a paired end illumina sequencer
  depth <- function(fragments, read_length_in, paired_in, numb_samples, genome_size, total_gb_sequenced = 90000000) {
    
    total_bases_sequencable <- get_bases_sequencable(fragments, read_length_in, paired_in)
    genome_coverage <- (total_bases_sequencable / genome_size) * 100
    total_seq_generated <- total_gb_sequenced / numb_samples
    avg_loci_coverage <-  total_seq_generated / total_bases_sequencable
    return(avg_loci_coverage)
  }
  #back to regularly scheduled programming: find the number of bases that would
  #actually be sequenced on an Illumina Hiseq 3000
  sequenceable <- get_bases_sequencable(fragments = x.size.sel, read_length_in = read_length_in,
                                        paired_in = paired)
  #get the genome coverage you can expect
  coverage <- get_genome_coverage(fragments = x.size.sel, paired_in = paired,
                                  read_length_in = read_length_in, genome_size = nt_size)
  #find the depth for a number of samples:
  print_depth <- depth(fragments = x.size.sel, read_length_in = read_length_in, 
                       paired_in = paired, numb_samples = num.samp, genome_size = nt_size,
                       total_gb_sequenced = max.output)
  
  list(bases_sequencable = sequenceable, genome_coverage = coverage, 
       depth = print_depth) 
  }

X_digest(seed, nt_size, gc_freq, enzyme1, enzyme2, 
         min.frag.size, max.frag.size, read_length_in, paired, 
         num.samp, max.output) 
