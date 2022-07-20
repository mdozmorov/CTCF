### =========================================================================
### CTCF is an AnnotationHub package that stores genomic coordinates of predicted
### CTCF binding sites in BED format. With strand orientation (directionality of
### binding).
### -------------------------------------------------------------------------
###

# CTCF contains 51 RData objects with genome-specific CTCF GRanges.
# See https://github.com/dozmorovlab/CTCF.dev for processing information

# The object names are structured as "<assembly>.<Database>.<original database
# name or label>", e.g., "T2T.JASPAR2022_CORE_vertebrates_non_redundant_v2"

### -------------------------------------------------------------------------
### PWM matrices download and preprocessing, HPC environment
# https://github.com/dozmorovlab/CTCF.dev/blob/main/scripts/downloads.sh

# # Downloaded May 20, 2022
#
# # MEME Motif Databases (updated 20 Mar 2022)
# # https://meme-suite.org/meme/doc/download.html
# wget https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz
# # Extract
# tar zxvfh motif_databases.12.23.tgz
# # Go there
# cd motif_databases
# # File to store found CTCF motifs
# FILEOUT=CTCF_motifs.txt
# # Find CTCF motif occurrences, output in the file
# for file in `find . -type f -name "*.meme"`; do
# 	# Output only if CTCF, not CTCFL is found
# 	if grep -iv CTCFL $file | grep -iq CTCF; then
# 		echo $file >> ${FILEOUT};
# 		grep -iv CTCFL $file | grep -i CTCF >> ${FILEOUT};
# 	fi;
# done
#
# # Find CTCF motif occurrences, output on screen
# for file in `find . -type f -name "*.meme"`; do
# 	# Output only if CTCF, not CTCFL is found
# 	if grep -iv CTCFL $file | grep -iq CTCF; then
# 		echo $file;
# 		grep -iv CTCFL $file | grep -i CTCF;
# 	fi;
# done
#
# # CTCFBSDB PWM matrix
# wget https://insulatordb.uthsc.edu/download/CTCFBSDB_PWM.mat
# transfac2meme -use_acc CTCFBSDB_PWM.mat > CTCFBSDB_PWM.meme


### -------------------------------------------------------------------------
### FIMO analysis of PWM matrices, HPC environment
# https://github.com/dozmorovlab/CTCF.dev/blob/main/scripts/fimo_all.qsub

# DIRIN=/home/sequencing/juicer/Mikhail/CTCF.dev
# MOTIF=( CIS-BP_2.00_Homo_sapiens.meme CTCFBSDB_PWM_corrected.meme HOCOMOCOv11_core_HUMAN_mono_meme_format.meme JASPAR2022_CORE_vertebrates_non-redundant_v2.meme jolma2013_corrected.meme SwissRegulon_human_and_mouse.meme )
#
# SUFFIX=( T2T hg38 hg19 )
# for ASSEMBLY in ${SUFFIX[@]}; do
# 	DIRGENOME=/home/sequencing/data/ExtData/UCSC/${ASSEMBLY}/CHR
# 	GENOME=( `ls ${DIRGENOME}` )
# 	for PWM in ${MOTIF[@]}; do
# 		for CHR in ${GENOME[@]}; do
# 			DIROUT=${DIRIN}/`basename ${PWM} .meme`_${ASSEMBLY}_`basename ${CHR} .fa`
# 			fimo -o ${DIROUT} --max-stored-scores 1000000 ${DIRIN}/PWMs/${PWM} ${DIRGENOME}/${CHR}
# 		done
# 	done
# done
#
# DIRIN=/home/sequencing/juicer/Mikhail/CTCF.dev
# MOTIF=( CIS-BP_2.00_Mus_musculus.meme CTCFBSDB_PWM_corrected.meme HOCOMOCOv11_core_MOUSE_mono_meme_format.meme JASPAR2022_CORE_vertebrates_non-redundant_v2.meme jolma2013_corrected.meme SwissRegulon_human_and_mouse.meme )
#
# SUFFIX=( mm9 mm10 mm39 )
# for ASSEMBLY in ${SUFFIX[@]}; do
# 	DIRGENOME=/home/sequencing/data/ExtData/UCSC/${ASSEMBLY}/CHR
# 	GENOME=( `ls ${DIRGENOME}` )
# 	for PWM in ${MOTIF[@]}; do
# 		for CHR in ${GENOME[@]}; do
# 			DIROUT=${DIRIN}/`basename ${PWM} .meme`_${ASSEMBLY}_`basename ${CHR} .fa`
# 			fimo -o ${DIROUT} --max-stored-scores 1000000 ${DIRIN}/PWMs/${PWM} ${DIRGENOME}/${CHR}
# 		done
# 	done
# done
#
# # For Human liftover
# DIRIN=/home/sequencing/juicer/Mikhail/CTCF.dev
# MOTIF=( MA0139.1.meme )
#
# SUFFIX=( hg18 hg19 hg38 T2T )
# for ASSEMBLY in ${SUFFIX[@]}; do
# 	DIRGENOME=/home/sequencing/data/ExtData/UCSC/${ASSEMBLY}/CHR
# 	GENOME=( `ls ${DIRGENOME}` )
# 	for PWM in ${MOTIF[@]}; do
# 		for CHR in ${GENOME[@]}; do
# 			DIROUT=${DIRIN}/`basename ${PWM} .meme`_${ASSEMBLY}_`basename ${CHR} .fa`
# 			fimo -o ${DIROUT} --max-stored-scores 1000000 ${DIRIN}/PWMs/${PWM} ${DIRGENOME}/${CHR}
# 		done
# 	done
# done
#
# # For Mouse liftover
# DIRIN=/home/sequencing/juicer/Mikhail/CTCF.dev
# MOTIF=( MA0139.1.meme )
#
# SUFFIX=( mm9 mm10 mm39 )
# for ASSEMBLY in ${SUFFIX[@]}; do
# 	DIRGENOME=/home/sequencing/data/ExtData/UCSC/${ASSEMBLY}/CHR
# 	GENOME=( `ls ${DIRGENOME}` )
# 	for PWM in ${MOTIF[@]}; do
# 		for CHR in ${GENOME[@]}; do
# 			DIROUT=${DIRIN}/`basename ${PWM} .meme`_${ASSEMBLY}_`basename ${CHR} .fa`
# 			fimo -o ${DIROUT} --max-stored-scores 1000000 ${DIRIN}/PWMs/${PWM} ${DIRGENOME}/${CHR}
# 		done
# 	done
# done

### -------------------------------------------------------------------------
### Download and Process https://screen.encodeproject.org/ data
# https://github.com/dozmorovlab/CTCF.dev/blob/main/02_EDA_SCREEN.Rmd

# Project folder path
dir_project <- "/Users/mdozmorov/Documents/Work/GitHub/CTCF.dev"
# Input files
# Human CTCF-bound cCRE (hg38)
fileNameIn1 <- file.path(dir_project, "data/GRCh38-CTCF.bed")
if (!file.exists(fileNameIn1)) {
  download.file("https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-CTCF.bed", fileNameIn1)
}
# Mouse CTCF-bound cCRE (mm10)
fileNameIn2 <- file.path(dir_project, "data/mm10-CTCF.bed")
if (!file.exists(fileNameIn2)) {
  download.file("https://api.wenglab.org/screen_v13/fdownloads/cCREs/mm10-CTCF.bed", fileNameIn2)
}
# Output files
fileNameOut1 <- file.path(dir_project, "RData", "hg38.SCREEN.GRCh38_CTCF.RData")
fileNameOut2 <- file.path(dir_project, "RData", "mm10.SCREEN.mm10_CTCF.RData")
# Read in data
mtx_human <- read_tsv(fileNameIn1, col_names = FALSE)
mtx_mouse <- read_tsv(fileNameIn2, col_names = FALSE)
# Function to convert SCREEN data to GRanges
mtx_to_gr <- function(mtx = mtx_human, genome_id = "hg38") {
  gr <- GRanges(seqnames = mtx$X1,
                IRanges(start = mtx$X2,
                        end = mtx$X3))

  gr$ID1 <- mtx$X4
  gr$ID2 <- mtx$X5
  gr$Type <- mtx$X6

  # Get chromosome info and match it to the chromosome order in ctcfBED
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
  # Subset to common chromosomes
  common_chromosomes <- intersect(chrom_data$chrom, seqlevels(gr))
  chrom_data <- chrom_data[chrom_data$chrom %in% common_chromosomes, ]
  gr <- keepSeqlevels(gr, common_chromosomes, pruning.mode = "tidy")
  # Match order
  chrom_data <- chrom_data[match(seqlevels(gr), chrom_data$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(gr), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(gr) <- chrom_data$size
  isCircular(gr) <- chrom_data$circular
  genome(gr)     <- genome_id

  return(gr)
}
# Convert SCREEN data to GRanges
hg38.SCREEN.GRCh38_CTCF <- mtx_to_gr(mtx = mtx_human, genome_id = "hg38")
mm10.SCREEN.mm10_CTCF <- mtx_to_gr(mtx = mtx_mouse, genome_id = "mm10")
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = "hg38.SCREEN.GRCh38_CTCF", file = fileNameOut1)
save(list = "mm10.SCREEN.mm10_CTCF", file = fileNameOut2)


### -------------------------------------------------------------------------
### Download and Process CTCFBSDB, predicted data
# https://github.com/dozmorovlab/CTCF.dev/blob/main/03_EDA_CTCFBSDB.Rmd
# Project folder path
dir_project <- "/Users/mdozmorov/Documents/Work/GitHub/CTCF.dev"

# Input files
# Predicted data
fileNameIn1 <- file.path(dir_project, "data/allcomp.txt.gz")
if (!file.exists(fileNameIn1)) {
  download.file("https://insulatordb.uthsc.edu/download/allcomp.txt.gz", fileNameIn1)
}
# Experimental data
fileNameIn2 <- file.path(dir_project, "data/CTCFBSDB_all_exp_sites_Sept12_2012.txt.gz")
# download.file("https://insulatordb.uthsc.edu/download/CTCFBSDB_all_exp_sites_Sept12_2012.txt.gz", fileNameIn2)
# Output files
fileNameOut1 <- file.path(dir_project, "RData", "hg18.CTCFBSDB.CTCF_predicted_human.RData")
fileNameOut2 <- file.path(dir_project, "RData", "mm8.CTCFBSDB.CTCF_predicted_mouse.RData")
# liftOver human
fileNameOut3 <- file.path(dir_project, "RData", "hg19.CTCFBSDB.CTCF_predicted_human.RData")
fileNameOut4 <- file.path(dir_project, "RData", "hg38.CTCFBSDB.CTCF_predicted_human.RData")
# liftOver mouse
fileNameOut5 <- file.path(dir_project, "RData", "mm9.CTCFBSDB.CTCF_predicted_mouse.RData")
fileNameOut6 <- file.path(dir_project, "RData", "mm10.CTCFBSDB.CTCF_predicted_mouse.RData")
# Read in data
mtx <- read_tsv(fileNameIn1)
# Function to convert CTCFBSDB data to GRanges
mtx_to_gr <- function(mtx = mtx_human, genome_id = "hg19") {
  gr <- GRanges(seqnames = sapply(mtx$`Chromosome Location`, function(x) strsplit(x, ":|-")[[1]][1]),
                IRanges(start = sapply(mtx$`Chromosome Location`, function(x) strsplit(x, ":|-")[[1]][2]) %>% as.numeric(),
                        end = sapply(mtx$`Chromosome Location`, function(x) strsplit(x, ":|-")[[1]][3]) %>% as.numeric()))

  gr$`5PrimeGene` <- mtx$`5' Flanking Gene`
  gr$`3PrimeGene` <- mtx$`3' Flanking Gene`

  # Get chromosome info and match it to the chromosome order in ctcfBED
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
  # Subset to common chromosomes
  common_chromosomes <- intersect(chrom_data$chrom, seqlevels(gr))
  chrom_data <- chrom_data[chrom_data$chrom %in% common_chromosomes, ]
  gr <- keepSeqlevels(gr, common_chromosomes, pruning.mode = "tidy")
  # Match order
  chrom_data <- chrom_data[match(seqlevels(gr), chrom_data$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(gr), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(gr) <- chrom_data$size
  isCircular(gr) <- chrom_data$circular
  genome(gr)     <- genome_id

  return(gr)
}
# Convert CTCFBSDB data to GRanges
hg18.CTCFBSDB.CTCF_predicted_human <- mtx_to_gr(mtx = mtx_human, genome_id = "hg18")
mm8.CTCFBSDB.CTCF_predicted_mouse <- mtx_to_gr(mtx = mtx_mouse, genome_id = "mm8")
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = "hg18.CTCFBSDB.CTCF_predicted_human", file = fileNameOut1)
save(list = "mm8.CTCFBSDB.CTCF_predicted_mouse", file = fileNameOut2)
# Function to LiftOver hg18-hg19-hg38, mm8-mm9-mm10
# Downloads a chain for the gf18_to_XX conversion and create a GRanges object
liftOver_custom <- function(URL = "https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz", gr_in = hg18.CTCFBSDB.CTCF_predicted_human, genome_id = "hg19"){
  fileChain <- file.path(dir_project, "data", basename(URL))
  if (!file.exists(str_remove(fileChain, ".gz"))) {
    download.file(URL, destfile = fileChain)
    R.utils::gunzip(fileChain)
  }

  ch <- import.chain(str_remove(fileChain, ".gz"))
  gr_converted <- liftOver(gr_in, ch) %>% unlist()
  # Get chromosome info and match it to the chromosome order in ctcfBED
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
  # Subset to common chromosomes
  common_chromosomes <- intersect(chrom_data$chrom, seqlevels(gr_converted))
  chrom_data <- chrom_data[chrom_data$chrom %in% common_chromosomes, ]
  gr_converted <- keepSeqlevels(gr_converted, common_chromosomes, pruning.mode = "tidy")
  # Match order
  chrom_data <- chrom_data[match(seqlevels(gr_converted), chrom_data$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(gr_converted), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(gr_converted) <- chrom_data$size
  isCircular(gr_converted) <- chrom_data$circular
  genome(gr_converted)     <- genome_id
  return(gr_converted)
}
# hg18.CTCFBSDB.CTCF_predicted_human to hg19 conversion
gr_out <- "hg19.CTCFBSDB.CTCF_predicted_human"
# Assign the converted GRanges to the variable name
assign(gr_out, liftOver_custom(URL = "https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz", gr_in = hg18.CTCFBSDB.CTCF_predicted_human, genome_id = "hg19"))
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = gr_out, file = fileNameOut3)
# hg18.CTCFBSDB.CTCF_predicted_human to hg38 conversion
gr_out <- "hg38.CTCFBSDB.CTCF_predicted_human"
# Assign the converted GRanges to the variable name
assign(gr_out, liftOver_custom(URL = "https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz", gr_in = hg18.CTCFBSDB.CTCF_predicted_human, genome_id = "hg38"))
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = gr_out, file = fileNameOut4)
# mm8.CTCFBSDB.CTCF_predicted_mouse to mm9 conversion
gr_out <- "mm9.CTCFBSDB.CTCF_predicted_mouse"
# Assign the converted GRanges to the variable name
assign(gr_out, liftOver_custom(URL = "https://hgdownload.cse.ucsc.edu/goldenpath/mm8/liftOver/mm8ToMm9.over.chain.gz", gr_in = mm8.CTCFBSDB.CTCF_predicted_mouse, genome_id = "mm9"))
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = gr_out, file = fileNameOut5)
# mm8.CTCFBSDB.CTCF_predicted_mouse to mm10 conversion
gr_out <- "mm10.CTCFBSDB.CTCF_predicted_mouse"
# Assign the converted GRanges to the variable name
assign(gr_out, liftOver_custom(URL = "https://hgdownload.cse.ucsc.edu/goldenpath/mm8/liftOver/mm8ToMm10.over.chain.gz", gr_in = mm8.CTCFBSDB.CTCF_predicted_mouse, genome_id = "mm10"))
# Save as RData object. subdir is the character name of the CTCF GRanges variable
save(list = gr_out, file = fileNameOut6)


### -------------------------------------------------------------------------
### Processing FIMO chromosome-specific results processed on an HPC cluster
# https://github.com/dozmorovlab/CTCF.dev/blob/main/04_FIMO_processing.Rmd
# Project folder path
dir_data <- "/Users/mdozmorov/Documents/Data/GoogleDrive/CTCF.dev/merlot"
dir_project <- "/Users/mdozmorov/Documents/Work/GitHub/CTCF.dev/"
# Motifs
motif_human  <- c("CIS_BP_2.00_Homo_sapiens", "CTCFBSDB_PWM", "HOCOMOCOv11_core_HUMAN_mono_meme_format", "JASPAR2022_CORE_vertebrates_non_redundant_v2", "Jolma2013", "SwissRegulon_human_and_mouse")
motif_mouse  <- c("CIS_BP_2.00_Mus_musculus", "CTCFBSDB_PWM", "HOCOMOCOv11_core_MOUSE_mono_meme_format", "JASPAR2022_CORE_vertebrates_non_redundant_v2", "Jolma2013", "SwissRegulon_human_and_mouse")
# Genomes
genome_human <- c("T2T", "hg38", "hg19")
genome_mouse <- c("mm10", "mm9", "mm39")
# Chromosomes
chromosome_human <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chromosome_mouse <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
# Function to combine database- and chromosome-specific FIMO results for a given genome assembly
combine_fimo <- function(motif = "SwissRegulon_human_and_mouse", genome_id = "T2T", chromosomes = chromosome_human) {
  ctcfGR <- c()
  for(chrom in chromosomes) {
    # print(chrom)
    fileNameIn <- file.path(dir_data, "results", paste(motif, genome_id, chrom, sep = "_"), "fimo.txt.gz")
    # Read "fimo.bed" created by "fimo.qsub"
    ctcfBED <- read_tsv(fileNameIn, col_types = c("ccddcdddc"))
    # Convert to GRanges object
    ctcf_chr_gr <- GRanges(seqnames = ctcfBED$`sequence name`,
                           IRanges(start = ctcfBED$start, end = ctcfBED$stop),
                           strand = ctcfBED$strand)
    # Add metadata
    ctcf_chr_gr$name = ctcfBED$`#pattern name`
    ctcf_chr_gr$score = ctcfBED$score
    ctcf_chr_gr$pvalue = ctcfBED$`p-value`
    ctcf_chr_gr$qvalue = ctcfBED$`q-value`
    ctcf_chr_gr$sequence = ctcfBED$`matched sequence`
    ctcfGR <- c(ctcfGR, ctcf_chr_gr)
  }
  # Combine all chromosomes
  ctcfGR <- do.call(c, as(ctcfGR, "GRangesList"))
  # Sort
  ctcfGR <- ctcfGR %>% sort()

  # Add seqinfo
  # Parse out genome ID from the folder name, to get hg19, hg38, mm9, or mm10
  if (genome_id == "T2T") {
    # Seqinfor for T2T genome
    chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = "GCA_009914755.4")
    chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
    chrom_data <- chrom_data[chrom_data$AssignedMolecule %in% seqlevels(ctcfGR), ]
    chrom_data <- chrom_data[match(seqlevels(ctcfGR), chrom_data$AssignedMolecule), ]
    # Check if chromosome order is the same
    if (!all.equal(seqlevels(ctcfGR), chrom_data$AssignedMolecule)) {
      print(paste("Chromosome order does not match for", genome_id, "genome."))
      break
    }
    # Assign seqinfo data
    seqlengths(ctcfGR) <- chrom_data$SequenceLength
    isCircular(ctcfGR) <- ifelse(is.na(chrom_data$circular), FALSE, TRUE)
    genome(ctcfGR)     <- "GCA_009914755.4"
  } else {
    # Get chromosome info and match it to the chromosome order in ctcfBED
    chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
    chrom_data <- chrom_data[chrom_data$chrom %in% seqlevels(ctcfGR), ]
    chrom_data <- chrom_data[match(seqlevels(ctcfGR), chrom_data$chrom), ]
    # Check if chromosome order is the same
    if (!all.equal(seqlevels(ctcfGR), chrom_data$chrom)) {
      print(paste("Chromosome order does not match for", genome_id, "genome."))
      break
    }
    # Assign seqinfo data
    seqlengths(ctcfGR) <- chrom_data$size
    isCircular(ctcfGR) <- chrom_data$circular
    genome(ctcfGR)     <- genome_id
  }

  ctcfGR_to_save <- paste(genome_id, motif, sep = ".")
  fileNameOut1 <- file.path(dir_project, "RData", paste0(ctcfGR_to_save, ".RData"))
  fileNameOut1.1 <- file.path(dir_project, "RData", paste0(ctcfGR_to_save, ".bed"))
  # Assign this object to the subfolder-specific variable name
  assign(ctcfGR_to_save, ctcfGR)
  # Save as RData object. ctcfGR_to_save is the character name of the CTCF GRanges variable
  save(list = ctcfGR_to_save, file = fileNameOut1)
  # load(file = fileNameOut1)
  # export.bed(eval(parse(text = ctcfGR_to_save)), fileNameOut1.1)
  ctcfDF_to_save <- as.data.frame(eval(parse(text = ctcfGR_to_save)))
  ctcfDF_to_save <- ctcfDF_to_save[, c("seqnames", "start", "end", "name", "score", "strand", "width", "pvalue", "qvalue", "sequence")]
  write_tsv(ctcfDF_to_save, fileNameOut1.1, col_names = FALSE)
}
# Human data for liftOver
for (gen in c("hg18", genome_human)) {
  for (mot in "MA0139.1") {
    print(paste(gen, mot, sep = "."))
    combine_fimo(motif = mot, genome_id = gen, chromosomes = chromosome_human)
  }
}
# Human data all
for (gen in genome_human) {
  for (mot in motif_human) {
    print(paste(gen, mot, sep = "."))
    combine_fimo(motif = mot, genome_id = gen, chromosomes = chromosome_human)
  }
}
# Mouse data for liftOver
for (gen in c(genome_mouse)) {
  for (mot in "MA0139.1") {
    print(paste(gen, mot, sep = "."))
    combine_fimo(motif = mot, genome_id = gen, chromosomes = chromosome_mouse)
  }
}
# Mouse data all
for (gen in genome_mouse) {
  for (mot in motif_mouse) {
    print(paste(gen, mot, sep = "."))
    combine_fimo(motif = mot, genome_id = gen, chromosomes = chromosome_mouse)
  }
}


