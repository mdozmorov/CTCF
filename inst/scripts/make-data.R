### =========================================================================
### CTCF is an AnnotationHub package that stores genomic coordinates of predicted
### CTCF binding sites in BED format. With strand orientation (directionality of
### binding).
### -------------------------------------------------------------------------
###

# CTCF contains 4 RData objects with genome-specific CTCF GRanges.

# The object names are structured as "CTCF_"<genome version abbreviation>, e.g.,
# "CTCF_hg19".  The input data include:

# # CTCF motif in MEME format. Downloaded 03/31/2021
# http://jaspar.genereg.net/matrix/MA0139.1/
# wget http://jaspar.genereg.net/api/v1/matrix/MA0139.1.meme

# # hg19 human genome. Downloaded 12/03/2015
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa > hg19.fa

# # hg38 human genome. Downloaded 12/03/2015
# wget hhttp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
# cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa > hg38.fa

# # mm9 genome. Downloaded 11/23/2016
# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
# cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrX.fa chrY.fa chrM.fa > mm9.fa

# # mm10 genome. Downloaded 02/19/2016
# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
# cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrX.fa chrY.fa chrM.fa > mm10.fa

# # The CTCF motif coordinates were determined by FIMO
# conda create -n meme -y
# conda install -c bioconda meme bedtools
# source activate meme
#
# # The GENOME and DIROUT variables were adjusted for each genome
# DIRIN=/home/sequencing/juicer/Mikhail/CTCF_fimo
# MOTIF=${DIRIN}/MA0139.1.meme
# GENOME=/home/sequencing/data/ExtData/UCSC/hg19/hg19.fa
# DIROUT=${DIRIN}/CTCF_hg19
# fimo -o ${DIROUT} ${MOTIF} ${GENOME}
# # Sorting using bedtools
# FILEIN=${DIROUT}/fimo.txt
# FILEOUT=${DIROUT}/fimo.bed
# # chr,start,stop,motif,score,strand,p-value,q-value,sequence
# cat ${FILEIN} | sed '1d' | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$1,$6,$5,$7,$8,$9}' | bedtools sort -i - > ${FILEOUT}
# conda deactivate

# The following example demonstrate how the CTCF coordinates were converted
# into RData objects

# Download results data from https://drive.google.com/drive/folders/19ZXr7IETfks0OdYlmuc1Hqe700Pw3jPc?usp=sharing
# Folder with results
dir_in <- "/Volumes/GoogleDrive/My Drive/CTCF_fimo"
# Subfolders named like "CTCF_hg19
subdirs <- list.dirs(path = dir_in, full.names = FALSE, recursive = FALSE)
# In each subfolder
for (subdir in subdirs) {
  # Read "fimo.bed" created by "fimo.qsub"
  ctcfBED <- read.delim(file.path(dir_in, subdir, "fimo.bed"), sep = "\t",
                        col.names = c("chr", "start", "stop", "motif", "score",
                                      "strand", "p-value","q-value","sequence"),
                        header = FALSE, stringsAsFactors = FALSE)
  # Convert to GRanges object
  ctcfGR <- GenomicRanges::makeGRangesFromDataFrame(ctcfBED, keep.extra.columns = TRUE)

  # Add seqinfo
  # Parse out genome ID from the folder name, to get hg19, hg38, mm9, or mm10
  genome_id <- sub("CTCF_", "", subdir)
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

  # Assign this object to the subfolder-specific variable name
  assign(subdir, ctcfGR)
  # Save as RData object. subdir is the character name of the CTCF GRanges variable
  save(list = subdir, file = paste0(subdir, ".RData"))
  # load(file = paste0(subdir, ".RData"))
}
