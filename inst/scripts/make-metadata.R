### =========================================================================
### CTCF metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
    Title = c("CTCF_hg19.RData",
              "CTCF_hg38.RData",
              "CTCF_mm9.RData",
              "CTCF_mm10.RData"),
    Description = c("hg19 genomic coordinates of CTCF binding motif MA0139.1, detected by FIMO",
                    "hg38 genomic coordinates of CTCF binding motif MA0139.1, detected by FIMO",
                    "mm9 genomic coordinates of CTCF binding motif MA0139.1, detected by FIMO",
                    "mm10 genomic coordinates of CTCF binding motif MA0139.1, detected by FIMO"),
    BiocVersion = rep("3.12", 4),
    Genome = c("hg19", "hg38", "mm9", "mm10"),
    SourceType = rep("RData", 4),
    SourceUrl = rep("https://drive.google.com/drive/folders/19ZXr7IETfks0OdYlmuc1Hqe700Pw3jPc?usp=sharing", 4),
    SourceVersion = rep("May 8 2021", 4),
    Species = c(rep("Homo sapiens", 2), rep("Mus musculus", 2)),
    TaxonomyId = c(rep(9606, 2), rep(10090, 2)),
    Coordinate_1_based = rep(TRUE, 4),
    DataProvider = rep(paste("UCSC", "Jaspar"), 4),
    Maintainer = rep("Mikhail Dozmorov <mikhail.dozmorov@gmail.com>", 4),
    RDataClass = rep("GRanges", 4),
    DispatchClass = rep("Rda", 4),
    RDataPath = c("CTCF/CTCF_hg19.RData",
                  "CTCF/CTCF_hg38.RData",
                  "CTCF/CTCF_mm9.RData",
                  "CTCF/CTCF_mm10.RData"),
    Tags = c(paste("hg19"),
             paste("hg38"),
             paste("mm9"),
             paste("mm10")),
    Notes = "")

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
