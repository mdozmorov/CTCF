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
    BiocVersion = c(rep("3.12", 4)),
    Genome = c("hg19", "hg38", "mm9", "mm10"),
    SourceType = c(rep("RData", 4)),
    SourceUrl = "https://drive.google.com/drive/folders/19ZXr7IETfks0OdYlmuc1Hqe700Pw3jPc?usp=sharing",
    SourceVersion = "May 8 2021",
    Species = c("Homo sapiens", "Mus musculus"),
    TaxonomyId = c(9606, 10090),
    Coordinate_1_based = TRUE,
    DataProvider = c("UCSC", "Jaspar"),
    Maintainer = "Mikhail Dozmorov <mikhail.dozmorov@gmail.com>",
    RDataClass = c(rep("GRanges", 4)),
    DispatchClass = c(rep("Rda", 4)),
    RDataPath = c("CTCF/CTCF_hg19.RData",
                  "CTCF/CTCF_hg38.RData",
                  "CTCF/CTCF_mm9.RData",
                  "CTCF/CTCF_mm10.RData"),
    Tags = "",
    Notes = "")

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
