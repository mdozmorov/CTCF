
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CTCF

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`CTCF` defines an AnnotationHub resource representing genomic
coordinates of FIMO-predicted CTCF binding sites with motif MA0139.1
(Jaspar).

-   Human (hg19, hg38) and mouse (mm9, mm10) genomes.
-   The binding sites were detected using the FIMO tool of the MEME
    suite using default settings.
-   Extra columns include motif name (MA0139.1), score, p-value,
    q-value, and the motif sequence.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `CTCF` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("CTCF")
```

## Example

``` r
suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
#> Warning: DEPRECATION: As of AnnotationHub (>2.23.2), default caching location has changed.
#>   Problematic cache: /Users/mdozmorov/Library/Caches/AnnotationHub
#>   See https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html#default-caching-location-update
#> snapshotDate(): 2021-05-18
query_data <- query(ah, "CTCF")
query_data
#> AnnotationHub with 466 records
#> # snapshotDate(): 2021-05-18
#> # $dataprovider: UCSC, Haemcode, UCSC Jaspar, Pazar
#> # $species: Homo sapiens, Mus musculus, NA
#> # $rdataclass: GRanges, BigWigFile
#> # additional mcols(): taxonomyid, genome, description,
#> #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#> #   rdatapath, sourceurl, sourcetype 
#> # retrieve records with, e.g., 'object[["AH22248"]]' 
#> 
#>             title                                             
#>   AH22248 | pazar_CTCF_Cui_20120522.csv                       
#>   AH22249 | pazar_CTCF_HEPG2_Schmidt_20120522.csv             
#>   AH22519 | wgEncodeAwgTfbsBroadDnd41CtcfUniPk.narrowPeak.gz  
#>   AH22521 | wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak.gz
#>   AH22524 | wgEncodeAwgTfbsBroadH1hescCtcfUniPk.narrowPeak.gz 
#>   ...       ...                                               
#>   AH28453 | CTCF_GSM918744_Immortalized_Erythroid.csv         
#>   AH95565 | CTCF_hg19.RData                                   
#>   AH95566 | CTCF_hg38.RData                                   
#>   AH95567 | CTCF_mm9.RData                                    
#>   AH95568 | CTCF_mm10.RData
```

The FIMO-predicted CTCF sites are named as
“CTCF\_<genome version abbreviation>”, e.g., “CTCF\_hg38”. Use
`query_data <- query(ah , "CTCF_hg38")` for a more targeted search.

We can check the details about the object.

``` r
query_data["AH95566"]
#> AnnotationHub with 1 record
#> # snapshotDate(): 2021-05-18
#> # names(): AH95566
#> # $dataprovider: UCSC Jaspar
#> # $species: Homo sapiens
#> # $rdataclass: GRanges
#> # $rdatadateadded: 2021-05-18
#> # $title: CTCF_hg38.RData
#> # $description: hg38 genomic coordinates of CTCF binding motif MA0139.1, det...
#> # $taxonomyid: 9606
#> # $genome: hg38
#> # $sourcetype: RData
#> # $sourceurl: https://drive.google.com/drive/folders/19ZXr7IETfks0OdYlmuc1Hq...
#> # $sourcesize: NA
#> # $tags: c("FunctionalAnnotation", "GenomicSequence", "hg38") 
#> # retrieve record with 'object[["AH95566"]]'
```

And retrieve the object.

``` r
CTCF_hg38 <- query_data[["AH95566"]]
#> loading from cache
CTCF_hg38
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> GRanges object with 56049 ranges and 5 metadata columns:
#>           seqnames            ranges strand |       motif     score   p.value
#>              <Rle>         <IRanges>  <Rle> | <character> <numeric> <numeric>
#>       [1]     chr1       11223-11241      - |    MA0139.1   24.4754  1.34e-09
#>       [2]     chr1       11281-11299      - |    MA0139.1   22.7377  1.01e-08
#>       [3]     chr1       24782-24800      - |    MA0139.1   17.3770  7.11e-07
#>       [4]     chr1       91420-91438      + |    MA0139.1   16.2951  1.41e-06
#>       [5]     chr1     104985-105003      - |    MA0139.1   16.7869  1.04e-06
#>       ...      ...               ...    ... .         ...       ...       ...
#>   [56045]     chrY 57044316-57044334      - |    MA0139.1   16.4590  1.27e-06
#>   [56046]     chrY 57189659-57189677      + |    MA0139.1   15.7541  1.95e-06
#>   [56047]     chrY 57203409-57203427      - |    MA0139.1   15.6393  2.09e-06
#>   [56048]     chrY 57215279-57215297      + |    MA0139.1   19.5738  1.53e-07
#>   [56049]     chrY 57215337-57215355      + |    MA0139.1   24.4754  1.34e-09
#>             q.value            sequence
#>           <numeric>         <character>
#>       [1]    0.0216 TCGCCAGCAGGGGGCGCCC
#>       [2]    0.0398 GCGCCAGCAGGGGGCGCTG
#>       [3]    0.2350 CGTCCAGCAGATGGCGGAT
#>       [4]    0.3080 GTGGCACCAGGTGGCAGCA
#>       [5]    0.2750 CCAACAGCAGGTGGCAGCC
#>       ...       ...                 ...
#>   [56045]    0.2990 TGGTCACCTGGGGGCACTA
#>   [56046]    0.3430 TGTCCTCTAGGGGTCAGCC
#>   [56047]    0.3510 CTGCCGCAAGGGGGCGCAT
#>   [56048]    0.1190 gcgccacgagggggcggtg
#>   [56049]    0.0216 tcgccagcagggggcgccc
#>   -------
#>   seqinfo: 24 sequences from hg38 genome
```

Note that the default q-value cutoff is 0.5. Looking at the q-value
distribution:

<img src="man/figures/CTCF_hg38_qvalue.png" width="100%" />

one may decide to use a more stringent cutoff. E.g., filtering by
q-value less than 0.3 filters out more than half of the predicted sites.
The remaining sites may be considered as high-confidence CTCF sites.

``` r
# Check length before filtering
length(CTCF_hg38)
#> [1] 56049
# Filter and check length after filtering
CTCF_hg38 <- CTCF_hg38[CTCF_hg38$q.value < 0.3]
length(CTCF_hg38)
#> [1] 25474
```

## CTCF GRanges for other organisms

``` r
# hg19 CTCF coordinates
CTCF_hg19 <- query_data[["AH95565"]]
# mm9 CTCF coordinates
CTCF_mm9 <- query_data[["AH95567"]]
# mm10 CTCF coordinates
CTCF_mm10 <- query_data[["AH95568"]]
```

See [inst/scripts/make-data.R](inst/scripts/make-data.R) how to create
the CTCF GRanges objects.

## Citation

Below is the citation output from using `citation('CTCF')` in R. Please
run this yourself to check for any updates on how to cite **CTCF**.

``` r
print(citation("CTCF"), bibtex = TRUE)
#> 
#> Dozmorov MG, Davis E, Mu W, Lee S, Triche T, Phanstiel D, Love M
#> (2021). _CTCF_. https://github.com/mdozmorov/CTCF/CTCF - R package
#> version 0.99.4, <URL: https://github.com/mdozmorov/CTCF>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {CTCF},
#>     author = {Mikhail G. Dozmorov and Eric Davis and Wancen Mu and Stuart Lee and Tim Triche and Douglas Phanstiel and Michael Love},
#>     year = {2021},
#>     url = {https://github.com/mdozmorov/CTCF},
#>     note = {https://github.com/mdozmorov/CTCF/CTCF - R package version 0.99.4},
#>   }
```

Please note that the `CTCF` was only made possible thanks to many other
R and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.13/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://mdozmorov.github.io/CTCF) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.13/biocthis)*.

## Code of Conduct

Please note that the CTCF project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
