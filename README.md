
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CTCF

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Genomic coordinates of CTCF binding sites with motif MA0139.1 (Jaspar),
in BED format, with strand orientation (directionality of binding).
Human (hg19, hg38) and mouse (mm9, mm10) genomes. The binding sites were
detected using the FIMO tool of the MEME suite using default settings.
Extra columns include motif name (MA0139.1), score, p-value, q-value,
and the motif sequence.

Experimental, to be submitted as an AnnotatiohHub Bioconductor package.

## Installation instructions

<!--
Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `CTCF` using from [Bioconductor](http://bioconductor.org/) the following code:


```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("CTCF")
```

And the development version from [GitHub](https://github.com/mdozmorov/CTCF) with:
-->

``` r
BiocManager::install("mdozmorov/CTCF")
```

## Example

``` r
# hg38 CTCF coordinates
download.file(url = "https://drive.google.com/uc?export=download&id=1s70dQy_zSSP8rsMa9VPEQS1DvUF-bA8o", destfile = "CTCF_hg38.RData")
load(file = "CTCF_hg38.RData")
CTCF_hg38
```

    > CTCF_hg38
    GRanges object with 56049 ranges and 5 metadata columns:
              seqnames            ranges strand |       motif     score       p.value   q.value            sequence
                 <Rle>         <IRanges>  <Rle> | <character> <numeric>     <numeric> <numeric>         <character>
          [1]     chr1       11223-11241      - |    MA0139.1   24.4754 0.00000000134    0.0216 TCGCCAGCAGGGGGCGCCC
          [2]     chr1       11281-11299      - |    MA0139.1   22.7377 0.00000001010    0.0398 GCGCCAGCAGGGGGCGCTG
          [3]     chr1       24782-24800      - |    MA0139.1   17.3770 0.00000071100    0.2350 CGTCCAGCAGATGGCGGAT
          [4]     chr1       91420-91438      + |    MA0139.1   16.2951 0.00000141000    0.3080 GTGGCACCAGGTGGCAGCA
          [5]     chr1     104985-105003      - |    MA0139.1   16.7869 0.00000104000    0.2750 CCAACAGCAGGTGGCAGCC
          ...      ...               ...    ... .         ...       ...           ...       ...                 ...
      [56045]     chrY 57044316-57044334      - |    MA0139.1   16.4590 0.00000127000    0.2990 TGGTCACCTGGGGGCACTA
      [56046]     chrY 57189659-57189677      + |    MA0139.1   15.7541 0.00000195000    0.3430 TGTCCTCTAGGGGTCAGCC
      [56047]     chrY 57203409-57203427      - |    MA0139.1   15.6393 0.00000209000    0.3510 CTGCCGCAAGGGGGCGCAT
      [56048]     chrY 57215279-57215297      + |    MA0139.1   19.5738 0.00000015300    0.1190 gcgccacgagggggcggtg
      [56049]     chrY 57215337-57215355      + |    MA0139.1   24.4754 0.00000000134    0.0216 tcgccagcagggggcgccc

Download the full data from the [Google Drive
folder](https://drive.google.com/drive/folders/19ZXr7IETfks0OdYlmuc1Hqe700Pw3jPc?usp=sharing)

See [inst/scripts/make-data.R](inst/scripts/make-data.R) how to create
the CTCF GRanges objects.

<!--
This is a basic example which shows you how to solve a common problem:


```r
library("CTCF")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
-->

## Citation

Below is the citation output from using `citation('CTCF')` in R. Please
run this yourself to check for any updates on how to cite **CTCF**.

``` r
print(citation("CTCF"), bibtex = TRUE)
#> 
#> Mikhail, Dozmorov, mikhail.dozmorov@gmail.com (2021). _CTCF_.
#> https://github.com/mdozmorov/CTCF/CTCF - R package version 0.0.1, <URL:
#> https://github.com/mdozmorov/CTCF>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {CTCF},
#>     author = {{Mikhail} and {Dozmorov} and {mikhail.dozmorov@gmail.com}},
#>     year = {2021},
#>     url = {https://github.com/mdozmorov/CTCF},
#>     note = {https://github.com/mdozmorov/CTCF/CTCF - R package version 0.0.1},
#>   }
```

Please note that the `CTCF` was only made possible thanks to many other
R and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `CTCF` project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.12/BiocCheck)*.
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
*[biocthis](https://bioconductor.org/packages/3.12/biocthis)*.

## Code of Conduct

Please note that the CTCF project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
