# psupplementary

Hello! Hopefully you're here because you're interested in the R package `psupertime`. `psupertime` is an R package which analyses single cell RNA-seq ("scRNAseq") data where groups of the cells have labels following a known or expected sequence (for example, samples from a time series experiment *day1*, *day2*, ..., *day5*). It uses *ordinal logistic regression* to identify a small set of genes which recapitulate the group-level sequence for individual cells. It can be used for discovery of relevant genes, for exploration of unlabelled data, and assessment of one dataset with respect to the labels known for another dataset. You can find the `psupertime` package [here](https://github.com/wmacnair/psupertime) and read the pre-print [here](https://www.biorxiv.org/content/10.1101/622001v1).

`psupplementary` is a package for replicating the analyses in the `psupertime` paper, and allowing users to play with `psupertime` themselves and see what it can do. The `psupertime` package has everything you need to do your own analysis; splitting the heavy datasets off into the `psupplementary` package keeps the main package light. 


## How to install / use

If you haven't already, you'll need to install `psupertime`, as follows:
```R
remotes::install_github('wmacnair/psupertime', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
library('psupertime')
```
(You may need to install the package `remotes`, with `install.packages('remotes')`. Installation took <90s on a Macbook Pro.)

Due to the large files, installing the `psupplementary` package is slightly more complicated. You need to first clone the package
```sh
cd /path/to/packages
git clone https://github.com/wmacnair/psupplementary.git
```
then you have two options.

*Fast option (~7 minutes):* run this line in R to install it without building the vignettes:
```R
devtools::install('path/to/packages/psupplementary')
```
*Slower option (~15 minutes):* run this line in R to install it with the vignettes:
```R
devtools::install('path/to/packages/psupplementary', build_vignettes=TRUE)
```
This gives you a couple of webpages which walk you through some of the analyses done in the paper.

Once installed, you can call `library('psupplementary')` to load up the package.

## Datasets included

There are six scRNAseq datasets included in this package:

- *acinar_sce*, consisting of aging acinar cells, with sequential labels corresponding to age of donor in years, stored in the variable 'donor_age' (from [here](https://www.sciencedirect.com/science/article/pii/S009286741731053X), GSE81547)
- *germ_sce*, consisting of developing human female germline cells, with sequential labels corresponding to age in weeks, stored in the variable 'time' (from [here](https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(17)30078-4), GSE86146)
- *beta_sce*, consisting of developing beta cells, with sequential labels corresponding to developmental stage, stored in the variable 'age' (from [here](https://www.cell.com/cell-metabolism/fulltext/S1550-4131(17)30208-5), GSE87375)
- *hesc_sce*, consisting of human embryonic stem cells, with sequential labels corresponding to embryonic day, stored in the variable 'esc_day' (from [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3929/), E-MTAB-3929)
- *mef_sce*, consisting of MEFs reprogrammed to neurons, with sequential labels corresponding to days since induction, stored in the variable 'time_point' (from [here](https://www.nature.com/articles/nature18323), GSE67310)
- *colon_sce*, consisting of human colon cells, where the sequential labels are derived from unsupervised clustering (from [here](https://www.cell.com/cell-systems/fulltext/S2405-4712(17)30449-0), GSE102698)

To load `this_sce`, just call `data(this_sce)`.

## Vignettes for replicating analyses

There are several vignettes included in this package:

- _Analysis of acinar cells labelled with donor ages_ (replicates Fig 1C, Fig 1D, Fig 1E, Supp Fig 01, Supp Fig 02, Supp Fig 03, Supp Fig 04, Supp Fig 05)
- _Exploratory data analysis of unlabelled colon data_ (replicates Fig 1F, Supp Fig 15, Supp Fig 16, Supp Fig 17, Supp Fig 18)

(Our intention in future is to allow replication of all figures in the manuscript. At present the code allows replication of all analysis relating to Figure 1, including multiple supplementary figures; replication of the remaining figures will follow!)

To view the vignettes, run this code:
```R
browseVignettes(package = 'psupplementary')
```

Probably the best thing to do is to run some of the analyses in the vignettes, then make copies of the code yourself to explore the possibilities of `psupertime` further.


## Suggestions

Please add any issues or requests to the _Issues_ page. All feedback enthusiastically received.

Thanks!

Will
