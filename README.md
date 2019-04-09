# psupplementary

Hello! Hopefully you're here because you're interested in the R package `psupertime`. `psupertime` is an R package which uses single cell RNAseq data, where the cells have labels following a known sequence (e.g. a time series), to identify a small number of genes which place cells in that known order. It can be used for discovery of relevant genes, for exploration of unlabelled data, and assessment of one dataset with respect to the labels known for another dataset. You can find `psupertime` [here](https://github.com/wmacnair/psupertime).

`psupplementary` is a package for replicating the analyses in the `psupertime` paper, and allowing users to play with `psupertime` themselves and see what it can do. The `psupertime` package has everything you need to do your own analysis; splitting the heavy data off into this supplementary package keeps the main package light.

## How to install / use

If you haven't already, you'll need to install `psupertime`, as follows:
```R
devtools::install_github('wmacnair/psupertime')
```

Then install the supplementary package, run the following lines in R:
```R
devtools::install_github('wmacnair/psupplementary')
library('psupplementary')
```

## Replicating analyses in the manuscript

To see what you can do, have a look at the vignettes available:
```R
browseVignettes(package = 'psupplementary')
```

Probably the best thing to do is to run some of the analyses in the vignettes, then make copies of the code yourself to explore the possibilities of `psupertime` further.

## Suggestions

Please add any issues or requests to the _Issues_ page. All feedback enthusiastically received.

Thanks!

Will
