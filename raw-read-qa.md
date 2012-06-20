# Raw Read Quality Assessment

It's extremely important to do quality assessment and improvement
steps before working with sequencing data. I highly recommend this one
take an iteractive approach to this, consisting of diagnostics with
`qrqc`, then improvement, then more diagnostics. This iterative
approach has two benefits:

1. It prevents the dangerous assumption that tools that work well
generally are working well with a particular dataset. 

2. Extremely pathological datasets stick out more in comparative
before/after quality improvement software diagnostics.



```r
opts_chunk$set(fig.width = 7, fig.height = 7)
opts_knit$set(base.url = "https://github.com/vsbuffalo/rna-seq-example/raw/master/")
```




## Load Required Packages and Set Number of Cores to Use



```r
library(mclapply)
```

```
## Error: there is no package called 'mclapply'
```

```r
library(qrqc)

options(mc.cores = 4)
```




## Raw Quality Reports



```r

raw.fastq.files <- list.files("data/raw-reads", pattern = ".*\\.fastq", 
    full.names = TRUE)
names(raw.fastq.files) <- basename(raw.fastq.files)
raw.fastq.summaries <- mclapply(raw.fastq.files, readSeqFile)
```




### Base Quality Reports



```r

baseQualPlot(raw.fastq.summaries)
```

```
## Error: could not find function "baseQualPlot"
```



