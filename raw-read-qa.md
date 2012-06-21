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
opts_chunk$set(fig.width = 7, fig.height = 7, cache = TRUE)
opts_knit$set(base.url = "https://github.com/vsbuffalo/rna-seq-example/raw/master/")
```




## Load Required Packages and Set Number of Cores to Use



```r
library(reshape)
library(ggplot2)
library(multicore)
library(qrqc)

options(mc.cores = 4)
```




## Raw Quality Reports



```r

raw.fastq.files <- list.files("data/raw-reads", pattern = ".*\\.fastq", 
    full.names = TRUE)
names(raw.fastq.files) <- basename(raw.fastq.files)

raw.fastq.summaries <- mclapply(raw.fastq.files, readSeqFile)

# Add in random reads (for comparison). This also helps with scaling the y
# axis
raw.fastq.summaries[["random"]] <- readSeqFile(system.file("extdata", 
    "random.fasta", package = "qrqc"), type = "fasta", hash.prop = 1)
```

```
## Warning: Some k-mer counts are infinite, meaning there was occurence of a
## k-mer over the maximum double size.
```




### Base Quality



```r

# omit random FASTQ file
qualPlot(raw.fastq.summaries[-length(raw.fastq.files)])
```

```
## Error: All items in list must have class FASTQSummary.
```




### Base Frequency



```r

basePlot(raw.fastq.summaries, type = "proportion")
```

![plot of chunk raw-base-frequency](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/raw-base-frequency.png) 


### K-mer Contaminant Plots



```r

kmerKLPlot(raw.fastq.summaries)
```

```
## Error: replacement has 1 rows, data has 0
```




### Entropy Contaminant Plots



```r

kmerEntropyPlot(raw.fastq.summaries)
```

![plot of chunk raw-entropy](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/raw-entropy.png) 


## Statistic of Processed Reads (after Sickle and Scythe)

These statistics were gathered from `scythe` and `sickle` output. This
provides a look how many reads were removed, the distribution of the
adapter contaminants.




```r

d <- read.table("data/scythe-data.txt", sep = "\t", header = TRUE)
```




This is a bit of a hodgepodge of a data file; it has original read
counts, counts at each stage, and a character column of Scythe output
indicating where the contaminants were found. First, we extract this
last column. I plan on changing Scythe soon (and perhaps Sickle) so
that output is nicer for downstream statistics.



```r
d.scythe <- local({
    # look at just adapter columns, and the id (file) col
    tmp <- d[, c(1, grep("adapter", colnames(d))), ]
    
    # split out and convert the values to numeric
    tmp.a1 <- lapply(strsplit(as.character(tmp[, 2]), ", "), as.numeric)
    tmp.a2 <- lapply(strsplit(as.character(tmp[, 3]), ", "), as.numeric)
    
    # append a position vector
    tmp.a1[[length(tmp.a1) + 1]] <- seq_along(tmp.a1[[1]])
    tmp.a2[[length(tmp.a2) + 1]] <- seq_along(tmp.a2[[1]])
    
    # create a long dataframe of the above data for each adapter, and rbind
    # both together
    a1 <- data.frame(adapter = "adapter 1", do.call(cbind, tmp.a1))
    colnames(a1)[2:6] <- c(as.character(tmp[, 1]), "position")
    a2 <- data.frame(adapter = "adapter 2", do.call(cbind, tmp.a2))
    colnames(a2)[2:6] <- c(as.character(tmp[, 1]), "position")
    rbind(a1, a2)
})

d.scythe <- melt(d.scythe, id.vars = c("adapter", "position"))

p <- ggplot(d.scythe) + geom_bar(aes(x = position, y = value, fill = adapter), 
    position = "dodge", stat = "identity")
p <- p + scale_y_continuous("count")
p
```

![plot of chunk remove-contaminant-col](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/remove-contaminant-col.png) 


Now we can look at how many reads Scythe found to be contaminated:



```r

d.scythe.trimmed <- melt(d[, c("file", "total", "uncontaminated")], 
    id.vars = c("file"))

p <- ggplot(d.scythe.trimmed) + geom_bar(aes(x = file, y = value, 
    fill = variable))
p <- p + scale_y_continuous("count")
p
```

![plot of chunk unnamed-chunk-1](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/unnamed-chunk-1.png) 


And how many reads Sickle trimmed:



```r

d.sickle <- melt(d[, c("file", "kept", "discarded")], id.vars = c("file"))

p <- ggplot(d.sickle) + geom_bar(aes(x = file, y = value, fill = variable))
p <- p + scale_y_continuous("count")
p
```

![plot of chunk unnamed-chunk-2](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/unnamed-chunk-2.png) 


## Post-Processed Quality Reports



```r
o
```

```
## Error: object 'o' not found
```

```r
processed.fastq.files <- list.files("data/improved-reads", pattern = ".*final\\.fastq", 
    full.names = TRUE)
names(processed.fastq.files) <- basename(processed.fastq.files)

processed.fastq.summaries <- mclapply(processed.fastq.files, readSeqFile)

# Add in random reads (for comparison). This also helps with scaling the y
# axis
processed.fastq.summaries[["random"]] <- readSeqFile(system.file("extdata", 
    "random.fasta", package = "qrqc"), type = "fasta", hash.prop = 1)
```




### Base Quality



```r

qualPlot(processed.fastq.summaries[-length(raw.fastq.files)])
```

```
## Error: All items in list must have class FASTQSummary.
```




### Base Frequency



```r

basePlot(processed.fastq.summaries, type = "proportion")
```

![plot of chunk processed-base-frequency](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/processed-base-frequency.png) 


### K-mer Contaminant Plots



```r

kmerKLPlot(processed.fastq.summaries)
```

```
## Warning: Stacking not well defined when ymin != 0
```

```
## Warning: Stacking not well defined when ymin != 0
```

```
## Warning: Stacking not well defined when ymin != 0
```

```
## Warning: Stacking not well defined when ymin != 0
```

![plot of chunk processed-kmer-kl](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/processed-kmer-kl.png) 


### Entropy Contaminant Plots



```r

kmerEntropyPlot(processed.fastq.summaries)
```

![plot of chunk processed-entropy](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/processed-entropy.png) 



### Sequence Length Plot



```r

seqlenPlot(processed.fastq.summaries)
```

![plot of chunk processed-length](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/processed-length.png) 


