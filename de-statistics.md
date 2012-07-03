# Differential Expression Statistics



```r
opts_chunk$set(fig.width = 7, fig.height = 7, cache = TRUE)
opts_knit$set(base.url = "https://github.com/vsbuffalo/rna-seq-example/raw/master/")
```




## Load Required Packages



```r
library(ggplot2)
library(DESeq)
library(IRanges)
```



## Design

Normally, I extract the design from the file names, but in this case,
the file names are just SRA identifiers, so this is done manually:



```r

design <- data.frame(file = c("SRR070570", "SRR070571", "SRR070572", 
    "SRR070573"), treatment = c("wildtype", "wildtype", "mutant", "mutant"), 
    stringsAsFactors = FALSE)
```




## Read in Count Data and Do EDA



```r

d <- read.table("results/raw-counts.txt", header = TRUE, row.names = 1)
```






```r

colSums(d)
```

```
## SRR070570 SRR070571 SRR070572 SRR070573 
##  11728960  11585090  13162350  12267736 
```




### MDS Plot

An MDS plot is a great diagnostic - both on raw counts and normalized
and variance stabilized counts. 



```r
mdsPlot <- function(counts, conds, useVST = FALSE) {
    if (useVST) {
        cds <- newCountDataSet(counts, conds)
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds, method = "blind")
        counts <- getVarianceStabilizedData(cds)
    }
    d <- dist(t(counts))
    
    mds.fit <- cmdscale(d, eig = TRUE, k = 2)
    
    mds.d <- data.frame(x1 = mds.fit$points[, 1], x2 = mds.fit$points[, 2], 
        labels = colnames(counts))
    
    mds.d$treatment <- as.factor(conds)
    
    ggplot(mds.d) + geom_text(aes(x = x1, y = x2, label = labels, colour = treatment))
}

mdsPlot(d, design$treatment[match(colnames(d), design$file)]) + opts(title = "MDS on Raw Counts")
```

![plot of chunk mdsplot](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/mdsplot1.png) 

```r
mdsPlot(d, design$treatment[match(colnames(d), design$file)], useVST = TRUE) + 
    opts(title = "MDS on Normalized and VST Counts")
```

![plot of chunk mdsplot](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/mdsplot2.png) 


What do we look for in MDS plots? Well, as with many diagnostic plots
it takes practice and is subjective. Personally, I like to see that:

 - Samples are linearly seperable (or close) by treatment. If not,
   it's *likely* (again, this is subjective) that there is either too
   much biological heterogeniety which leads to high variability
   across replicates. If all samples clumped, or replicates per
   treatment all over in the MDS plot, this *could* cause problems in
   variance estimation and lead to few differentially expressed genes
   being found (or those that are found are unreliable). Either more
   biology replicates are needed, or further experimental controls are
   needed (or both!).

 - Replicates are somewhat clustered (especially if there are over
   four). Don't over interpret the absolute distances in MDS plots:
   remember, things are scaled to fit the plot window well, which can
   make small actual distances look larger. Consider relative and
   pairwise distances.

In this plot, we see there is more distance between wild type
replicates than mutant replicates.

### Normalization

There are many RNA-seq normalization strategies. The most common
(sadly, but that's just my opinion) is RPKM: reads per kilobase of
exon per million mapped reads. This is a normalization scheme that
attempts to control for varying transcript lengths (longer
transcripts, more reads) and differences in library sizes. One can
think of this as using a global scaling factor: total mapped reads (or
in our example, the column sums).

The problem with this approach is that some genes are extremely
expressed and drive total lane counts. We also know from the
mean-variance relationship of count data that highly expressed genes
also have the highest variance (worsened by overdispersion!), so tiny
(likely) differences in the expression of these highly expressed genes
across replicates can lead to drastic changes in total lane counts and
thus the scaling factor. 

To illustrate this, let's look at the top 2% of highly expressed genes
in `SRR070570`:



```r

q98 <- quantile(d$SRR070570, probs = 0.98)
sum(d$SRR070570[d$SRR070570 >= q98])/sum(d$SRR070570)
```

```
## [1] 0.899
```




So the top 2% of gene's counts make up nearly 90% of total lane
counts. 

Other approaches try to be more robust. DESeq's approach is to use the
geometric mean across samples to build a reference distribution. Then
the median of the genewise ratios of counts to the reference
distribution is used as a per-lane scaling factor. This is used
directly in fitting process.

An MA-plot can be used (with caution - this is not a microarray
analysis!) to look at how things are distributed.



```r

makeMAData <- function(d, conds, a.name, b.name) {
    samples.a <- d[, conds == a.name]
    samples.b <- d[, conds == b.name]
    
    lfc <- log2(rowMeans(samples.a)/rowMeans(samples.b))
    all.zero <- rowMeans(samples.a) == 0 | rowMeans(samples.b) == 0
    lfc[all.zero] <- (log2(rowMeans(samples.a + 1)/rowMeans(samples.b + 1)))[all.zero]
    mean.exp <- rowMeans(d)
    mean.exp[all.zero] <- 0.5
    ma.d <- data.frame(lfc, mean.exp, adjusted = all.zero)
    ma.d
}
ma.d <- makeMAData(d, design$treatment[match(colnames(d), design$file)], 
    "mutant", "wildtype")
p <- ggplot(ma.d) + geom_point(aes(x = mean.exp, y = lfc, color = adjusted), 
    size = 0.8) + scale_x_log10("mean expression") + scale_y_continuous("log2 fold change (mutant/wildtype)")
p <- p + geom_smooth(aes(x = mean.exp, y = lfc), se = FALSE)
p
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using
## gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the
## smoothing method.
```

![plot of chunk ma-plot](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/ma-plot.png) 


We *usually* operator under the null hypothesis that most genes are
not differentially expressed (although, depending on what you're
comparing, this belief should be updated), so we would expect the blue
lowess curve to be close to the y=0 line. 

We can see how much DESeq's normalization procedure would change
things by running the first few steps of a DE analysis and outputing
normalized counts:



```r

cds <- newCountDataSet(d, design$treatment[match(colnames(d), design$file)])
cds <- estimateSizeFactors(cds)
norm.counts <- counts(cds, normalized = TRUE)
ma.d$title <- "Before Normalization"
ma.d.norm <- makeMAData(norm.counts, design$treatment[match(colnames(norm.counts), 
    design$file)], "mutant", "wildtype")
ma.d.norm$title <- "After Normalization"

ma.d.all <- rbind(ma.d, ma.d.norm)

p <- ggplot(ma.d.all) + geom_point(aes(x = mean.exp, y = lfc, color = adjusted), 
    size = 0.8) + scale_x_log10("mean expression") + scale_y_continuous("log2 fold change (mutant/wildtype)")
p <- p + geom_smooth(aes(x = mean.exp, y = lfc), se = FALSE)
p + facet_wrap(~title)
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using
## gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the
## smoothing method.
```

```
## geom_smooth: method="auto" and size of largest group is >=1000, so using
## gam with formula: y ~ s(x, bs = "cs"). Use 'method = x' to change the
## smoothing method.
```

![plot of chunk deseq-norm](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/deseq-norm.png) 


We can see (at least according to DESeq's normalization scheme) that
not much normalization was needed.

### Mean and Variance

Count data is often modeled using the Poisson distribution, which is a
distribution with one parameter, equal to the mean and variance. The
assumption that the mean equals the variance is often not valid: there
is what is known as **overdispersion**, or variance > mean. This is
usually the case when replicates are really more accurately draws from
a mixture model (often the case when you have heterogeneous biological
individuals). Technical replicates often do satisfy the Poisson model
(since the source of error is not heterogeneity, but rather techincal
noise). 



```r

mean.wt <- rowMeans(d[, design$file[design$treatment == "wildtype"]])
mean.mut <- rowMeans(d[, design$file[design$treatment == "mutant"]])
var.wt <- apply(d[, design$file[design$treatment == "wildtype"]], 
    1, var)
var.mut <- apply(d[, design$file[design$treatment == "mutant"]], 
    1, var)

mv.d <- data.frame(mean = c(mean.wt, mean.mut), variance = c(var.wt, 
    var.mut), sample = as.vector(Rle(c("wildtype", "mutant"), c(nrow(d), nrow(d)))))

ggplot(mv.d) + geom_point(aes(x = mean, y = variance, color = sample), 
    size = 1) + scale_x_log10() + scale_y_log10() + geom_abline(slope = 1, intercept = 0, 
    color = "blue")
```

![plot of chunk mean-var](https://github.com/vsbuffalo/rna-seq-example/raw/master/figure/mean-var.png) 


## Differential Expressio Analysis

We can conduct an initial differential expression analysis with
`estimateDispersions` and `nbinomTest`:



```r

cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "mutant", "wildtype")
```



