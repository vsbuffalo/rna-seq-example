# An Example Differential Expression Analysis Using RNA-Seq Data

This is an example analysis of RNA-seq data using open source tools R
and Bioconductor. It starts with raw reads downloaded from the Short
Read Archive (SRA), does quality assessment and improvement, mapping,
and analysis. 

## Data

The data for this comes from Zhang et al., 2011, *Genome-wide mapping
of the HY5-mediated genenetworks in Arabidopsis that involve both
transcriptional and post-transcriptional regulation*
(http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2010.04426.x/full). It
is composed of two Arabidopsis thaliana ecotype Col-0 sets, one hy5
mutant and one wild type. The SRA accession is SRX029582.

 - `GEO:GSM613465`: wild type samples
   - `SRR070570`: replicate 1
   - `SRR070571`: replicate 2
 - `GEO:GSM613466`: hy5 mutant samples
   - `SRR070572`: replicate 1
   - `SRR070573`: replicate 2
 
These were sequenced on an Illumina Genome Analyzer after oligo(dT)
selection and random hexamer priming.

### Downloading

Read data was acquired from the SRA via `wget`, and checked with:

    cd data/raw-reads/
    md5sum -c *sra.md5

The TAIR10 Arabidopsis Genome was downloaded with: 

    wget -nd ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/*fas

`-nd` tells `wget` not to mirror the directory structure.

### FASTQ Extraction

FASTQ data was extracted with NCBI's SRA Toolkit's `fastq-dump`, used
with `xargs` and `find` to run this in parallel:

    find . -name "*sra" | xargs -n1 -P4 /share/apps/sratoolkit/fastq-dump

### Process a Genome for GSNAP

    cd data/athaliana-genome
    gmap_build -d athaliana10 -C *.fas
    
The `-C` option tries to extract chromosome information from the FASTA
header.


### Raw Sequence Quality Assessment

First, we assess the quality with `qrqc`; see `raw-read-qa.Rmd`. 

### Scythe

`sycthe` is a 3'-end quality trimmer. Illumina adapters are removed
from the raw reads:

    find data/raw-reads/ -name "*.fastq" | xargs -n1 -I{} basename {} .fastq | xargs -n1 -P4 -I{} /share/apps/scythe/scythe -q sanger -a /share/apps/scythe/solexa_adapters.JNF.fa -o data/improved-reads/{}-trimmed.fastq data/raw-reads/{}.fastq


### Sickle

`sickle` is a tool for trimming low-quality bases off of the 5'-end
and 3'-end of reads.

## Aligning Reads with GSNAP

I use a script (`map-gsnap.sh`) to process each alignment file. The
purpose of using the script is that I can use it with `xargs`. Because
`gsnap` returns results to standard out, we need to wrap the call in a
script to allow this to be parallelizable.

    find data/improved-reads/ -name "*final.fastq" | xargs -n1 -P4 bash map-gsnap.sh

I am using `-N1` with `gsnap`, which allows for novel
splicing. `gsnap` also allows for known splicing junctions to be used.

## Converting SAM files to BAM files

Another task for `xargs`. As before, I use basename to keep just the
unique file name as a key, and run `samtools` with the correct
directory and extension. This allows me to change the extension too.

    find data/alignments/ -name "*sam" | xargs -n1 -I{} basename {} .sam | xargs -n1 -I{} -P4 samtools view -b -S -o data/alignments/{}.bam data/alignments/{}.sam



