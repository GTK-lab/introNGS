---
title: "Performing read alignments using Bowtie2"
teaching: 10
exercises: 20
questions:
- "How are reads aligned to the reference genome?"
objectives:
- "Understand the different options available in *bowtie2*".
- "Be able to perform an alignment of reads obtained to a reference genome." 
keypoints:


---
# Read alignment using *Bowtie2* 
## Choosing between *bowtie* and *bowtie2* 
Two of the most commonly used alignment software packages are *bowtie* and *bowtie2*. Although relatively similar in usage, both software packages are used in different scenarios. This is summarized in the table below. 

> ## Setting up 
>
> For this practical session, we will be using *bowtie2*. Install *bowtie2* on your local machines using the `apt-get` package manager. Although the version on `apt-get` might not be the latest, it will suffice for the purposes of this exercise. 

## Overview of alignment using *bowtie2* 

### Generating the *bowtie2* index file
Before we can use *bowtie2* to align the reads, we first need to provide the *index files*. These files contain information about the reference genome, such as the sequence and the coordinates. Although *bowtie2* contains some pre-compiled index files, we will walk through how to generate the index files. 

The index files are generated using a single command, `bowtie2-build`. We can find out the list of arguments and its usage as follows:

~~~
bowtie2-build -h
~~~
{: .bash}

This will yield the following:
~~~
Bowtie 2 version 2.2.5 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
    reference_in            comma-separated list of files with ref sequences
    bt2_index_base          write bt2 data to files with this dir/basename
*** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***
Options:
    -f                      reference files are Fasta (default)
    -c                      reference sequences given on cmd line (as
                            <reference_in>)
    --large-index           force generated index to be 'large', even if ref
                            has fewer than 4 billion nucleotides
    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
    -p/--packed             use packed strings internally; slower, less memory
    --bmax <int>            max bucket sz for blockwise suffix-array builder
    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
    --dcv <int>             diff-cover period for blockwise (default: 1024)
    --nodc                  disable diff-cover (algorithm becomes quadratic)
    -r/--noref              don't build .3/.4 index files
    -3/--justref            just build .3/.4 index files
    -o/--offrate <int>      SA is sampled every 2^<int> BWT chars (default: 5)
    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
    --seed <int>            seed for random number generator
    -q/--quiet              verbose output (for debugging)
    -h/--help               print detailed description of tool and its options
    --usage                 print this usage message
    --version               print version information and quit
~~~ 
{: .output}

The usage information tells us that `bowtie2-build` requires only 2 arguments: `reference_in` and `bt2_index_base`. `reference_in` basically is the *fasta* file containing the sequence of the entire genome. This has been provided for todays' class, and should be downloaded. `bt2_index_base` is the prefix of the index file. What happens at the end of running `bowtie2-build` is that a series of files ending with `.bt2` (the index files) will be generated. 

> ## A note on index files
>
> Although *bowtie* and *bowtie2* are similar, their index files cannot be used interchangeably. *bowtie* index files end with `.bt` whilst *bowtie2* index files end with `.bt2`. 

> ## Try it
>
> Download the fasta file containing the human genome. Thereafter, create the *bowtie2* index files. This can take sometime, so be patient with it.
{: .challenge}


