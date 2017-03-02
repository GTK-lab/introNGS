---
title: "Read alignment using bowtie2"
teaching: 15
exercises: 30
questions:
- "How do we perform alignment of reads, starting from a fastq file?"
objectives:
- "Understand the difference between bowtie and bowtie2, and choose the appropriate software to perform alignment."
- "Be able to align raw reads, starting from a raw fastq file." 
- "Understand the range of options available in bowtie2 for alignment." 
keypoints:

---

# Read alignment using *Bowtie2* 
## Choosing between *bowtie* and *bowtie2* 
Two of the most commonly used alignment software packages are *bowtie* and *bowtie2*. Although relatively similar in usage, both software packages are used in different scenarios. This is summarized in the table below. 

> ## Setting up 
>
> For this practical session, we will be using *bowtie2*. Install *bowtie2* on your local machines using the `apt-get` package manager. Although the version on `apt-get` might not be the latest, it will suffice for the purposes of this exercise. 

## Overview of alignment using *bowtie2* 

## Generating the *bowtie2* index file
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

## Aligning reads to the reference genome
While `bowtie2-build` works hard to build the index files (this can take up to 15 minutes), we will discuss the usage of `bowtie2` which is the workhorse of alignment. Typing `bowtie2 -h` will yield a **very long** page full of arguments that we can provide to `bowtie2` when we perform the alignment. 

~~~
jeremy@atlas:~$ bowtie2
No index, query, or output file specified!
Bowtie 2 version 2.2.5 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
Usage: 
  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <sam>      File for SAM output (default: stdout)

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.

Options (defaults in parentheses):

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
  --qseq             query input files are in Illumina's qseq format
  -f                 query input files are (multi-)FASTA .fa/.mfa
  -r                 query input files are raw one-sequence-per-line
  -c                 <m1>, <m2>, <r> are sequences themselves, not files
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
  --phred33          qualities are Phred+33 (default)
  --phred64          qualities are Phred+64
  --int-quals        qualities encoded as space-delimited integers

 Presets:                 Same as:
  For --end-to-end:
   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

  For --local:
   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

 Alignment:
  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
  -L <int>           length of seed substrings; must be >3, <32 (22)
  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
  --dpad <int>       include <int> extra ref chars on sides of DP table (15)
  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
  --ignore-quals     treat all quality values as 30 on Phred scale (off)
  --nofw             do not align forward (original) version of read (off)
  --norc             do not align reverse-complement version of read (off)
  --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to
                     scan for the optimal seeded alignments
  --end-to-end       entire read must align; no clipping (on)
   OR
  --local            local alignment; ends might be soft clipped (off)

 Scoring:
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
  --rdg <int>,<int>  read gap open, extend penalties (5,3)
  --rfg <int>,<int>  reference gap open, extend penalties (5,3)
  --score-min <func> min acceptable alignment score w/r/t read length
                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)

 Reporting:
  (default)          look for multiple alignments, report best, with MAPQ
   OR
  -k <int>           report up to <int> alns per read; MAPQ not meaningful
   OR
  -a/--all           report all alignments; very slow, MAPQ not meaningful

 Effort:
  -D <int>           give up extending after <int> failed extends in a row (15)
  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)

 Paired-end:
  -I/--minins <int>  minimum fragment length (0)
  -X/--maxins <int>  maximum fragment length (500)
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
  --no-mixed         suppress unpaired alignments for paired reads
  --no-discordant    suppress discordant alignments for paired reads
  --no-dovetail      not concordant when mates extend past each other
  --no-contain       not concordant when one mate alignment contains other
  --no-overlap       not concordant when mates overlap at all

 Output:
  -t/--time          print wall-clock time taken by search phases
  --un <path>           write unpaired reads that didn't align to <path>
  --al <path>           write unpaired reads that aligned at least once to <path>
  --un-conc <path>      write pairs that didn't align concordantly to <path>
  --al-conc <path>      write pairs that aligned concordantly at least once to <path>
  (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
  --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
  --quiet            print nothing to stderr except serious errors
  --met-file <path>  send metrics to file at <path> (off)
  --met-stderr       send metrics to stderr (off)
  --met <int>        report internal counters & metrics every <int> secs (1)
  --no-unal          supppress SAM records for unaligned reads
  --no-head          supppress header lines, i.e. lines starting with @
  --no-sq            supppress @SQ header lines
  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                     Note: @RG line only printed when --rg-id is set.
  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.

 Performance:
  -p/--threads <int> number of alignment threads to launch (1)
  --reorder          force SAM output order to match order of input reads
  --mm               use memory-mapped I/O for index; many 'bowtie's can share

 Other:
  --qc-filter        filter out reads that are bad according to QSEQ filter
  --seed <int>       seed for random number generator (0)
  --non-deterministic seed rand. gen. arbitrarily instead of using read attributes
  --version          print version information and quit
  -h/--help          print this usage message
~~~
{: .output}

Now that is some seriously long list of options! Thankfully, we will not be required to provide parameters for all the options. In fact, for most practical purposes, we use the defaults for most of the arguments. 

Lets start by looking at the usage. The usage of `bowtie2` is simple, and takes the following form: 

~~~
bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
~~~
{: .bash}

Remember the index files we generated earlier on? Here, we need to tell `bowtie2` where these files are located, and the prefix of the files. In this way, we can (in theory) have the index files of many different genomes stored in the same location. This is specified in the `-x` argument. The next set of arguments `-1 <m1> -2 <m2>| -R <r>` basically tell `bowtie2` where the fastq files are. As you will recall, NGS can yield either paired-end reads or unpaired-end reads. If we are working with paired-end reads, the results will be returned in 2 files (one for each strand). If such is the case, we will use `-1 <name> -2 <name>`. However, for unpaired end reads, we will only need to use `-U <name>`. Finally, the last option `-S <sam>` specifies the output SAM files. This needs to be specified explicitly. 

> ## Find out
>
> What happens if no SAM file is specified with the `-S` option? 