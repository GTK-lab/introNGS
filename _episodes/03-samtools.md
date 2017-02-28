---
title: "Samtools"
teaching: 15
exercises: 15
questions:
- "What are Samtools?"
- "What is sam format?"
- "How to "
objectives:
- "Explain SAM and BAM format."
- "Learn basic usage of Samtools."
keypoints:

---

## What is Samtools

Samtools is a powerful tools that primarily concerned with manipulating SAM files. It allow users to get suitable and desirable input for the downstream analysis.

## SAM format

Now let's take a closer look at SAM file.

SAM files consist of two types of lines: headers and alignments.

* Headers begin with @, and provide meta-data regarding the entire alignment file.

* Alignments begin with any character except @, and describe a single alignment of a sequence read against the reference genome.


Each alignment in SAM format consist of 11 fields, and may also be followed by various optional fields. The 11 fields are:

| Field name   |Description    |Example |
| ------------ |---------------|--------|
| QNAME        | Unique identifier of the read derived from the original FASTQ file. | |
| FLAG | A single integer value (e.g. 16), which encodes multiple elements of meta-data regarding a read and its alignment.| |
| RNAME | Reference genome identifier | |
| POS | 	Left-most position within the reference genome where the alignment occurs.(**1-based** ) | |
| MAPQ | Quality of the genome mapping. The MAPQ field uses a Phred-scaled probability value to indicate the probability that the mapping to the genome is incorrect. Higher values indicate increased confidence in the mapping. | |
| CIGAR | A compact string that (partially) summarizes the alignment of the raw sequence read to the reference genome. Three core abbreviations are used: M for alignment match; I for insertion; and D for Deletion.  | |
| RNEXT | Reference genome identifier where the mate pair aligns. Only applicable when processing paired-end sequencing data. A value of * indicates that information is not available.| |
| PNEXT | Position with the reference genome, where the second mate pair aligns. As with RNEXT, this field is only applicable when processing paired-end sequencing data. A value of 0 indicates that information is not available. | |
| TLEN | Template Length. Only applicable for paired-end sequencing data, TLEN is the size of the original DNA or RNA fragment, determined by examining both of the paired-mates, and counting bases from the left-most aligned base to the right-most aligned base. A value of 0 indicates that TLEN information is not available. | |
| SEQ | The raw sequence, as originally defined in the FASTQ file.| |
| QUAL | The Phred quality score for each base, as originally defined in the FASTQ file. | |

### FLAG

The flag field in SAM in a *combination* of bitwise FLAGs. It can provide us with many useful information.

|FLAG| Description|
|---|------------|
|1 or 0x1|Template having multiple segments in sequencing|
|2 or 0x2|Each segment properly aligned according to the aligner|
|4 or 0x4|Segment unmapped|
|8 or 0x8|Next segment in the template unmapped|
|16 or 0x10|SEQ being reverse complemented|
|32 or 0x20|SEQ of the next segment in the template being reversed|
|64 or 0x40|The first segment in the template|
|128 or 0x80|The last segment in the template|
|256 or 0x100|Secondary alignment|
|512 or 0x200|Not passing quality controls|
|1024 or 0x400|PCR or optical duplicate|
|2048 or 0x800|Supplementary alignment|

Those 0x1,0x2... are hexadecimal numbers (base 16) as opposed to decimal numbers (base 10).

If you see a alignment with a FLAG of 163, we can transform it into a binary string: 10100011. Map this string from right to left to the table above, we can get:

|Bit|Description|
|1||Template having multiple segments in sequencing|
|1|Each segment properly aligned according to the aligner|
|0|Segment unmapped|
|0|Next segment in the template unmapped|
|0|SEQ being reverse complemented|
|1|SEQ of the next segment in the template being reversed|
|0|The first segment in the template|
|1|The last segment in the template|
|0|Secondary alignment|
|0|Not passing quality controls|
|0|PCR or optical duplicate|
|0|Supplementary alignment|

Consider 1 as TRUE and 0 as FALSE, it can tell us that this alignment:
1. comes from a pair-end sequencing data
2. is mapped in a proper pair
3. has a mate read that is mapped reversely
4. is the second read in the pair


## Converting SAM to BAM

To do anything meaningful with alignment data from Bowtie2 or other aligners (which produce text-based SAM output), we need to first convert the SAM to its binary counterpart, BAM format. The binary format is much easier for computer programs to work with. However, it is consequently very difficult for humans to read.You can do:

~~~
samtools view -b -S -o sample.bam sample.sam
~~~
{: .bash}

* `-b`: Indicate the output to be BAM format
* `-S`: Indicate the input to be SMA format
* `-o`: Specifies the output file name

Or you can also use a pipe line:

~~~
samtools view -b -S sample.sam > sample.bam
~~~
{: .bash}

Even if BAM format is unreadable, we can use view command to display all alignments.

~~~
samtools view sample.bam | head
~~~
{: .bash}

> ## Exercise
>
> Try to converting your SAM file to BAM file
>
>
{: .challenge}

View command of Samtools also allows us to display specified alignments. As we discussed earlier, the FLAG field in the BAM format encodes several key pieces of information regarding how an alignment aligned to the reference genome. We can exploit this information to isolate specific types of alignments that we want to use in our analysis.

~~~
samtools view -f 4 sample.bam | head
~~~
{: .bash}

* `-f`: Extract those alignments that match the specified SAM flag. In the case above, we are only looking for alignments with flag value of 4(read fails to map to the reference genome.)

Now, we can also exclude alignments that match:

~~~
samtools view -F 4 sample.bam | head
~~~
{: .bash}

You can also try out `-c` option, which output number of alignments instead of alignments themselves.

~~~
samtools view -c -f 4 sample.bam | head
~~~
{: .bash}

`-q` option can perform a mapping quality score filter of the alignments.

~~~
samtools view -c -q 42 sample.bam | head
~~~
{: .bash}

This would output the number of alignments that have a mapping quality score of 42 or higher.

> ## Exercise
>
> Count the number of reversely-mapped reads with a quality score cutoff 42 in your BAM file.
>
> Hint: flag 0x10 = SEQ being reverse complemented
>
> > ## Solution
> >
> > samtools view -q 42 -F 0x4 -f 0x10 sample.bam | wc -l
> >
> > or: samtools view -c -q 42 -F 0x4 -f 0x10 sample.bam
> >
> {: .solution}
{: .challenge}


## Sorting a bam file

When you align FASTQ files with all current sequence aligners, the alignments produced are in random order with respect to their position in the reference genome. In other words, the BAM file is in the order that the sequences occurred in the input FASTQ files. You can observe that using Samtools view command.

There are two options for sorting BAM files: by read name (`-n`), and by genomic location (default).

~~~
samtools sort sample.bam sample.sorted.bam
~~~
{: .bash}

Now use the view command again to see the order.

## Indexing

Indexing a genome sorted BAM file allows one to quickly extract alignments overlapping particular genomic regions. Moreover, indexing is required by genome viewers such as IGV so that the viewers can quickly display alignments in each genomic region to which you navigate.

~~~
samtools index sample.sorted.bam
~~~
{: .bash}

This will create a BAM index file with `.bai` extension. Note that if you attempt to index a BAM file which has not been sorted, you will receive an error and indexing will fail.

After indexing, we are allowed to see all alignments with a specific loci.

~~~
samtools view sample.bam chr1:200000-500000
~~~
{: .bash}

> ## Exercise
>
> Sort your BAM file and create an index file of the sorted BAM file. Then, count the number of reversely-mapped reads between chr3:123456 and chr3:124456
>
> > ## Solution
> >
> > samtools view -F 0X4 -f 0X10 sample.sorted.bam chr3:123456-124456 | wc -l
> >
> > or: samtools view -c -F 0X4 -f 0X10 sample.sorted.bam chr3:123456-124456
> >
> {: .solution}
{: .challenge}
