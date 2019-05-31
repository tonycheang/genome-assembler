# &#129516; Genome Assembler

> Given a set of error-prone short strings (reads) drawn from a long unknown parent string (genome), we would like to reconstruct the parent string as accurately as possible.

This assembler handles circular strings, characteristic of bacteria genomes.  It takes in either read-pairs (two reads known to be a distance apart) or non-paired reads and outputs contigs (stretches of the parent string with consensus).

## :world_map: Table of Contents

* [Installation and Contents](#-installation-and-contents)
* [Usage](#-usage)
    * [Examples](#examples)
* [Assessment of Assemblies](#bar_chart-assessment-of-assemblies)
* [Main Programs](#bookmark_tabs-main-programs)
    * [generate_reads.py](#generate_readspy)
    * [assemble.py](#assemblepy)
* [Acknowledgements](#pray-acknowledgements)

## &#128230; Installation and Contents

Download the repository. Run the scripts with Python 3.7.

The repository contains the assembler, three genomes drawn from NCBI, and a substring generator to simulate input for the assembler. The example outputs provided have been assessed [here](http://quast.bioinf.spbau.ru/reports/30\_May\_2019\_23:09:09\_000000/report.html).

## &#128736; Usage

First-time assembly of a genome will not have a reference to draw from. However, for our purposes—to see the assembler in action—it’s easiest to use both the read-generator along with the assembler.

The read-generator takes a string from standard input and outputs the number of lines and a list of either unpaired-reads or read-pairs to standard output. The assembler will write its output to a file in ./output or to the standard output if using the stdout flag.

### Examples

The following command generates 19,000 paired-reads, each read of length 100 from the genome of N. Deltocephalinicola, breaks them into substrings of length 28 (kmers) for reconstruction, and filters out infrequent kmers, which contain errors.

```Bash
cat ./reference_genomes/n_deltocephalinicola.txt |
> python generate_reads.py --paired 100 19000 |
> python assemble.py --kmer_length 28 --filter 3
```

The results by default are written to a file in ./output titled with the starting timestamp and constants used in reconstruction.

> May\_30\_09:57:10\_2019\_k6\_f3\_e2.FASTQ

```Bash
>Time started: Thu May 30 16:05:05 2019
>Number of contigs: 15
>CONTIG1
TTAAAAAACATAAATTTAATTATATTTAAAAAAAATAAAA
...
```

> **_Note:_** This assembly will take a few seconds. Assembling E. coli from 700,000 read-pairs of length 100 each read takes about 10 minutes and 6.5 gigabytes of memory.

Short strings can also be used for fun or to test the behavior of the program quickly.

```Bash
echo "It_was_many_and_many_a_year_ago_" |
> python generate_reads.py --paired 8 60 |
> python assemble.py --stdout --kmer_length 6

>Time started:  Wed May 29 21:24:40 2019
>Number of contigs:  1
>CONTIG1
o_It_was_many_and_many_a_year_ag
>Time finished:  Wed May 29 21:24:40 2019
```

Unpaired-reads can be used with the assembler too.

```Bash
$ echo "It_was_many_and_many_a_year_ago_" |
> python generate_reads.py 8 120 |
> python assemble.py --stdout --kmer_length 6

>Time started:  Wed May 29 21:48:57 2019
>Number of contigs:  3
>CONTIG1
nd_many
>CONTIG2
_year_ago_It_was_many
>CONTIG3
_a
>Time finished:  Wed May 29 21:48:57 2019
```

> Resolving repeats can pose a problem if the kmer length doesn't bridge repeated regions.

For best results, the product between the read length and the number of reads should be about 30 for best results (halved for read-pairs, i.e. average coverage depth per character should be 30). The following attempt at reconstructing a string of length 24 only uses a coverage depth of 10.

```Bash
$ echo "In_a_kingdom_by_the_sea_" |
> python generate_reads.py 10 24 |
> python assemble.py --stdout --kmer_length 6 --filter 2

>Time started:  Wed May 29 21:43:23 2019
>Number of contigs:  1
>CONTIG1
ingdom_by_the_sea_In
>Time finished:  Wed May 29 21:43:23 2019
```

>With less coverage, substitution errors in the reads may not be handled properly, leading to erroneous contigs or less-complete reconstruction.

<a href="#-genome-assembler">top</a>

## :bar_chart: Assessment of Assemblies

Output files can be assessed using [QUAST](http://quast.bioinf.spbau.ru/ "Quality Assessment Tool for Genome Assembly").

To assess, first add the files to the list of assemblies.

<p align="center">
<img src="https://user-images.githubusercontent.com/18232816/58675344-c3d4f000-8308-11e9-9bcd-a169cf328c8f.png" style="max-height:200px" title="Upload">
</p>

Then check 'another genome' and upload the corresponding reference (.fasta file in ./reference_genomes).

<p align="center">
<img src="https://user-images.githubusercontent.com/18232816/58675346-c6374a00-8308-11e9-8812-77462590ad0d.png" style="max-height:150px" title="Add a Reference Genome">
</p>

A successful assessment will produce the following table.

![Table of statistics](https://user-images.githubusercontent.com/18232816/58674624-828f1100-8305-11e9-865d-591f5a081516.png)

>The NGA50 value boxed in blue is a typical measure of quality of reconstruction. This value may vary depending on the input data, the parameters of assembly, and the complexity of the parent genome.

<a href="#-genome-assembler">top</a>

## :bookmark_tabs: Main Programs

### generate\_reads.py

Generates random substrings of fixed length from a circular string with the option for paired-substrings drawn uniformly from an average distance away. Introduces a substitution error with 1% chance. Takes in a string from standard input. Outputs to standard output the number of lines in the first line and substrings ("ss") or substring-pairs ("ss1|ss2|mean\_dist”) in subsequent lines.

```Bash
# Sample input:
$ echo "That_a_maiden_there_lived_whom_you_may_know_" |
> python generate_reads.py --paired --distance 4 10 5

# Sample output:
5
_That_a_ma|t_a_maiden|4
here_lived|_lived_who|4
y_know_Tha|ow_That_a_|4
now_That_a|That_a_mai|4
you_may_kn|may_know_T|4
```

More information can be found using the help flag.

```Bash
$ python generate_reads.py -h
usage: generate_reads.py [-h] [-p] [-d [DISTANCE]] [-e [ERROR]]
                         length num_reads
...
```

<a href="#-genome-assembler">top</a>

### assemble.py

Generates contigs of a parent string given a set of substrings or paired substrings. Takes in the number of reads or read-pairs in the first line, then a read ("read") or read-pair ("read1|read2|mean\_dist") on each subsequent line. Outputs a file containing contigs to the subdirectory ./output by default.

```Bash
# Sample input:
$ python assemble.py --kmer_length 5 --filter 1
5
nnabe_Lee;|_Lee;By_th|5
By_the_nam|e_name_of_|5
abe_Lee;By|ee;By_the_|5
e_of_Annab|Annabe_Lee|5
e_name_of_|e_of_Annab|5
```

The file created will look like this.

```FASTQ
# May_30_09:57:10_2019_k5_f1_e2.FASTQ
>Time started: Wed May 29 22:13:55 2019
>Number of contigs: 2
>CONTIG1
he_
>CONTIG2
me_of_Annabe_Lee;By
>Time finished: Wed May 29 22:14:04 2019
```

More information can be found using the help flag.

```Bash
$ python assemble.py -h
usage: assemble.py [-h] [-t] [-m] [-c] [-s] [-k KMER_LENGTH]
                   [-f FILTER_THRESHOLD] [-e ERROR]
...
```

<a href="#-genome-assembler">top</a>

## :pray: Acknowledgements

I'd like to thank

* the staff of Coursera's [Data Structures and Algorithms Specialization](https://www.coursera.org/specializations/data-structures-algorithms) for introducing me to this fascinating problem,
* the staff of [EE 372: Data Science for High-Throughput Sequencing](http://data-science-sequencing.github.io/) for making their lecture notes publicly available and readable,
* Lahnemann, Borkhardt, and McHardy for their [paper](https://www.ncbi.nlm.nih.gov/pubmed/26026159) on handling sequencing errors,
* Medvedev et al. for their [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3216098/) clarifying how to create a paired De Bruijn graph
* Ryan Wick for the great [visualization](https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size) of the effects of kmer size,
* the creators of [SPAdes](http://cab.spbu.ru/software/spades/) for helping me understand the shortcomings of this assembler,
* and lastly, my friends for all their moral support the past few months while I bumbled my way through learning.