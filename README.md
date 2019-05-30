# &#129516; Genome Assembler

## Problem Statement

Given a set of error-prone short strings (reads) drawn from a long unknown parent string (genome), we would like to reconstruct the parent string as accurately as possible.

This assembler handles circular strings, which is characteristic of bacteria genomes.  It takes in either read-pairs (two reads known to be a distance apart) or non-paired reads and outputs contigs (stretches of the parent string with consensus).

&nbsp;

## :open_file_folder: Contents

The repository contains the assembler, three genomes drawn from NCBI, and a substring generator to simulate input for the assembler. The output files can be assessed using [QUAST](http://quast.bioinf.spbau.ru/ "Quality Assessment Tool for Genome Assembly"). The example outputs provided have been assessed [here](http://quast.bioinf.spbau.ru/reports/30\_May\_2019\_23:09:09\_000000/report.html).

&nbsp;

## :floppy\_disk: Installation

Download the repository. Run the scripts with Python 3.7.

&nbsp;

## :hammer: Usage

First-time assembly of a genome will not have a reference to draw from. However, for our purposes—to see the assembler in action—it’s easiest to use both the read-generator along with the assembler. 

The read-generator takes a string from standard input and outputs the number of lines and a list of either unpaired-reads or read-pairs to standard output. The assembler will write its output to a file in ./output or to the standard output if using the stdout flag.

&nbsp;

### Examples

The following program generates 60 read-pairs, each read of length 8 from the echo’ed string, breaks them into substrings of length 6 (kmer\_length) for assembly, and reconstructs a single contig.

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

Unpaired-reads can also be used with the assembler, but resolving repeats can pose a problem if the kmer length doesn't bridge the repeat.

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

The product between the read length and the number of reads should be about 30 for best results (halved for read-pairs, i.e. average coverage per character should be 30). Otherwise substitution errors in the reads may not be handled properly, leading to erroneous contigs or less-complete reconstruction. The following attempt at reconstructing a string of length 24 only uses a coverage of 10.

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

Finally, these programs can simulate reconstructing a real genome. The results by default are written to a file in ./output with the starting timestamp and constants used in reconstruction (e.g. May\_30\_09:57:10\_2019\_k6\_f3\_e2.FASTQ)

```Bash
cat ./reference_genomes/n_deltocephalinicola.txt |
> python generate_reads.py --paired 100 19000 |
> python assemble.py --kmer_length 28 --filter 3
```

Note that this assembly will take a few seconds. Assembling E. coli from 700,000 read-pairs of length 100 each read takes about 10 minutes and 6.5 gigabytes of memory.

&nbsp;

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

&nbsp;

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