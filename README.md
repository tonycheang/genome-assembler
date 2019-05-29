# &#129516; Genome Assembler

## Problem

Given a set of error-prone short strings drawn from a long unknown parent string, we would like to reconstruct the parent string as accurately as possible. This assembler outputs contigs—consensus stretches in the parent string—assuming a circular string, characteristic of bacteria genomes. 

## Contents

The repository contains the assembler, three genomes drawn from the NIH, and a substring generator to simulate input for the assembler. The output files (extension FASTQ) can be compared using the NGA50 values from [QUAST], Quality Assessment Tool for Genome Assembly. Example files provided in ./output.

## Installation

Download the repository.