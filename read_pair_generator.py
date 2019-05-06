#python3

import sys
from random import seed, randint

""" Read (or Read Pair) Generator
    Generates a set of random reads given a genome.
"""


def read_input():
    return sys.stdin.readline().strip()


def generate_reads(genome, read_len=100, num_reads=34000, paired=False, d=0, delta=3):
    seed()
    reads = list()
    for _ in range(num_reads):
        start_pos = randint(0, len(genome))
        last_pos = start_pos+read_len
        read = genome[start_pos:min(last_pos, len(genome))]
        if last_pos > len(genome):
            read += genome[0:last_pos % len(genome)]

        assert len(read) == read_len, "READ GENERATED IS SHORT"

        if not paired:
            reads.append(read)
        else:
            paired_start = (
                start_pos + d + randint(-delta, delta)) % len(genome)
            paired_last_pos = paired_start + read_len
            paired_read = genome[paired_start:min(
                paired_last_pos, len(genome))]
            if paired_last_pos > len(genome):
                paired_read += genome[0:paired_last_pos % len(genome)]

            assert len(paired_read) == read_len, "READ GENERATED IS SHORT"
            reads.append("|".join([read, paired_read, str(d)]))
    return reads


def output_reads(reads):
    print(str(len(reads)))
    for read in reads:
        print(read)


def main():
    genome = read_input()
    reads = generate_reads(genome, paired=True,
                           read_len=100, num_reads=34000, d=400, delta=3)
    output_reads(reads)


if __name__ == "__main__":
    main()

'''
Test string (length 30):
ACTGGGAAACCCTGGCCATGCCAGTACGAG
'''
