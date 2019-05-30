#python3

import sys
import argparse
from random import seed, randint

""" Substring (or Substring Pair) Generator
    Generates a set of random substrings given a circular string.
    Output to stdout:
        Number of lines on the first line
        Each subsequent line follows this formatting:
            Paired - "substring1|substring2|mean_distance"
            Unpaired - "substring"
"""


def read_args():
    parser = argparse.ArgumentParser(
        description="Generates random substrings of fixed length from a circular string \
            with the option for paired-substrings drawn uniformly from an average distance away. \
            Introduces a substitution error with 1% chance. \
            Takes in a string from standard input. \
            Outputs to standard output the number of lines in the first line and \
            substrings (\"ss\") or substring-pairs (\"ss1|ss2|mean_dist\") in subsequent lines.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('length', type=int, help="length of generated substrings")
    parser.add_argument('num_reads', type=int, help="number of single substrings or paired-substrings to generate")
    parser.add_argument('-p', '--paired', action='store_true',
                        help="allows for paired substrings to be drawn from a distance away")
    parser.add_argument('-d', '--distance', type=int, nargs='?',
                        help="mean distance to draw paired substring from", default=125)
    parser.add_argument('-e', '--error', type=int, nargs='?',
                        help="maximum error around mean distance to paired substring, i.e. distance +/- delta", default=0)
   
    return parser.parse_args()


def read_string_from_stdin():
    return sys.stdin.readline().strip()


def generate_reads(string, read_len=100, num_reads=34000, paired=False, d=0, delta=3):
    seed()
    BASES = ['A', 'T', 'C', 'G']
    reads = list()
    for _ in range(num_reads):
        start_pos = randint(0, len(string))
        last_pos = start_pos+read_len
        read = string[start_pos:min(last_pos, len(string))]
        if last_pos > len(string):
            read += string[0:last_pos % len(string)]

        # Introduce a substitution error with a 1% chance
        if (randint(0, 100//read_len - 1) == 0):
            read = list(read)
            read[randint(0, len(read) - 1)] = BASES[randint(0, 3)]
            read = "".join(read)

        assert len(read) == read_len, "STRING GENERATED IS SHORT"

        if not paired:
            reads.append(read)
        else:
            paired_start = (
                start_pos + d + randint(-delta, delta)) % len(string)
            paired_last_pos = paired_start + read_len
            paired_read = string[paired_start:min(
                paired_last_pos, len(string))]
            if paired_last_pos > len(string):
                paired_read += string[0:paired_last_pos % len(string)]

            assert len(paired_read) == read_len, "SUBSTRING GENERATED IS SHORT"
            reads.append("|".join([read, paired_read, str(d)]))
    return reads


def output_reads(reads):
    print(str(len(reads)))
    for read in reads:
        print(read)


def main():
    args = read_args()
    string = read_string_from_stdin()
    reads = generate_reads(string, paired=args.paired,
                           read_len=args.length, num_reads=args.num_reads, d=args.distance, delta=args.error)
    output_reads(reads)


if __name__ == "__main__":
    main()
