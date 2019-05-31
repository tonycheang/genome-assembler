import argparse
import sys
import time

from debruijn_graph import DeBruijnGraph, PairedDeBruijnGraph
from debug_graph import DebugDeBruijnGraph, DebugCMSDeBruijnGraph
from debug_graph import DebugPairedDeBruijnGraph, DebugCMSPairedDeBruijnGraph
from pathlib import Path


class IOHandler:
    @staticmethod
    def read_args():
        parser = argparse.ArgumentParser(
            description="Generates contigs (consensus regions) of a parent string \
                given a set of substrings or paired substrings. Takes in the number \
                of reads or read-pairs in the first line, then a read (\"read\") or \
                read-pair (\"read1|read2|mean_dist\") on each subsequent line. Outputs \
                a file containing contigs to the subdirectory ./output by default.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-t', '--time', action='store_true',
                            help="prints time at each stage during runtime to track program progression")
        parser.add_argument('-m', '--memory', action='store_true',
                            help="prints size of major data structures during runtime")
        parser.add_argument('-c', '--count_min_sketch', action='store_true',
                            help="uses a probabilistic data structure instead of a dictionary")
        parser.add_argument('-s', '--stdout', action='store_true',
                            help="switches output to stdout instead of file write")

        parser.add_argument('-k', '--kmer_length', type=int,
                            help="the k-value used to break down reads")
        parser.add_argument('-f', '--filter_threshold', type=int,
                            help="filter threshold for erroneous kmers")
        parser.add_argument('-e', '--error', type=int,
                            help="allowed error in paired distance for reconstruction with paired-reads")

        return parser.parse_args()

    @staticmethod
    def read_input():
        ''' Returns a list of reads or read-pairs, not yet broken down,
            indications of whether the reads are paired, their distance if so,
            and the number of bases total.        
        '''
        n = int(sys.stdin.readline().strip())
        paired = False
        reads = list()
        d = 0

        num_bases = 0

        # Check for read-pair vs non-paired input format
        first_line = sys.stdin.readline().strip().split('|')
        if len(first_line) > 1:
            read1, read2, d = first_line
            num_bases += len(read1) + len(read2)
            reads.append((read1, read2))
            paired = True
            for _ in range(n-1):
                read1, read2, d = sys.stdin.readline().strip().split('|')
                reads.append((read1, read2))
                num_bases += len(read1) + len(read2)
        else:
            reads.append(first_line[0])
            num_bases += len(first_line[0])
            for _ in range(n-1):
                read = sys.stdin.readline().strip()
                reads.append(read)
                num_bases += len(read)

        return (reads, paired, int(d), num_bases)

    @staticmethod
    def write_stdout(contigs, consants, start_time):
        print(">Time started: ", time.strftime("%c", time.localtime(start_time)))
        print(">Number of contigs: ", len(contigs))
        for i, contig in enumerate(contigs):
            print(">CONTIG" + str(i+1))
            print(contig)
        print(">Time finished: ", time.strftime("%c", time.localtime()))

    @staticmethod
    def write_FASTQ(contigs, constants, start_time):
        output_dir = Path("./output")
        output_dir.mkdir(exist_ok=True)

        formatted_start_time = time.strftime("%b_%d_%H:%M:%S_%Y", time.localtime(start_time))
        formatted_constants = "k{0}_f{1}_e{2}".format(*constants)
        filename = formatted_start_time + "_" + formatted_constants + ".FASTQ"

        output_file = Path("./output/" + filename)
        with output_file.open(mode='x') as file:
            file.write(">Time started: " + time.strftime("%c", time.localtime(start_time)) + "\n")
            file.write(">Number of contigs: " + str(len(contigs)) + "\n")
            for i, contig in enumerate(contigs):
                file.write(">CONTIG" + str(i+1) + "\n")
                file.write(contig)
                file.write("\n")
            file.write(">Time finished: " + time.strftime("%c", time.localtime()))


class DebugIOHandler(IOHandler):
    @staticmethod
    def read_input(print_runtime=True, print_syssizeof=False, start_time=0):
        if print_runtime:
            print("\n>--- STARTING ASSEMBLY PROGRAM AT T = 0.00 ---")

        reads, paired, d, num_bases = IOHandler.read_input()

        if print_runtime:
            print(">FINISHED READING INPUT AT T = {:.2f}".format(time.time()-start_time))
        if print_syssizeof:
            total_string_memory = 0
            for read in reads:
                total_string_memory += sys.getsizeof(read)
            print(">SIZE OF READ CONTAINER: {:,}".format(sys.getsizeof(reads)))
            print(">SIZE OF ALL READ STRINGS: {:,}".format(total_string_memory))

        return (reads, paired, int(d), num_bases)

    @staticmethod
    def write_FASTQ(contigs, constants, start_time, timeit):
        if timeit:
            print(">--- WRITING TO FILE AT T = {:.2f} ---".format(time.time()-start_time))
        IOHandler.write_FASTQ(contigs, constants, start_time)


def assemble_with_options(print_runtime=False, print_syssizeof=False, using_count_min_sketch=False,
                          start_time=time.time(), k=None, hamming_dist=None, paired_error=None):

    reads, paired, _, _ = DebugIOHandler.read_input(
        print_runtime=print_runtime, start_time=start_time, print_syssizeof=print_syssizeof)

    graph = None
    if using_count_min_sketch:
        if paired:
            graph = DebugCMSPairedDeBruijnGraph(
                print_runtime=print_runtime, start_time=start_time, print_syssizeof=print_syssizeof,
                reads=reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)
        else:
            graph = DebugCMSDeBruijnGraph(
                print_runtime=print_runtime, start_time=start_time, print_syssizeof=print_syssizeof,
                reads=reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)
    else:
        if paired:
            graph = DebugPairedDeBruijnGraph(
                print_runtime=print_runtime, start_time=start_time, print_syssizeof=print_syssizeof,
                reads=reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)
        else:
            graph = DebugDeBruijnGraph(
                print_runtime=print_runtime, start_time=start_time, print_syssizeof=print_syssizeof,
                reads=reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)

    contigs = graph.enumerate_contigs()
    constants = (graph.KMER_LEN, graph.HAMMING_DIST, graph.ALLOWED_PAIRED_DIST_ERROR)
    return contigs, constants


def assemble_without_options(k=None, hamming_dist=None, paired_error=None):
    reads, paired, _, _ = IOHandler.read_input()
    graph = None
    if paired:
        graph = PairedDeBruijnGraph(reads, k=k, hamming_dist=hamming_dist,
                                    paired_error=paired_error)
    else:
        graph = DeBruijnGraph(reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)
    contigs = graph.enumerate_contigs()
    constants = (graph.KMER_LEN, graph.HAMMING_DIST, graph.ALLOWED_PAIRED_DIST_ERROR)
    return contigs, constants


def main():
    start_time = time.time()
    args = IOHandler.read_args()

    contigs = list()
    constants = tuple()

    if args.time or args.memory or args.count_min_sketch:
        contigs, constants = assemble_with_options(
            print_runtime=args.time, print_syssizeof=args.memory,
            using_count_min_sketch=args.count_min_sketch, start_time=start_time,
            k=args.kmer_length, hamming_dist=args.filter_threshold, paired_error=args.error)
    else:
        contigs, constants = assemble_without_options(
            k=args.kmer_length, hamming_dist=args.filter_threshold, paired_error=args.error)

    if args.stdout:
        IOHandler.write_stdout(contigs, constants, start_time)
    else:
        DebugIOHandler.write_FASTQ(contigs, constants, start_time, timeit=args.time)

    if args.time:
        print(">--- PROGRAM FINISHED AT T = {:.2f} ---".format(time.time()-start_time))


if __name__ == "__main__":
    main()
    # import cProfile
    # cProfile.run('main()')
