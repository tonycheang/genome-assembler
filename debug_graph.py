import time
import sys
from collections import deque, defaultdict
from countminsketch import CountMinSketch
from debruijn_graph import DeBruijnGraph, PairedDeBruijnGraph
from debruijn_graph import CMSDeBruijnGraph, CMSPairedDeBruijnGraph


class DebugGraph:
    ''' Class to compose with various DeBruijnGraph classes for debugging,
        via multiple inheritance. Not meant to instantiate on its own.
    '''

    def __init__(self, print_syssizeof=False, print_runtime=False, start_time=0, **kwargs):
        self.print_syssizeof = print_syssizeof
        self.print_runtime = print_runtime
        self.start_time = start_time
        super().__init__(**kwargs)

    def enumerate_contigs(self):
        if self.print_runtime:
            print(
                "\n>--- STARTING TO ENUMERATE CONTIGS AT T = {:.2f} ---".format(time.time()-self.start_time))
        contigs = super().enumerate_contigs()
        if self.print_runtime:
            print(
                ">FINISHED ENUMERATING CONTIGS AT T = {:.2f} ---".format(time.time()-self.start_time))
        return contigs

    def _build_graph(self, kmer_counts, reads):
        if self.print_runtime:
            print(
                "\n>--- STARTING TO BUILD GRAPH AT T = {:.2f} ---".format(time.time()-self.start_time))
        super()._build_graph(kmer_counts, reads)
        if self.print_runtime:
            print(">FINISHED BUILDING GRAPH AT T = {:.2f}".format(time.time()-self.start_time))
        if self.print_syssizeof:
            node_mem = 0
            container_mem = sys.getsizeof(self.nodes)
            for key, node in self.nodes.items():
                container_mem += sys.getsizeof(key)
                node_mem += sys.getsizeof(node)
            print(">SIZE OF GRAPH CONTAINER: {:,}".format(container_mem))
            print(">SIZE OF ALL NODES: {:,}".format(node_mem))

    def _count_kmers(self, k, reads):
        if self.print_runtime:
            print(
                "\n>--- STARTING TO COUNT KMERS AT T = {:.2f} ---".format(time.time()-self.start_time))

        kmer_counts_dict = super()._count_kmers(k, reads)

        if self.print_runtime:
            print(">FINISHED COUNTING KMERS AT T = {:.2f}".format(time.time()-self.start_time))
        if self.print_syssizeof:
            string_mem = 0
            container_mem = sys.getsizeof(kmer_counts_dict)
            if kmer_counts_dict:
                container_mem += sum(map(lambda x: sys.getsizeof(x[0]), kmer_counts_dict.items()))
                string_mem += sum(map(lambda x: sys.getsizeof(x[1]), kmer_counts_dict.items()))

            print(">SIZE OF COUNTS CONTAINER: {:,}".format(container_mem))
            print(">SIZE OF STRINGS IN COUNTS: {:,}".format(string_mem))
        return kmer_counts_dict

    def _make_sketch(self, kmer_counts_dict: defaultdict) -> CountMinSketch:
        if self.print_runtime:
            print(
                "\n>--- STARTING TO MAKE COUNTMIN SKETCH AT T = {:.2f} ---".format(time.time()-self.start_time))

        # Read the dictionary into a compressed data structure
        NUM_ROWS = 10
        kmer_counts = CountMinSketch(NUM_ROWS)
        for i, (kmer, count) in enumerate(kmer_counts_dict.items()):
            if self.print_runtime and i % 50000 == 0:
                print(">Processed {0} kmers by time T={1:.2f}".format(
                    i, time.time()-self.start_time))
            kmer_counts.update(kmer, count)

        if self.print_runtime:
            print(">FINISHED MAKING COUNTMIN SKETCH AT T = {:.2f}".format(
                time.time()-self.start_time))
        if self.print_syssizeof:
            print(">SIZE OF COUNTMIN SKETCH: {:,}".format(sys.getsizeof(kmer_counts)))
        return kmer_counts


class DebugDeBruijnGraph(DebugGraph, DeBruijnGraph):
    """ Delegates its instantiation to DebugGraph. Composes DebugGraph with DeBruijnGraph. """
    pass


class DebugCMSDeBruijnGraph(DebugGraph, CMSDeBruijnGraph):
    """ Delegates its instantiation to DebugGraph. Composes DebugGraph with CMSDeBruijnGraph. """
    pass


class DebugPairedDeBruijnGraph(DebugGraph, PairedDeBruijnGraph):
    """ Delegates its instantiation to DebugGraph. Composes DebugGraph with PairedDeBruijnGraph. """
    pass


class DebugCMSPairedDeBruijnGraph(DebugGraph, CMSPairedDeBruijnGraph):
    """ Delegates its instantiation to DebugGraph. Composes DebugGraph with CMSPairedDeBruijnGraph. """
    pass
