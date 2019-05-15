# python3

import sys
from abc import ABC, abstractmethod
from collections import deque, defaultdict
from itertools import tee, combinations
from functools import lru_cache
import time

import cProfile
import tracemalloc

""" Implementing an Assembler

    --- Goal ---

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. Instead of outputting an Eulerian path--unique
    construction can be rare--this assembler outputs contigs (non-branching segments).

    Data Sets Tested Against:
    N. Deltocephalinicola                - t = 34000,   len = 100,    coverage: 30
    N. Deltocephalinicola (read-pairs)   - t = 19679,   len = varies, coverage: unknown (likely around 30)
    E. Coli O104                         - t = 1400000, len = 100,    coverage: 25
    E. Coli O104          (read-pairs)   - t = 700000,  len = 100,    coverage: 25

    --- Approach ---

        - Building the Graph -

    (1) Filter kmer or kmer-pairs with low coverage. (Hamming Distance approach)
            Can use Count-min sketch to store counts using log(n) space if space becomes a concern.
            We can ignore bubbles and tips after this.
    (2) Build two nodes for each edge to create an unconnected graph.
    (3) Merge nodes (a|b) and (a|b') when b and b' overlap by more than length 2*DELTA, 
            where DELTA is the misalign distance. This creates a connected simple graph.
    (4) Each node at this point must record whether they are a was_branching node or not.
            Otherwise, once the last edges remain from popping, it will be counted as part of another contig.
    
        - Building Contigs -
    
    (5) Traverse each edge forward to obtain a non-was_branching segment. This is a contig.
            If starting at a branch (i.e. in/outdegree > 1) or at the beginning of a tip,
            take any node forward and only go forward until next branch.
                Keep doing this for any remaining edges out of the original node.

            Traversal backward needs the correct forward edge removed. Maybe store in a dict?
    
        - Infering Multiplicity of Contigs -
        (SKIP for now. Infer only if necessary. Contigs of repeated regions will be unrepresented.)

    (5.5) On was_branching nodes, record where each contig joins another.
            The was_branching nodes with contigs as edges will make up the circulation flow graph.
    (6) Build a graph where each edge is a contig that connects at the correct location to other contigs.
            Long contigs (>1000 BP) will have max and min flow of 1. They are correct.
    (7) Run a max circulation flow through the network to infer multiplicity. 
            O(V[E**2]) doable for ~300 contigs, even 10**3 contigs.
            O(VE) algorithm exists if runtime becomes a concern.

    Do not store all k-1mers when counting them. Re-slice when building graph.
        This should save a lot of string overhead and runtime is negligible.
    CountMin Sketch to estimate counts of each kmer to store in O(log(n)) space
        Profile to make sure the counts array is actually large...
"""


class Node:
    def __init__(self, data):
        # Data contains the k-1mer/prefix/suffix
        self.data = data
        # Edges are the string, or the key into the constructor
        self.edges = dict()
        # Bool may be set to True at the end of building the graph
        self.was_branching = False
        self.num_edges_in = 0

    @property
    def outdegree(self):
        ''' Outdegree changes as edges get popped '''
        return len(self.edges)

    @property
    def indegree(self):
        ''' Indegree stays set after graph is built '''
        return self.num_edges_in

    def pop_edge(self):
        ''' Gets and removes an edge. '''
        return self.edges.popitem()

    def append_edge(self, edge):
        self.edges[edge] = True

    def __repr__(self):
        ''' For debugging small examples. '''
        return "Data: {0} | Edges: {2}".format(self.data, self.edges)
    
    def __sizeof__(self):
        if hasattr(self, "total_mem"):
            return self.total_mem
        total_mem = super().__sizeof__()
        total_mem += sys.getsizeof(self.data)
        total_mem += sys.getsizeof(self.edges)
        for _, edge in self.edges:
            total_mem += sys.getsizeof(edge)
        total_mem += sys.getsizeof(self.was_branching)
        total_mem += sys.getsizeof(self.num_edges_in)
        self.total_mem = total_mem
        self.total_mem += sys.getsizeof(self.total_mem)
        return total_mem


class PairedNode(Node):
    def __init__(self, data, paired_data):
        self.paired_data = paired_data
        super().__init__(data)
    
    def __repr__(self):
        ''' For debugging small examples. '''
        return "Data: {0} | Pair: {1} | Edges: {2}".format(self.data, self.paired_data, self.edges)
    
    def __sizeof__(self):
        if hasattr(self, "total_mem"):
            return self.total_mem
        total_mem = super().__sizeof__() + sys.getsizeof(self.paired_data)
        self.total_mem = total_mem
        return total_mem


class AbstractDeBruijnGraph(ABC):
    @abstractmethod
    def enumerate_contigs(self):
        pass

    @abstractmethod
    def _get_longest_contig(self, start_node_data, visited):
        pass

    @abstractmethod
    def _build_graph(self, reads):
        pass

    @staticmethod
    @abstractmethod
    def _count_kmers(reads):
        pass

    @staticmethod
    @abstractmethod
    def _break_read_into_k_minus_one_mers(k, read):
        ''' Breaks into list of prefix_pairs and suffix_pairs. '''
        pass

    @staticmethod
    def _pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


class DeBruijnGraph(AbstractDeBruijnGraph):
    ''' Implementation relies on k-1mer being a single key into self.nodes.
        This is unlike the PairedDeBruijnGraph, and allows for kmer_counts
        to hold the entire edge, removing the need to keep track of reads
        or all the k-1mers in a separate structure.
    '''

    KMER_LEN = 29
    HAMMING_DIST = 3
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D,
    # the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads):
        self.num_edges = 0
        # Indexed by data/prefix_paired/suffix_paired.
        self.nodes = dict()
        kmer_counts = DeBruijnGraph._count_kmers(reads)
        del reads
        self._build_graph(kmer_counts)

    def enumerate_contigs(self) -> list:
        contigs = list()
        for _, node in self.nodes.items():
            if node.outdegree > 0:
                while node.outdegree > 0 and (node.was_branching or node.indegree == 0):          #
                    contig = self._get_longest_contig(node)
                    contigs.append(contig)
                if self.num_edges == 0:
                    return contigs
        return contigs

    def _get_longest_contig(self, cur_node) -> str:
        ''' Finds the longest contig by moving both forward and backward until
            nodes with branches are found.
        '''
        contig_pieces = deque()
        # Ensures movement forward at least one edge.
        edge, _ = cur_node.pop_edge()
        self.num_edges -= 1
        contig_pieces.append(edge[-1])
        cur_node = self.nodes[edge]

        while cur_node.outdegree > 0 and not cur_node.was_branching:
            edge, _ = cur_node.pop_edge()
            self.num_edges -= 1
            # Traversal forward adds the last character of next k-1mer
            contig_pieces.append(edge[-1])
            cur_node = self.nodes[edge]
        
        return "".join(contig_pieces)

    def _build_graph(self, kmer_counts):
        ''' Builds the path constructor with coverage considerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        for edge, count in kmer_counts.items():
            prefix, suffix = edge[0], edge[1]
            if count > DeBruijnGraph.HAMMING_DIST:
                if prefix not in self.nodes:
                    self.nodes[prefix] = Node(prefix)
                if suffix not in self.nodes:
                    self.nodes[suffix] = Node(suffix)

                self.nodes[prefix].append_edge(suffix)
                self.nodes[suffix].num_edges_in += 1
                self.num_edges += 1
        
        # Mark was_branching nodes before graph gets consumed.
        for _, node in self.nodes.items():
            if node.outdegree > 1 or node.indegree > 1:
                node.was_branching = True

    @staticmethod
    def _count_kmers(reads):
        kmer_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = DeBruijnGraph._break_read_into_k_minus_one_mers(
                DeBruijnGraph.KMER_LEN, cur_read)

            for prefix, suffix in AbstractDeBruijnGraph._pairwise(k_minus_one_mers):
                kmer_counts[(prefix, suffix)] += 1
        return kmer_counts

    @staticmethod
    def _break_read_into_k_minus_one_mers(k, read, paired=False):
        ''' Non-paired output format for 'ACTGAC', k=4: ['ACT', 'CTG', 'TGA', 'GAC'] '''
        return [sys.intern(read[i:i+k-1]) for i in range(len(read)-(k-2))]


class PairedDeBruijnGraph(AbstractDeBruijnGraph):
    ''' Implementation relies on idea of dual-key with paired-reads.
        That is, self.nodes is a defaultdict of dicts; two keys are necessary to identify a node.
        A tuple would have sufficed if paired reads were drawn with a perfect distance,
        but since it's inexact, it's helpful to be able to iterate through the secondary dict
        for key comparison.
    '''

    KMER_LEN = 27
    HAMMING_DIST = 5
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D, the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads, d):
        self.num_edges = 0
        # Indexed by (data, paired_data)/(prefix, paired_prefix)/(suffix, paired_suffix)
        self.nodes = defaultdict(dict)
        self.dist = d
        kmer_counts, broken_read_pairs = self._count_kmers(reads)
        del reads
        self._build_graph(kmer_counts, broken_read_pairs)

    def enumerate_contigs(self) -> list:
        # What if graph is perfectly circular?
        contigs = list()
        for _, nodes in self.nodes.items():
            for _, node in nodes.items():
                while node.outdegree > 0 and (node.was_branching or node.indegree == 0):          #
                    contig = self._get_longest_contig(node)
                    contigs.append(contig)
                if self.num_edges == 0:
                    return contigs
        return contigs

    def _get_longest_contig(self, cur_node) -> str:
        ''' Builds up the contig in the specified direction until either a branching
            node is encountered or there are no more edges to traverse.
        '''

        contig_pieces = deque()
        # Ensures movement forward at least one edge.
        edge, _ = cur_node.pop_edge()
        self.num_edges -= 1
        key, paired_key = edge
        contig_pieces.append(key[-1])
        cur_node = self.nodes[key][paired_key]

        while cur_node.outdegree > 0 and not cur_node.was_branching:
            edge, _ = cur_node.pop_edge()
            self.num_edges -= 1
            key, paired_key = edge
            # Traversal forward adds the last character of next k-1mer
            contig_pieces.append(key[-1])
            cur_node = self.nodes[key][paired_key]

        return "".join(contig_pieces)

    def _build_graph(self, kmer_counts, broken_read_pairs) -> tuple:

        for read_pairs in broken_read_pairs:
            for prefix_paired, suffix_paired in AbstractDeBruijnGraph._pairwise(read_pairs):
                # Indiscriminately filters any edge that has an erroneous k-mer
                if kmer_counts[(prefix_paired[0], suffix_paired[0])] > PairedDeBruijnGraph.HAMMING_DIST and \
                        kmer_counts[(prefix_paired[1], suffix_paired[1])] > PairedDeBruijnGraph.HAMMING_DIST:

                    # Check whether the each node already exists, accounting for inexact distance between paired reads.
                    # self._find_matching_node() will set the node variables to existing  nodes if they exist.
                    prefix_found, prefix_node = self._find_matching_node(prefix_paired)
                    suffix_found, suffix_node = self._find_matching_node(suffix_paired)

                    # Create and add any nodes that don't already exist.
                    # The check for not found ensures we do not overwrite existing matching nodes with new nodes.
                    if not prefix_found:
                        prefix_node = PairedNode(
                            data=sys.intern(prefix_paired[0]), paired_data=sys.intern(prefix_paired[1]))
                        self.nodes[prefix_paired[0]][prefix_paired[1]] = prefix_node
                    if not suffix_found:
                        suffix_node = PairedNode(
                            data=sys.intern(suffix_paired[0]), paired_data=sys.intern(suffix_paired[1]))
                        self.nodes[suffix_paired[0]][suffix_paired[1]] = suffix_node

                    # Case: Either not node found, therefore edge cannot already exist.
                    if not (prefix_found and suffix_found):
                        prefix_node.append_edge((suffix_node.data, suffix_node.paired_data))
                        suffix_node.num_edges_in += 1
                        self.num_edges += 1
                    # Case: Both nodes found but edge between them doesn't exist.
                    elif (suffix_node.data, suffix_node.paired_data) not in prefix_node.edges:
                        prefix_node.append_edge((suffix_node.data, suffix_node.paired_data))
                        suffix_node.num_edges_in += 1
                        self.num_edges += 1
                    # Case: Both nodes found, and an edge between them exists. Don't add a new edge.
                    # Ensures simple graph created. Self-edges not accounted for.
                    else:
                        continue

        # Mark was_branching nodes before graph gets consumed.
        for _, nodes in self.nodes.items():
            for _, node in nodes.items():
                if node.outdegree > 1 or node.indegree > 1:
                    node.was_branching = True

    def _find_matching_node(self, paired_strings) -> tuple:
        ''' Checks whether self.nodes has a node (A|B) node with matching A
            and a matching B. Matching A must be exact, while matching B can
            misalign by 2*DELTA or ALLOWED_PAIRED_DIST_ERROR
        '''
        if paired_strings[0] in self.nodes:
            # Early exit if possible.
            if paired_strings[1] in self.nodes[paired_strings[0]]:
                return (True, self.nodes[paired_strings[0]][paired_strings[1]])

            for paired_key, potential_match_node in self.nodes[paired_strings[0]].items():
                # Check both alignments, since paired-reads are drawn from uniform distribution AROUND d.
                if PairedDeBruijnGraph._find_longest_overlap_brute(paired_key, paired_strings[1]) != 0 or \
                        PairedDeBruijnGraph._find_longest_overlap_brute(paired_strings[1], paired_key) != 0:
                    return (True, potential_match_node)
        return (False, None)

    @staticmethod
    @lru_cache(maxsize=512)
    def _find_longest_overlap_brute(pattern, text) -> int:
        for start_text_pos in range(len(text) - PairedDeBruijnGraph.ALLOWED_PAIRED_DIST_ERROR):
            len_possible_overlap = min(len(text) - start_text_pos, len(pattern))
            for pattern_pos in range(len_possible_overlap):
                text_pos = start_text_pos + pattern_pos
                if text[text_pos] != pattern[pattern_pos]:
                    break
            else:
                return len_possible_overlap
        return 0

    @staticmethod
    def _count_kmers(reads) -> tuple:
        ''' Generates a structure to check whether a kmer appears frequently enough.
            Currently uses a dict but can switch over a count-min sketch if space is a concern.
        '''

        kmer_counts = defaultdict(int)
        broken_read_pairs = list()
        for cur_read_pair in reads:
            # Need to count all kmers for error correction / filtering.
            paired_k_minus_one_mers = PairedDeBruijnGraph._break_read_into_k_minus_one_mers(
                PairedDeBruijnGraph.KMER_LEN, cur_read_pair)

            broken_read_pairs.append(paired_k_minus_one_mers)

            for prefix_paired_paired, suffix_paired_paired in AbstractDeBruijnGraph._pairwise(paired_k_minus_one_mers):
                # Each kmer_pair generates two kmers, each of which we'll count for coverage.
                kmer_counts[(prefix_paired_paired[0], suffix_paired_paired[0])] += 1
                kmer_counts[(prefix_paired_paired[1], suffix_paired_paired[1])] += 1
        return (kmer_counts, broken_read_pairs)

    @staticmethod
    def _break_read_into_k_minus_one_mers(k, read) -> list:
        ''' Paired output format for ('ACTGAC', 'TCGATC'), k=4: 
            [('ACT', 'TCG'), ('CTG', 'CGA'), ('TGA', 'GAT'), ('GAC', 'ATC')]
        '''
        return [(sys.intern(read[0][i:i+k-1]), sys.intern(read[1][i:i+k-1])) for i in range(len(read[0])-(k-2))]


class DebugPairedDeBruijnGraph(PairedDeBruijnGraph):
    def __init__(self, reads, d, print_snapshot=False, print_syssizeof=False, print_runtime=False, start_time=0):
        self.print_syssizeof = print_syssizeof
        self.print_snapshot = print_snapshot
        self.print_runtime = print_runtime
        self.start_time = start_time
        super().__init__(reads, d)

    def enumerate_contigs(self):
        if self.print_runtime:
            print("\n--- STARTING TO ENUMERATE CONTIGS AT T = {:.2f} ---".format(time.time()-self.start_time))
        contigs = super().enumerate_contigs()
        if self.print_runtime:
            print("FINISHED ENUMERATING CONTIGS AT T = {:.2f} ---".format(time.time()-self.start_time))
        if self.print_snapshot:
            print_memory_snapshot("AFTER ENUMERATE CONTIGS")
        return contigs

    def _build_graph(self, kmer_counts, broken_read_pairs):
        if self.print_runtime:
            print("\n--- STARTING TO BUILD GRAPH AT T = {:.2f} ---".format(time.time()-self.start_time))
        super()._build_graph(kmer_counts, broken_read_pairs)
        if self.print_runtime:
            print("FINISHED BUILDING GRAPH AT T = {:.2f}".format(time.time()-self.start_time))
        if self.print_syssizeof:
            node_mem = 0
            container_mem = sys.getsizeof(self.nodes)
            for _, node_dict in self.nodes.items():
                container_mem += sys.getsizeof(node_dict)
                for _, node in node_dict.items():
                    node_mem += sys.getsizeof(node)
            print("SIZE OF GRAPH CONTAINER: {:,}".format(container_mem))
            print("SIZE OF ALL NODES: {:,}".format(node_mem))
        if self.print_snapshot:
            print_memory_snapshot("AFTER BUILDING GRAPH")
            

    def _count_kmers(self, reads):
        if self.print_runtime:
            print("\n--- STARTING TO COUNT KMERS AT T = {:.2f} ---".format(time.time()-self.start_time))
        kmer_counts, broken_read_pairs = super()._count_kmers(reads)
        if self.print_runtime:
            print("FINISHED COUNTING KMERS AT T = {:.2f}".format(time.time()-self.start_time))
        if self.print_syssizeof:
            string_mem = 0
            container_mem = sys.getsizeof(kmer_counts)
            print("SIZE OF COUNTS CONTAINER WITHOUT TUPLES: {:,}".format(container_mem))
            for pair in kmer_counts.items():
                container_mem += sys.getsizeof(pair)
                string_mem += sys.getsizeof(pair[0])
                string_mem += sys.getsizeof(pair[1])
            print("SIZE OF COUNTS CONTAINER: {:,}".format(container_mem))
            print("SIZE OF STRINGS IN COUNTS: {:,}".format(string_mem))
            container_mem = sys.getsizeof(broken_read_pairs)
            string_mem = 0
            for read_pairs in broken_read_pairs:
                container_mem += sys.getsizeof(read_pairs)
                for pair in read_pairs:
                    container_mem += sys.getsizeof(pair)
                    string_mem += sys.getsizeof(pair[0])
                    string_mem += sys.getsizeof(pair[1])
            print("SIZE OF BROKEN-READ-PAIRS CONTAINER: {:,}".format(container_mem))
            print("SIZE OF STRINGS IN BROKEN-READ-PAIRS: {:,}".format( string_mem))
        if self.print_snapshot:
            print_memory_snapshot("AFTER COUNTING KMERS")
        return kmer_counts, broken_read_pairs


class IOHandler:
    @staticmethod
    def read_input():
        ''' Returns a list of reads or read-pairs, not yet broken down,
            and indications of which (paired?, d).        
        '''
        n = int(sys.stdin.readline().strip())
        paired = False
        reads = list()
        d = 0

        # Check for read-pair vs non-paired input format
        first_line = sys.stdin.readline().strip().split('|')
        if len(first_line) > 1:
            read1, read2, d = first_line
            reads.append((read1, read2))
            paired = True
            for _ in range(n-1):
                read1, read2, d = sys.stdin.readline().strip().split('|')
                reads.append((read1, read2))
        else:
            reads.append(first_line[0])
            for _ in range(n-1):
                read = sys.stdin.readline().strip()
                reads.append(read)

        return (reads, paired, int(d))

    @staticmethod
    def print_FASTA(contigs):
        for i, contig in enumerate(contigs):
            print(">CONTIG" + str(i+1))
            print(contig)


class DebugIOHandler(IOHandler):
    @staticmethod
    def read_input(print_runtime=True, print_snapshot=False, print_syssizeof=False, start_time=0):
        if print_runtime:
            print("\n--- STARTING ASSEMBLY PROGRAM AT T = 0.00 ---")

        reads, paired, d = IOHandler.read_input()

        if print_runtime:
            print("FINISHED READING INPUT AT T = {:.2f}".format(time.time()-start_time))
        if print_syssizeof:
            total_string_memory = 0
            for read in reads:
                total_string_memory += sys.getsizeof(read)
            print("SIZE OF READ CONTAINER: {:,}".format(sys.getsizeof(reads)))
            print("SIZE OF ALL READ STRINGS: {:,}".format(total_string_memory))
        if print_snapshot:
            print_memory_snapshot("AFTER READING INPUT")

        return (reads, paired, int(d))


def print_memory_snapshot(location):
    print("\n")
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')
    # top_stats = snapshot.statistics('traceback')

    print("[Top 5] MEMORY USE for " + location)
    for stat in top_stats[:5]:
        print(stat)
        # for line in stat.traceback.format():
        #     print(line)
    return


def profile_assembler(print_runtime=False, print_snapshot=False, print_syssizeof=False):

    if print_snapshot:
        tracemalloc.start(25)
    start_time = time.time()

    reads, paired, d = DebugIOHandler.read_input(
        print_runtime=print_runtime, print_snapshot=print_snapshot, start_time=start_time,
        print_syssizeof=print_syssizeof)
    graph = None
    if paired:
        graph = DebugPairedDeBruijnGraph(reads=reads, d=d, print_runtime=print_runtime,
                                         start_time=start_time, print_syssizeof=print_syssizeof,
                                         print_snapshot=print_snapshot)
    else:
        graph = DeBruijnGraph(reads)
    contigs = graph.enumerate_contigs()
    DebugIOHandler.print_FASTA(contigs)

    if print_runtime:
        print("--- PROGRAM FINISHED AT T = {:.2f} ---".format(time.time()-start_time))


def main():
    reads, paired, d = IOHandler.read_input()
    graph = None
    if paired:
        graph = PairedDeBruijnGraph(reads, d)
    else:
        graph = DeBruijnGraph(reads)
    contigs = graph.enumerate_contigs()
    IOHandler.print_FASTA(contigs)


if __name__ == "__main__":
    sys.setswitchinterval(10000)
    # main()
    profile_assembler(print_runtime=True, print_syssizeof=True, print_snapshot=False)
    # cProfile.run('main()')
    #cProfile.run('profile_assembler(print_runtime=True, print_syssizeof=False, print_snapshot=False)')
