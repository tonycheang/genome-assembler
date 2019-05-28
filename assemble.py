# python3

from abc import ABC, abstractmethod
from collections import deque, defaultdict
from itertools import tee, combinations
from functools import lru_cache

from random import randint

from array import array
from math import log

import sys
import argparse
import time
import cProfile

''' Implementing an Assembler

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

    Maybe instead of checking for kmer alignment/repeats, we check for read alignment/repeats
    Then we... don't we need to count all the kmers anyway?
    Or maybe we can store the reads instead in the CountMinSketch or whatever we use to check fidelity.
    
'''

""" ----- STORAGE STRUCTURES ----- """


class CountMinSketch:
    # Various arrays of primes for testing purposes.
    # 20 primes cenetered around 6*10**7
    primes_6_10_7 = [59999879, 59999881, 59999887, 59999917, 59999957,
                     59999971, 59999981, 59999983, 59999993, 59999999,
                     60000011, 60000013, 60000023, 60000047, 60000049,
                     60000067, 60000071, 60000091, 60000103, 60000113]

    # 20 primes cenetered around 10**7
    primes_1_10_7 = [9999889, 9999901, 9999907, 9999929, 9999931,
                     9999937, 9999943, 9999971, 9999973, 9999991,
                     10000019, 10000079, 10000103, 10000121, 10000139,
                     10000141, 10000169, 10000189, 10000223, 10000229]

    # 20 primes cenetered around 5*10**6
    primes_5_10_6 = [4999871, 4999879, 4999889, 4999913, 4999933,
                     4999949, 4999957, 4999961, 4999963, 4999999,
                     5000011, 5000077, 5000081, 5000087, 5000101,
                     5000111, 5000113, 5000153, 5000161, 5000167]

    def __init__(self, num_rows):
        assert num_rows < len(
            CountMinSketch.primes), "Requested number of rows in CountMinSketch exceeds number of primes in list"
        self.num_rows = num_rows
        # Store in arrays with short ints to save space
        self.hash_values = [array('H', [0] * CountMinSketch.primes_1_10_7[i])
                            for i in range(num_rows)]

    def update(self, string, amount):
        string_hash = self._hash(string)
        for row in range(self.num_rows):
            self.hash_values[row][string_hash %
                                  len(self.hash_values[row])] += amount  # * self._get_sign(string_hash, row)

    def estimate(self, string):
        string_hash = self._hash(string)
        est_from_min = min(self.hash_values[row][string_hash %
                                                 len(self.hash_values[row])] for row in range(self.num_rows))
        return est_from_min

    @staticmethod
    @lru_cache(maxsize=512)
    def _hash(data, seed=0):
        ''' MurmurHash3 x86_32 port from Java by Maurus Decimus via StackOverflow'''
        c1 = 0xcc9e2d51
        c2 = 0x1b873593

        length = len(data)
        h1 = seed
        roundedEnd = (length & 0xfffffffc)  # round down to 4 byte block
        for i in range(0, roundedEnd, 4):
            # little endian load order
            k1 = (ord(data[i]) & 0xff) | ((ord(data[i + 1]) & 0xff) << 8) | \
                ((ord(data[i + 2]) & 0xff) << 16) | (ord(data[i + 3]) << 24)
            k1 *= c1
            k1 = (k1 << 15) | ((k1 & 0xffffffff) >> 17)  # ROTL32(k1,15)
            k1 *= c2

            h1 ^= k1
            h1 = (h1 << 13) | ((h1 & 0xffffffff) >> 19)  # ROTL32(h1,13)
            h1 = h1 * 5 + 0xe6546b64

        # tail
        k1 = 0

        val = length & 0x03
        if val == 3:
            k1 = (ord(data[roundedEnd + 2]) & 0xff) << 16
        # fallthrough
        if val in [2, 3]:
            k1 |= (ord(data[roundedEnd + 1]) & 0xff) << 8
        # fallthrough
        if val in [1, 2, 3]:
            k1 |= ord(data[roundedEnd]) & 0xff
            k1 *= c1
            k1 = (k1 << 15) | ((k1 & 0xffffffff) >> 17)  # ROTL32(k1,15)
            k1 *= c2
            h1 ^= k1

        # finalization
        h1 ^= length

        # fmix(h1)
        h1 ^= ((h1 & 0xffffffff) >> 16)
        h1 *= 0x85ebca6b
        h1 ^= ((h1 & 0xffffffff) >> 13)
        h1 *= 0xc2b2ae35
        h1 ^= ((h1 & 0xffffffff) >> 16)

        return h1 & 0xffffffff

    def __getitem__(self, key):
        # To allow working similar to dict
        return self.estimate(key)

    def __sizeof__(self):
        if hasattr(self, "total_mem"):
            return self.total_mem
        total_mem = 0
        total_mem += sys.getsizeof(self.num_rows)
        total_mem += sys.getsizeof(self.hash_values)
        for row in self.hash_values:
            total_mem += row.buffer_info()[1] * row.itemsize
        self.total_mem = total_mem
        self.total_mem += sys.getsizeof(self.total_mem)
        return total_mem


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
        # print(self.edges)
        if self.edges:
            if isinstance(next(iter(self.edges)), tuple):
                for _, edge in self.edges:
                    total_mem += sys.getsizeof(edge)
            elif isinstance(next(iter(self.edges)), str):
                for edge in self.edges:
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


""" ----- REGULAR GRAPHS ----- """


class AbstractDeBruijnGraph(ABC):
    @abstractmethod
    def enumerate_contigs(self):
        pass

    @abstractmethod
    def _get_longest_contig(self, start_node):
        pass

    @abstractmethod
    def _build_graph(self, reads:list):
        pass

    @staticmethod
    @abstractmethod
    def _count_kmers(k:int, reads:list):
        pass

    @staticmethod
    @abstractmethod
    def _break_read_into_k_minus_one_mers(k:int, read:str):
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
    ALLOWED_PAIRED_DIST_ERROR = 2

    def __init__(self, reads: list, k=None, hamming_dist=None, paired_error=None):
        # Set constants if provided, otherwise defaults to class values
        if k is not None:
            self.KMER_LEN = k
        if hamming_dist is not None:
            self.HAMMING_DIST = hamming_dist
        if paired_error is not None:
            self.ALLOWED_PAIRED_DIST_ERROR = paired_error

        self.num_edges = 0
        # Indexed by data/prefix_paired/suffix_paired.
        self.nodes = dict()
        kmer_counts_dict = self._count_kmers(self.KMER_LEN, reads)
        self._build_graph(kmer_counts_dict, reads)

    def enumerate_contigs(self) -> list:
        contigs = list()
        total_len = 0
        for _, node in self.nodes.items():
            while node.outdegree > 0 and (node.was_branching or node.indegree == 0):          #
                contig = self._get_longest_contig(node)
                contigs.append(contig)
                total_len += len(contig)
            if self.num_edges == 0:
                return contigs

        # Handles when some sub contig is perfectly circular
        # All non-circular segments will have been found by iteration above.
        for _, nodes in self.nodes.items():
            while node.outdegree > 0:          #
                contig = self._get_longest_contig(node)
                contigs.append(contig)
                total_len += len(contig)
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

    def _build_graph(self, kmer_counts, reads: list):
        ''' Builds the path constructor with coverage considerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        for i, read in enumerate(reads):
            # if i % 50000 == 0:
            #     print("Processing read no.", i)
            k_minus_one_mers = self._break_read_into_k_minus_one_mers(self.KMER_LEN, read)
            for prefix, suffix in AbstractDeBruijnGraph._pairwise(k_minus_one_mers):
                # Check for novel edges
                if prefix in self.nodes and suffix in self.nodes:
                    if suffix in self.nodes[prefix].edges:
                        continue

                if kmer_counts[prefix] > self.HAMMING_DIST and \
                        kmer_counts[suffix] > self.HAMMING_DIST:

                    if prefix not in self.nodes:
                        self.nodes[prefix] = Node(sys.intern(prefix))
                    if suffix not in self.nodes:
                        self.nodes[suffix] = Node(sys.intern(suffix))

                    self.nodes[prefix].append_edge(sys.intern(suffix))
                    self.nodes[suffix].num_edges_in += 1
                    self.num_edges += 1

        # Mark was_branching nodes before graph gets consumed.
        for _, node in self.nodes.items():
            if node.outdegree > 1 or node.indegree > 1:
                node.was_branching = True

    @staticmethod
    def _count_kmers(k: int, reads: list) -> defaultdict:
        kmer_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = DeBruijnGraph._break_read_into_k_minus_one_mers(k, cur_read)

            for k_minus_one_mer in k_minus_one_mers:
                kmer_counts[k_minus_one_mer] += 1
        return kmer_counts

    @staticmethod
    def _break_read_into_k_minus_one_mers(k: int, read: str, paired=False) -> list:
        ''' Non-paired output format for 'ACTGAC', k=4: ['ACT', 'CTG', 'TGA', 'GAC'] '''
        return [read[i:i+k-1] for i in range(len(read)-(k-2))]


class CMSDeBruijnGraph(DeBruijnGraph):
    def __init__(self, reads: list, k=None, hamming_dist=None, paired_error=None):
        # Set constants if provided, otherwise defaults to class values
        if k is not None:
            self.KMER_LEN = k
        if hamming_dist is not None:
            self.HAMMING_DIST = hamming_dist
        if paired_error is not None:
            self.ALLOWED_PAIRED_DIST_ERROR = paired_error

        self.num_edges = 0
        # Indexed by data/prefix_paired/suffix_paired.
        self.nodes = dict()
        kmer_counts_dict = self._count_kmers(self.KMER_LEN, reads)
        kmer_counts_sketch = self._make_sketch(kmer_counts_dict)
        del kmer_counts_dict
        self._build_graph(kmer_counts_sketch, reads)

    @staticmethod
    def _make_sketch(kmer_counts_dict: defaultdict) -> CountMinSketch:
        # Read the dictionary into a compressed data structure to allow deleting kmer_counts_dict
        NUM_ROWS = 10
        kmer_counts = CountMinSketch(NUM_ROWS)
        for kmer, count in kmer_counts_dict.items():
            kmer_counts.update(kmer, count)
        return kmer_counts


class PairedDeBruijnGraph(AbstractDeBruijnGraph):
    ''' Implementation relies on idea of dual-key with paired-reads.
        That is, self.nodes is a defaultdict of dicts; two keys are necessary to identify a node.
        A tuple would have sufficed if paired reads were drawn with a perfect distance,
        but since it's inexact, it's helpful to be able to iterate through the secondary dict
        for key comparison.
    '''

    # 28 for N. Delto, ?? for E. Coli
    KMER_LEN = 23
    HAMMING_DIST = 3
    ALLOWED_PAIRED_DIST_ERROR = 2

    def __init__(self, reads: list, k=None, hamming_dist=None, paired_error=None):
        # Set constants if provided, otherwise defaults to class values
        if k is not None:
            self.KMER_LEN = k
        if hamming_dist is not None:
            self.HAMMING_DIST = hamming_dist
        if paired_error is not None:
            self.ALLOWED_PAIRED_DIST_ERROR = paired_error

        self.num_edges = 0
        # Indexed by (data, paired_data)/(prefix, paired_prefix)/(suffix, paired_suffix)
        self.nodes = defaultdict(dict)
        kmer_counts_dict = self._count_kmers(self.KMER_LEN, reads)
        self._build_graph(kmer_counts_dict, reads)

    def enumerate_contigs(self) -> list:
        contigs = list()
        total_len = 0
        for _, nodes in self.nodes.items():
            for _, node in nodes.items():
                while node.outdegree > 0 and (node.was_branching or node.indegree == 0):          #
                    contig = self._get_longest_contig(node)
                    contigs.append(contig)
                    total_len += len(contig)
                if self.num_edges == 0:
                    return contigs

        # Handles when some sub contig is perfectly circular
        # All non-circular segments will have been found by iteration above.
        for _, nodes in self.nodes.items():
            for _, node in nodes.items():
                while node.outdegree > 0:          #
                    contig = self._get_longest_contig(node)
                    contigs.append(contig)
                    total_len += len(contig)
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

    def _build_graph(self, kmer_counts, reads: list):

        for i, read in enumerate(reads):
            read_pairs = self._break_read_into_k_minus_one_mers(self.KMER_LEN, read)
            for prefix_paired, suffix_paired in self._pairwise(read_pairs):
                # Indiscriminately filters any edge that has an erroneous k-1mer
                if kmer_counts[prefix_paired[0]] > self.HAMMING_DIST and \
                        kmer_counts[suffix_paired[0]] > self.HAMMING_DIST and \
                        kmer_counts[prefix_paired[1]] > self.HAMMING_DIST and \
                        kmer_counts[suffix_paired[1]] > self.HAMMING_DIST:

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
                        prefix_node.append_edge(
                            (sys.intern(suffix_node.data), sys.intern(suffix_node.paired_data)))
                        suffix_node.num_edges_in += 1
                        self.num_edges += 1
                    # Case: Both nodes found but edge between them doesn't exist.
                    elif (suffix_node.data, suffix_node.paired_data) not in prefix_node.edges:
                        prefix_node.append_edge(
                            (sys.intern(suffix_node.data), sys.intern(suffix_node.paired_data)))
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
    def _find_longest_overlap_brute(pattern: str, text: str) -> int:
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
    def _count_kmers(k: int, reads: list) -> defaultdict:
        ''' Generates a structure to check whether a kmer appears frequently enough.
            Does not save the k-1mers to trade runtime for space.
            Currently uses a dict but can switch over a count-min sketch if space is a concern.
        '''
        kmer_counts_dict = defaultdict(int)
        for i, cur_read_pair in enumerate(reads):
            # Need to count all kmers for error correction / filtering.
            paired_k_minus_one_mers = PairedDeBruijnGraph._break_read_into_k_minus_one_mers(
                k, cur_read_pair)

            # Count k-1mers instead of kmer of form (k-1, k-1) to save space.
            # The filter principle stays the same to remove errors.
            for paired_k_minus_one_mer in paired_k_minus_one_mers:
                kmer_counts_dict[paired_k_minus_one_mer[0]] += 1
                kmer_counts_dict[paired_k_minus_one_mer[1]] += 1

        return kmer_counts_dict

    @staticmethod
    def _break_read_into_k_minus_one_mers(k: int, read: str) -> list:
        ''' Paired output format for ('ACTGAC', 'TCGATC'), k=4: 
            [('ACT', 'TCG'), ('CTG', 'CGA'), ('TGA', 'GAT'), ('GAC', 'ATC')]
        '''
        return [(read[0][i:i+k-1], read[1][i:i+k-1]) for i in range(len(read[0])-(k-2))]


class CMSPairedDeBruijnGraph(PairedDeBruijnGraph):
    def __init__(self, reads: list, k=None, hamming_dist=None, paired_error=None):
        # Set constants if provided, otherwise defaults to class values
        if k is not None:
            self.KMER_LEN = k
        if hamming_dist is not None:
            self.HAMMING_DIST = hamming_dist
        if paired_error is not None:
            self.ALLOWED_PAIRED_DIST_ERROR = paired_error

        self.num_edges = 0
        # Indexed by (data, paired_data)/(prefix, paired_prefix)/(suffix, paired_suffix)
        self.nodes = defaultdict(dict)
        kmer_counts_dict = self._count_kmers(self.KMER_LEN, reads)
        kmer_counts_sketch = self._make_sketch(kmer_counts_dict)
        del kmer_counts_dict
        self._build_graph(kmer_counts_sketch, reads)

    @staticmethod
    def _make_sketch(kmer_counts_dict: defaultdict) -> CountMinSketch:
        # Read the dictionary into a compressed data structure to allow deleting kmer_counts_dict
        NUM_ROWS = 8
        kmer_counts = CountMinSketch(NUM_ROWS)
        for kmer, count in kmer_counts_dict.items():
            kmer_counts.update(kmer, count)
        return kmer_counts


""" ----- DEBUG GRAPHS ----- """


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


""" ----- IO HANDERS ----- """


class IOHandler:
    @staticmethod
    def read_args():
        parser = argparse.ArgumentParser(
            description="Generates contigs of a parent string given a set of substrings or \
                paired substrings and writes them to a file.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-t', '--time', action='store_true',
                            help="prints time at each stage during runtime to track program progression")
        parser.add_argument('-m', '--memory', action='store_true',
                            help="prints size of major data structures during runtime")
        parser.add_argument('-c', '--count_min_sketch', action='store_true',
                            help="uses a probabilistic data structure instead of a dictionary")

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
            and indications of which (paired?, d).        
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
    def write_to_FASTQ(contigs):
        # MODIFY THIS TO WRITE TO A FILE
        # Check for subfolder to insert
        # Check for file and append appropriate info
        # Append K value for K-mer, Hamming Distance, Allowed Error
        for i, contig in enumerate(contigs):
            print(">CONTIG" + str(i+1))
            print(contig)


class DebugIOHandler(IOHandler):
    @staticmethod
    def read_input(print_runtime=True, print_snapshot=False, print_syssizeof=False, start_time=0):
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

    return contigs


def assemble_without_options(k=None, hamming_dist=None, paired_error=None):
    reads, paired, _, _ = IOHandler.read_input()
    graph = None
    if paired:
        graph = PairedDeBruijnGraph(reads, k=k, hamming_dist=hamming_dist,
                                    paired_error=paired_error)
    else:
        graph = DeBruijnGraph(reads, k=k, hamming_dist=hamming_dist, paired_error=paired_error)
    contigs = graph.enumerate_contigs()
    return contigs


def main():
    start_time = time.time()
    args = IOHandler.read_args()

    print(args)

    contigs = list()
    if args.time or args.memory or args.count_min_sketch:
        print("WITH OPTIONS")
        contigs = assemble_with_options(
            print_runtime=args.time, print_syssizeof=args.memory,
            using_count_min_sketch=args.count_min_sketch, start_time=start_time,
            k=args.kmer_length, hamming_dist=args.filter_threshold, paired_error=args.error)
    else:
        print("WITHOUT OPTIONS")
        contigs = assemble_without_options(
            k=args.kmer_length, hamming_dist=args.filter_threshold, paired_error=args.error)
    IOHandler.write_to_FASTQ(contigs)

    if args.time:
        print(">--- PROGRAM FINISHED AT T = {:.2f} ---".format(time.time()-start_time))


if __name__ == "__main__":
    main()
    # cProfile.run('main()')
