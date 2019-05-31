import sys

from abc import ABC, abstractmethod
from collections import deque, defaultdict
from countminsketch import CountMinSketch
from debruijn_node import Node, PairedNode
from functools import lru_cache
from itertools import tee


class AbstractDeBruijnGraph(ABC):
    @abstractmethod
    def enumerate_contigs(self):
        pass

    @abstractmethod
    def _get_longest_contig(self, start_node):
        pass

    @abstractmethod
    def _build_graph(self, reads: list):
        pass

    @staticmethod
    @abstractmethod
    def _count_kmers(k: int, reads: list):
        pass

    @staticmethod
    @abstractmethod
    def _break_read_into_k_minus_one_mers(k: int, read: str):
        ''' Breaks into list of prefix_pairs and suffix_pairs. '''
        pass
    
    @staticmethod
    def valid_allowed_error(KMER_LEN, ALLOWED_PAIRED_DIST_ERROR):
        return KMER_LEN > ALLOWED_PAIRED_DIST_ERROR

    @staticmethod
    def _pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


class DeBruijnGraph(AbstractDeBruijnGraph):
    ''' Regular DeBruijn Graph with prefix-suffix pairs as edges. '''

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
        
        if not self.valid_allowed_error(self.KMER_LEN, self.ALLOWED_PAIRED_DIST_ERROR):
            raise ValueError("Allowed error must be less than the kmer length.")

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
        ''' Finds the longest contig by moving both forward from a branching node to another. 
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
        
        if not self.valid_allowed_error(self.KMER_LEN, self.ALLOWED_PAIRED_DIST_ERROR):
            raise ValueError("Allowed error must be less than the kmer length.")

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
        A tuple into a dict would have sufficed if paired reads were drawn with a perfect distance,
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
        
        if not self.valid_allowed_error(self.KMER_LEN, self.ALLOWED_PAIRED_DIST_ERROR):
            raise ValueError("Allowed error must be less than the kmer length.")

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
        
        if not self.valid_allowed_error(self.KMER_LEN, self.ALLOWED_PAIRED_DIST_ERROR):
            raise ValueError("Allowed error must be less than the kmer length.")

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
