# python3

import sys
from abc import ABC, abstractmethod
from collections import deque, defaultdict
from itertools import tee, combinations
import time

import cProfile
import tracemalloc

# sys.setswitchinterval(10000)

""" Implementing an Assembler

    --- Goal ---

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. Instead of outputting an Eulerian path--unique
    construction can be rare--this assembler outputs contigs (non-was_branching segments).

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
    
    (5) Traverse each edge forward and backward to obtain a non-was_branching segment. This is a contig.
            If starting at a branch (i.e. in/outdegree > 1), take any node forward and only go 
                forward until next branch.
            If starting at a middle section, take nodes forward and backward until a was_branching node.
            If starting at a dead end, work backward until a was_branching node.

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

    --- Preliminary Test Results ---
    (Time from time module and memory from Activity Monitor)
    For 100,000 reads of length 100 each from N. Deltocephalinicola:
        STARTING ASSEMBLY PROGRAM AT T= 0
        FINISHED READING INPUT AT T= 1.0616950988769531
        STARTING TO BUILD GRAPH AT T= 1.06174898147583
        FINISHED COUNTING KMERS AT T= 13.7449312210083
        FINISHED ADDING EDGES AT T= 146.9288830757141
        FINISHED BUILDING GRAPH AT T= 146.99736309051514
        STARTING TO ENUMERATE CONTIGS AT T= 147.64608526229858
        FINISHED ENUMERATING CONTIGS AT T= 147.82946515083313
        PROGRAM FINISHED AT T= 147.83732104301453
        Memory Consumption: ~2 GB
    For 200,000 reads of length 100 each from N. Deltocephalinicola:
        STARTING ASSEMBLY PROGRAM AT T= 0
        FINISHED READING INPUT AT T= 1.9468772411346436
        STARTING TO BUILD GRAPH AT T= 1.9469482898712158
        FINISHED COUNTING KMERS AT T= 29.306518077850342
        FINISHED ADDING EDGES AT T= 300.31926918029785
        FINISHED BUILDING GRAPH AT T= 300.39726734161377
        STARTING TO ENUMERATE CONTIGS AT T= 301.8081531524658
        FINISHED ENUMERATING CONTIGS AT T= 302.0127410888672
        PROGRAM FINISHED AT T= 302.01618099212646
        Memory Consumption: ~3.8 GB
    For 200,000 reads of length 100 each from N. Deltocephalinicola:
    (With slight modification to early exit with hash! Would some associative container help too?)
        STARTING ASSEMBLY PROGRAM AT T= 0
        FINISHED READING INPUT AT T= 1.9912371635437012
        STARTING TO BUILD GRAPH AT T= 1.9912981986999512
        FINISHED COUNTING KMERS AT T= 26.330377101898193
        FINISHED BUILDING GRAPH AT T= 265.23331904411316
        STARTING TO ENUMERATE CONTIGS AT T= 266.47291707992554
        FINISHED ENUMERATING CONTIGS AT T= 266.6549780368805
        PROGRAM FINISHED AT T= 266.6550042629242
    For 400,000 reads of length 100 each from N. Deltocephalinicola:
        STARTING ASSEMBLY PROGRAM AT T= 0
        FINISHED READING INPUT AT T= 4.09502387046814
        STARTING TO BUILD GRAPH AT T= 4.0950767993927
        FINISHED COUNTING KMERS AT T= 58.88226580619812
        FINISHED BUILDING GRAPH AT T= 560.3410959243774
        STARTING TO ENUMERATE CONTIGS AT T= 565.3487329483032
        FINISHED ENUMERATING CONTIGS AT T= 565.5448567867279
        PROGRAM FINISHED AT T= 565.5449018478394
        Memory Consumption: ~7.4 GB

    Idea for saving memory: don't store backward edges. Rely on forward starts on branching nodes.
    What about processing each read through the graph as it comes in instead of storing an array?
    CountMin Sketch to estimate counts of each kmer to store in O(log(n)) space
        Profile to make sure the counts array is actually large...
    Debug class for DeBruijn Graph using multiple inheritance?
"""


class Node:
    def __init__(self, data):
        # data contains the prefix_paired/suffix_paired string
        self.data = data
        # Edges are the string, or the key into the constructor
        self.edges = dict()
        self.was_branching = False
        self.num_edges_in = 0

    @property
    def outdegree(self):
        return len(self.edges)

    @property
    def indegree(self):
        return self.num_edges_in

    @property
    def has_available_edges(self):
        return bool(self.edges)

    @property
    def has_available_reverse_edges(self):
        return bool(self.reverse_edges)

    def pop_edge(self):
        ''' Gets and removes a forward edge. '''
        return self.edges.popitem()

    def append_edge(self, edge):
        self.edges[edge] = True

    def __repr__(self):
        ''' For debugging small examples. '''
        return "Data: {0} | Pair: {1} | Edges: {2} | Reverse Edges: {3}" \
            .format(self.data, self.paired_data, self.edges, self.reverse_edges)


class PairedNode(Node):
    def __init__(self, data, paired_data):
        self.paired_data = paired_data
        super().__init__(data)


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
        ''' Breaks into list of prefix_pairedes and suffix_pairedes. '''
        pass

    @staticmethod
    def _pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


class DeBruijnGraph(AbstractDeBruijnGraph):

    KMER_LEN = 20
    HAMMING_DIST = 5
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D,
    # the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads):
        self.num_edges = 0
        # Indexed by data/prefix_paired/suffix_paired.
        self.nodes = dict()
        self._build_graph(reads)

    def enumerate_contigs(self):
        contigs = list()
        for node_data in self.nodes:
            if self.nodes[node_data].outdegree > 0:
                found, contig = self._get_longest_contig(node_data)
                if found:
                    contigs.append(contig)
        return contigs

    def _get_longest_contig(self, start_node_data, visited):
        ''' Finds the longest contig by moving both forward and backward until
            nodes with branches are found.
        '''
        chain = deque()
        cur_node = self.nodes[start_node_data]

        # Only traverse backward if the starting node is non-branching.
        if cur_node.outdegree == 1:
            # HOW TO TRAVERSE BACKWARD???
            pass
        # If the node is branching, pick some edge.
        elif cur_node.outdegree > 1:
            chain.append(cur_node)
            cur_node_data = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_data]
        # Otherwise there are no outward edges. Let another node reach this one.
        else:
            return (False, None)

        # Traverse forward.
        while cur_node.outdegree == 1:
            chain.append(cur_node)
            cur_node_data = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_data]

        # DOUBLE CHECK THIS. Feel like it's different.
        return (True, [node.data[-1] for node in chain])

    def _build_graph(self, reads):
        ''' Builds the path constructor with coverage considerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        kmer_counts = DeBruijnGraph._count_kmers(reads)
        for edge, count in kmer_counts.items():
            prefix, suffix = edge[0], edge[1]
            if count > DeBruijnGraph.HAMMING_DISTANCE:
                if prefix not in self.nodes:
                    self.nodes[prefix_paired] = Node(prefix)
                if suffix not in self.nodes:
                    self.nodes[suffix] = Node(suffix)

                self.nodes[prefix].append_edge(suffix)
                self.nodes[suffix].append_reverse_edge(prefix)
                self.num_edges += 1

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
        return [read[i:i+k-1] for i in range(len(read)-(k-2))]


class PairedDeBruijnGraph(AbstractDeBruijnGraph):

    KMER_LEN = 27
    HAMMING_DIST = 6
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D, the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads, d):
        self.num_edges = 0
        # Indexed by (data, paired_data)/(prefix, paired_prefix)/(suffix, paired_suffix)
        self.nodes = defaultdict(dict)
        self.dist = d
        kmer_counts, broken_read_pairs = PairedDeBruijnGraph._count_kmers(reads)
        del reads
        self._build_graph(kmer_counts, broken_read_pairs)

    def enumerate_contigs(self):
        contigs = list()
        for _, nodes in self.nodes.items():
            for _, node in nodes.items():
                while node.outdegree > 0 and (node.was_branching or node.indegree == 0):          #
                    contig = self._get_longest_contig(node)
                    contigs.append(contig)
                if self.num_edges == 0:
                    return contigs
        return contigs

    def _get_longest_contig(self, cur_node):
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

    def _build_graph(self, kmer_counts, broken_read_pairs):

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
                            data=prefix_paired[0], paired_data=prefix_paired[1])
                        self.nodes[prefix_paired[0]][prefix_paired[1]] = prefix_node
                    if not suffix_found:
                        suffix_node = PairedNode(
                            data=suffix_paired[0], paired_data=suffix_paired[1])
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

    def _find_matching_node(self, paired_strings):
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
    def _find_longest_overlap_brute(pattern, text):
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
    def _count_kmers(reads):
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
    def _break_read_into_k_minus_one_mers(k, read):
        ''' Paired output format for ('ACTGAC', 'TCGATC'), k=4: 
            [('ACT', 'TCG'), ('CTG', 'CGA'), ('TGA', 'GAT'), ('GAC', 'ATC')]
        '''
        return [(read[0][i:i+k-1], read[1][i:i+k-1]) for i in range(len(read[0])-(k-2))]


class DebugPairedDeBruijnGraph(PairedDeBruijnGraph):
    def __init__(self, reads, d, print_memory=False, print_runtime=False, start_time=0):
        self.print_memory = print_memory
        self.print_runtime = print_runtime
        self.start_time = start_time
        super().__init__(reads, d)

    def enumerate_contigs(self):
        if self.print_runtime:
            print("STARTING TO ENUMERATE CONTIGS AT T=", time.time()-self.start_time)
        contigs = super().enumerate_contigs()
        if self.print_runtime:
            print("FINISHED ENUMERATING CONTIGS AT T=", time.time()-self.start_time)
        if self.print_memory:
            print_memory_snapshot("AFTER ENUMERATE CONTIGS")
        return contigs

    def _build_graph(self, kmer_counts, broken_read_pairs):
        if self.print_runtime:
            print("STARTING TO BUILD GRAPH AT T=", time.time()-self.start_time)
        super()._build_graph(kmer_counts, broken_read_pairs)
        if self.print_runtime:
            print("FINISHED BUILDING GRAPH AT T=", time.time()-self.start_time)
        if self.print_memory:
            print("SIZE OF GRAPH: ", sys.getsizeof(self.nodes))
            print_memory_snapshot("AFTER BUILDING GRAPH")

    def _count_kmers(self, reads):
        kmer_counts, broken_read_pairs = super()._count_kmers(reads)
        if self.print_memory:
            count_mem = sys.getsizeof(kmer_counts)
            for pair in kmer_counts.items():
                count_mem += sys.getsizeof(pair)
                count_mem += sys.getsizeof(pair[0])
                count_mem += sys.getsizeof(pair[1])
            print("SIZE OF COUNTS: ", count_mem)
            total_mem = sys.getsizeof(broken_read_pairs)
            for read_pairs in broken_read_pairs:
                total_mem += sys.getsizeof(read_pairs)
                for pair in read_pairs:
                    total_mem += sys.getsizeof(pair)
                    total_mem += sys.getsizeof(pair[0])
                    total_mem += sys.getsizeof(pair[1])
            print("SIZE OF TOTAL K-1MER-PAIRS: ", total_mem)
            print_memory_snapshot("AFTER COUNTING KMERS")
        if self.print_runtime:
            print("FINISHED COUNTING KMERS AT T=", time.time()-self.start_time)
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
    def read_input(print_runtime=True, print_memory=False, start_time=0):
        if print_runtime:
            print("STARTING ASSEMBLY PROGRAM AT T=", 0)

        reads, paired, d = IOHandler.read_input()

        if print_runtime:
            print("FINISHED READING INPUT AT T=", time.time()-start_time)
        if print_memory:
            print("SIZE OF READS: ", sys.getsizeof(reads))
            total_string_memory = 0
            for read in reads:
                total_string_memory += sys.getsizeof(read)
            print("SIZE OF ALL READS: ", total_string_memory)
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


def profile_assembler(print_runtime=False, print_memory=False):

    if print_memory:
        tracemalloc.start(25)
    start_time = time.time()

    reads, paired, d = DebugIOHandler.read_input(
        print_runtime=print_runtime, print_memory=print_memory, start_time=start_time)
    graph = None
    if paired:
        graph = DebugPairedDeBruijnGraph(reads=reads, d=d, print_runtime=print_runtime,
                                         start_time=start_time, print_memory=print_memory)
    else:
        graph = DeBruijnGraph(reads)
    contigs = graph.enumerate_contigs()
    DebugIOHandler.print_FASTA(contigs)

    if print_runtime:
        print("PROGRAM FINISHED AT T=", time.time()-start_time)


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
    main()
    # profile_assembler(print_runtime=False, print_memory=False)
    # cProfile.run('main()')
    # cProfile.run('profile_assembler(print_runtime=True, print_memory=False)')
