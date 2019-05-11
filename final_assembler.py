#python3

import sys
from abc import ABC, abstractmethod
from collections import deque, defaultdict
from itertools import tee, combinations

""" Implementing an Assembler

    --- Goal ---

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. Instead of outputting an Eulerian path--unique
    construction can be rare--this assembler outputs contigs (non-branching segments).

    Data Sets Tested Against:
    N. deltocephalinicola                - t = 34000,   len = 100,    coverage: 30
    N. deltocephalinicola (read-pairs)   - t = 19679,   len = varies, coverage: unknown (likely around 30)
    E. Coli O104                         - t = 1400000, len = 100,    coverage: 25
    E. Coli O104          (read-pairs)   - t = 700000,  len = 100,    coverage: 25

    --- Approach ---

        - Building the Graph -

    (1) Filter kmers with low coverage. (Hamming distance approach)
            Can use Count-min sketch to store counts using log(n) space if space becomes a concern.
            DETAILS OF FILTERING FOR ERRONEOUS READ-PAIRS?
            We can ignore bubbles and tips after this.
    (2) Build two nodes for each edge to create an unconnected graph.
    (3) Merge nodes (a|b) and (a|b') when b and b' overlap by more than length 2*DELTA, 
            where DELTA is the misalign distance. This creates a connected simple graph.
    (4) Each node at this point must record whether they are a branching node or not.
            Otherwise, once the last edges remain from popping, it will be counted as part of another contig.
    
        - Building Contigs -
    
    (5) Traverse each edge forward and backward to obtain a non-branching segment. This is a contig.
            On branching nodes, record where each contig joins another.
            The branching nodes with contigs as edges will make up the circulation flow graph.
    
        - Infering Multiplicity of Contigs -
        (Skip for now. Infer only if necessary. Contigs of repeated regions will be unrepresented.)
    
    (6) Build a graph where each edge is a contig that connects at the correct location to other contigs.
            Long contigs (>1000 BP) will have max and min flow of 1. They are correct.
    (7) Run a max circulation flow through the network to infer multiplicity. 
            O(V[E**2]) doable for ~300 contigs even if this algorithm might be off by a factor of two or three.
            O(VE) algorithm exists if runtime becomes a concern.
"""


class Node:
    def __init__(self, data):
        # data contains the prefix_paired/suffix_paired string
        self.data = data
        # Edges are the string, or the key into the constructor
        self.edges = list()
        self.reverse_edges = list()
        self.branching = False

    @property
    def out_degree(self):
        return len(self.edges)

    @property
    def in_degree(self):
        return len(self.reverse_edges)

    @property
    def has_available_edges(self):
        return bool(self.edges)

    def pop_edge(self):
        ''' Gets and removes an edge. '''
        edge = self.edges.pop()
        if not self.edges:
            self.notify_tours_no_available_edges()
        return edge

    def append_edge(self, edge):
        self.edges.append(edge)

    def append_reverse_edge(self, edge):
        self.reverse_edges.append(edge)


class PairedNode(Node):
    def __init__(self, paired_data, pair_dist, data):
        self.paired_data = paired_data
        self.pair_dist = pair_dist
        super().__init__(data)

    # Consider whether this goes here.
    @staticmethod
    def overlaps(text, pattern, threshold):
        pass


class AbstractDeBruijnGraph(ABC):
    @abstractmethod
    def enumerate_contigs(self):
        pass

    @abstractmethod
    def __get_longest_contig(self, start_node_data, visited):
        pass

    @abstractmethod
    def __build_graph(self, reads):
        pass

    @staticmethod
    @abstractmethod
    def __count_kmers(reads):
        pass

    @staticmethod
    @abstractmethod
    def __break_read_into_k_minus_one_mers(k, read):
        ''' Breaks into list of prefix_pairedes and suffix_pairedes. '''
        pass

    @staticmethod
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


class DeBruijnGraph(AbstractDeBruijnGraph):

    KMER_LEN = 20
    HAMMING_DIST = 5
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D, the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads):
        self.num_edges = 0
        # Indexed by data/prefix_paired/suffix_paired.
        self.nodes = dict()
        self.__build_graph(reads)

    def enumerate_contigs(self):
        contigs = list()
        # Does a for-loop suffice? What if multi-edges allowed?
        for node_data in self.nodes:
            if self.nodes[node_data].out_degree > 0:
                found, contig = self._get_longest_contig(node_data)
                if found:
                    contigs.append(contig)
        return contigs

    def __get_longest_contig(self, start_node_data, visited):
        ''' Finds the longest contig by moving both forward and backward until
            nodes with branches are found.
        '''
        chain = deque()
        cur_node = self.nodes[start_node_data]

        # Only traverse backward if the starting node is non-branching.
        if cur_node.out_degree == 1:
            # HOW TO TRAVERSE BACKWARD???
            pass
        # If the node is branching, pick some edge.
        elif cur_node.out_degree > 1:
            chain.append(cur_node)
            cur_node_data = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_data]
        # Otherwise there are no outward edges. Let another node reach this one.
        else:
            return (False, None)

        # Traverse forward.
        while cur_node.out_degree == 1:
            chain.append(cur_node)
            cur_node_data = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_data]

        # DOUBLE CHECK THIS. Feel like it's different.
        return (True, [node.data[-1] for node in chain])

    def __build_graph(self, reads):
        ''' Builds the path constructor with coverage consdataerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        kmer_counts = DeBruijnGraph.__count_kmers(reads)
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
    def __count_kmers(reads):

        kmer_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = DeBruijnGraph._break_read_into_k_minus_one_mers(
                DeBruijnGraph.KMER_LEN, cur_read)

            for prefix, suffix in AbstractDeBruijnGraph.pairwise(k_minus_one_mers):
                kmer_counts[(prefix, suffix)] += 1
        return kmer_counts

    @staticmethod
    def __break_read_into_k_minus_one_mers(k, read, paired=False):
        ''' Non-paired output format for 'ACTGAC', k=4: ['ACT', 'CTG', 'TGA', 'GAC'] '''
        return [read[i:i+k-1] for i in range(len(read)-(k-2))]


class PairedDeBruijnGraph(AbstractDeBruijnGraph):

    KMER_LEN = 3  # 20
    HAMMING_DIST = 3
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D, the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads, d):
        self.num_edges = 0
        # Indexed by (data, paired_data)/(prefix, paired_prefix)/(suffix, paired_suffix)
        self.nodes = dict()
        self.dist = d
        self._build_graph(reads)

    def enumerate_contigs(self):
        pass

    def __get_longest_contig(self):
        pass

    def __build_graph(self, reads):
        kmer_counts, broken_read_pairs = PairedDeBruijnGraph.__count_kmers(
            reads)

        for read_pairs in broken_read_pairs:
            for prefix_paired, suffix_paired in AbstractDeBruijnGraph.pairwise(read_pairs):
                # Indiscriminately filters any edge that has an erroneous k-mer
                if kmer_counts[(prefix_paired[0], suffix_paired[0])] > PairedDeBruijnGraph.HAMMING_DIST and \
                        kmer_counts[(prefix_paired[1], suffix_paired[1])] > PairedDeBruijnGraph.HAMMING_DIST:

                    if prefix_paired not in self.nodes:
                        self.nodes[prefix_paired] = Node(
                            prefix_paired[0], prefix_paired[1])
                    if suffix_paired not in self.nodes:
                        self.nodes[suffix_paired] = Node(
                            suffix_paired[0], suffix_paired[1])

                    self.nodes[prefix_paired].append_edge(suffix_paired)
                    self.nodes[suffix_paired].append_reverse_edge(
                        prefix_paired)
                    self.num_edges += 1

    def __merge_nodes(self, node, node_to_remove):
        pass

    @staticmethod
    def __count_kmers(reads):
        ''' Generates a structure to check whether a kmer appears frequently enough.
            Currently uses a dict but can switch over a count-min sketch if space is a concern.
        '''

        kmer_counts = defaultdict(int)
        broken_read_pairs = list()
        for cur_read_pair in reads:
            paired_k_minus_one_mers = PairedDeBruijnGraph.__break_read_into_k_minus_one_mers(
                PairedDeBruijnGraph.KMER_LEN, cur_read_pair)

            broken_read_pairs.append(paired_k_minus_one_mers)

            for prefix_paired_paired, suffix_paired_paired in AbstractDeBruijnGraph.pairwise(paired_k_minus_one_mers):
                # Each kmer_pair generates two kmers, each of which we'll count for coverage.
                kmer_counts[(prefix_paired_paired[0],
                             suffix_paired_paired[0])] += 1
                kmer_counts[(prefix_paired_paired[1],
                             suffix_paired_paired[1])] += 1
        return (kmer_counts, broken_read_pairs)

    @staticmethod
    def __break_read_into_k_minus_one_mers(k, read):
        ''' Paired output format for ('ACTGAC', 'TCGATC'), k=4: 
            [('ACT', 'TCG'), ('CTG', 'CGA'), ('TGA', 'GAT'), ('GAC', 'ATC')]
        '''
        return [(read[0][i:i+k-1], read[1][i:i+k-1]) for i in range(len(read[0])-(k-2))]


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
            print(">CONTIG" + str(i))
            print(contig)


def main():
    reads, paired, d = IOHandler.read_input()
    graph = None
    if paired:
        graph = PairedDeBruijnGraph(reads, d)
    else:
        graph = DeBruijnGraph(reads)
    # contigs = graph.enumerate_contigs()
    # IOHandler.print_FASTA(contigs)


if __name__ == "__main__":
    main()
