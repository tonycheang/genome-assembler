#python3

import sys
from collections import deque, defaultdict
from itertools import tee, combinations

""" Implementing an Assembler

    -- Goal --

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. Instead of outputting an Eulerian path--unique
    construction can be rare--this assembler outputs contigs (non-branching segments).

    Data Sets Tested Against:
    N. deltocephalinicola                - t = 34000,   len = 100,    coverage: 30
    N. deltocephalinicola (read-pairs)   - t = 19679,   len = varies, coverage: unknown (likely around 30)
    E. Coli O104                         - t = 1400000, len = 100,    coverage: 25
    E. Coli O104          (read-pairs)   - t = 700000,  len = 100,    coverage: 25

    -- Approach --

    (1) Filter kmers with low coverage. (Hamming distance approach)
            Can use Count-min sketch to store counts using log(n) space if space becomes a concern.
            DETAILS OF FILTERING FOR ERRONEOUS READ-PAIRS?
            We can ignore bubbles and tips after this.
    (2) Build two nodes for each edge to create an unconnected graph.
    (3) Merge nodes (a|b) and (a|b') when b and b' overlap by more than length-DELTA, 
            where DELTA is the acceptable misalign distance. This creates a connected simple graph.
    (3) ENUMERATE CONTIGS AND THE JOINING SPOTS - HOW (start, ending, node)?
            When odd branching encountered or unbalanced node...
    (4) Build a graph where each edge is a contig that connects at the correct location to other contigs.
            DETAILS?
    (5) Run a max circulation flow through the network to infer multiplicity. 
            O(V[E**2]) doable for ~300 contigs even if this algorithm might be off by a factor of two or three.
            O(VE) algorithm exists if runtime becomes a concern.
"""


class DeBruijnGraph:
    """ Single-use class. If you enumerate_contigs, you will not be able to build_eulerian_path. """

    # Constants that control how the graph is built.
    KMER_LEN = 20
    HAMMING_DIST = 5
    # This is 2*DELTA where DELTA is the maximum disance from the true mean D, the distance between paried kmers.
    ALLOWED_PAIRED_DIST_ERROR = 6

    def __init__(self, reads, paired=False, d=0):
        self.num_edges = 0
        self.nodes = dict()         # Indexed by data/prefix/suffix.
        self.is_paired_graph = paired
        if self.is_paired_graph:
            self.paired_dist = d
            # We recieve half the number of reads with read-pairs
            DeBruijnGraph.HAMMING_DIST = DeBruijnGraph.HAMMING_DIST//2
        self._build_graph(reads)

    def enumerate_contigs(self):
        contigs = list()
        # Does a for-loop suffice? What if multi-edges allowed?
        for node_id in self.nodes:
            if self.nodes[node_id].out_degree > 0:
                found, contig = self._get_longest_contig(node_id)
                if found:
                    contigs.append(contig)
        return contigs

    def _get_longest_contig(self, start_node_id, visited):
        ''' Finds the longest contig by moving both forward and backward until
            nodes with branches are found.
        '''
        chain = deque()
        cur_node = self.nodes[start_node_id]

        # Only traverse backward if the starting node is non-branching.
        if cur_node.out_degree == 1:
            # HOW TO TRAVERSE BACKWARD???
            pass
        # If the node is branching, pick some edge.
        elif cur_node.out_degree > 1:
            chain.append(cur_node)
            cur_node_id = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_id]
        # Otherwise there are no outward edges. Let another node reach this one.
        else:
            return (False, None)

        # Traverse forward.
        while cur_node.out_degree == 1:
            chain.append(cur_node)
            cur_node_id = cur_node.pop_edge()
            cur_node = self.nodes[cur_node_id]

        # DOUBLE CHECK THIS. Feel like it's different.
        return (True, [node.id[-1] for node in chain])

    def build_eulerian_path(self):
        ''' Builds tours to save in class, then reconstructs based off build order. '''
        self.tours_with_available_edges = list()    # Tours used to build Eulerian path.
        self.tours_in_order = deque()

        built = self.__build_all_tours()
        if not built:
            return (False, None)

        eulerian_path = list()
        unfinished_tours = list()

        cur_tour = self.tours_in_order.popleft()

        # Handles single tour case
        if not self.tours_in_order:
            eulerian_path = [node.id[-1] for node in cur_tour.path]
            return (True, eulerian_path)

        next_tour = self.tours_in_order.popleft()

        while self.tours_in_order or unfinished_tours or cur_tour.has_remaining_path():
            # Either exhaust the current tour or find the starting point of the next (-1 if no next)
            while cur_tour.path and (next_tour == -1 or cur_tour.peek_next_path_node() != next_tour.start_v):
                cur_node = cur_tour.get_next_path_node()
                eulerian_path.append(cur_node.id[-1])

            # If we didn't exhaust the current tour, then record it and then start the next one
            if cur_tour.has_remaining_path():
                unfinished_tours.append(cur_tour)

            # If the next one is meant to join the current one, then move onto the next one
            if next_tour != -1 and next_tour.tour_to_join == cur_tour:
                cur_tour = next_tour
                # Handles when the cur_tour is the last tour
                if self.tours_in_order:
                    next_tour = self.tours_in_order.popleft()
                else:
                    next_tour = -1
            elif unfinished_tours:
                # Otherwise, work on the most recent unfinished one and don't change next_tour
                cur_tour = unfinished_tours.pop()

        return (True, eulerian_path)

    def __build_all_tours(self):
        ''' Exhausts all available edges building tours and keeps track of build order. '''

        prev_tour = -1
        start_v = self.nodes[next(iter(self.nodes))]

        while self.num_edges > 0:

            found, tour = self.__build_tour(start_v, prev_tour)
            if not found:
                return False
            self.tours_in_order.append(tour)

            if tour.has_available_edges():
                self.tours_with_available_edges.append(tour)
                prev_tour = tour
            elif self.tours_with_available_edges:
                prev_tour = self.tours_with_available_edges.pop()

                # Edges may be used up in the construction of other tours
                while self.tours_with_available_edges and not prev_tour.has_available_edges():
                    prev_tour = self.tours_with_available_edges.pop()

                # May have more edges even after using the next one, so add back on top
                if prev_tour.has_available_edges():
                    self.tours_with_available_edges.append(prev_tour)

            else:
                assert self.num_edges == 0, "Should have no more edges"
                break

            if prev_tour != -1 and prev_tour.has_available_edges():
                start_v = prev_tour.get_node_with_available_edge()
        return True

    def __build_tour(self, start_v, tour_to_join):
        ''' Constructs a single tour from remaining edges if possible.
            Returns (False, None) if not possible.
        '''
        circular_path = Tour(start_v, tour_to_join)
        cur_v = self.nodes[start_v.id]
        started = False

        while cur_v != start_v or not started:
            started = True
            # Check for available edges out from the current node
            if self.nodes[cur_v.id].edges:
                next_v_id = self.nodes[cur_v.id].pop_edge()
                self.num_edges -= 1
                circular_path.add_node(cur_v)
                # If the next vertex is the starting one, we don't add it to our path
                cur_v = self.nodes[next_v_id]
            else:
                # Executes when we hit a dead end (including empty adj dict)
                # that's not the starting vertex
                return (False, None)
        return (True, circular_path)

    def _build_graph(self, reads):
        ''' Builds the path constructor with coverage considerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        if self.is_paired_graph:
            edge_counts = DeBruijnGraph._count_paired_edges(reads)
            self._connect_paired_graph(edge_counts)
        else:
            edge_counts = DeBruijnGraph._count_edges(reads)
            self._connect_graph(edge_counts)

    def _connect_graph(self, edge_counts):
        for edge, count in edge_counts.items():
            prefix, suffix = edge[0], edge[1]
            if count > DeBruijnGraph.HAMMING_DISTANCE:
                if prefix not in self.nodes:
                    self.nodes[prefix] = Node(prefix)
                if suffix not in self.nodes:
                    self.nodes[suffix] = Node(suffix)

                self.nodes[prefix].add_edge(suffix)
                self.nodes[suffix].add_reverse_edge(prefix)
                self.num_edges += 1

    @staticmethod
    def _count_edges(reads):
        def pairwise(iterable):
            "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

        edge_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = DeBruijnGraph._break_read_into_k_minus_one_mers(
                DeBruijnGraph.KMER_LEN, cur_read, paired=False)

            for prefix, suffix in pairwise(k_minus_one_mers):
                edge_counts[(prefix, suffix)] += 1
        return edge_counts

    def _connect_paired_graph(self, edge_counts):
        pass

    @staticmethod
    def _count_paired_edges(reads):
        def pairwise(iterable):
            "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

        return list()

    @staticmethod
    def _break_read_into_k_minus_one_mers(k, read, paired=False):
        ''' Breaks into list of prefixes and suffixes. 
            Non-paired output format for 'ACTGAC', k=4: ['ACT', 'CTG', 'TGA', 'GAC']
            Paired output format for ('ACTGAC', 'TCGATC' , 3), k=4: 
                [('ACT', 'TCG', 3), ('CTG', 'CGA', 3), ('TGA', 'GAT', 3), ('GAC', 'ATC', 3)]
        '''
        if paired:
            return [(read[0][i:i+k-1], read[1][i:i+k-1]) for i in range(len(read[0])-(k-2))]
        else:
            return [read[i:i+k-1] for i in range(len(read)-(k-2))]


class Node:
    def __init__(self, id_, paired_data=None, d_to_pair=None):
        # PAIRED DATA MIGHT NEED UNQIUE ID FOR HASH TABLE? DEPENDS ON MERGING...
        # ID contains the prefix/suffix string
        self.id = id_
        # Edges are the string, or the key into the constructor
        self.edges = list()
        self.reverse_edges = list()
        self.tours = list()
        if paired_data is not None:
            self.paired_data = paired_data
            self.d_to_pair = d_to_pair

    @property
    def out_degree(self):
        return len(self.edges)

    @property
    def in_degree(self):
        return len(self.reverse_edges)

    def pop_edge(self):
        ''' Gets and removes an edge. '''
        edge = self.edges.pop()
        if not self.edges:
            self.notify_tours_no_available_edges()
        return edge
    
    def merge(self, other):
        pass

    def has_available_edges(self):
        return bool(self.edges)

    def add_edge(self, edge):
        self.edges.append(edge)

    def add_reverse_edge(self, edge):
        self.reverse_edges.append(edge)

    def notify_tours_no_available_edges(self):
        ''' Used to tell tours not to keep the node in its avaible_edges dict '''
        for tour in self.tours:
            tour.remove_node_from_available_edges(self)

    def subscribe_tour(self, tour):
        self.tours.append(tour)


class Tour:
    def __init__(self, start, tour_to_join):
        self.tour_to_join = tour_to_join
        self.start_v = start
        self.path = deque()
        self.available_edges = dict()

    def add_node(self, node):
        self.path.append(node)
        if node.has_available_edges() and not node in self.available_edges:
            node.subscribe_tour(self)
            self.available_edges[node] = node.edges

    def remove_node_from_available_edges(self, node):
        del self.available_edges[node]

    def has_available_edges(self):
        return bool(self.available_edges)

    def has_remaining_path(self):
        return bool(self.path)

    def get_node_with_available_edge(self):
        return next(iter(self.available_edges))

    def get_next_path_node(self):
        return self.path.popleft()

    def peek_next_path_node(self):
        return self.path[0]


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
    # graph = DeBruijnGraph(reads, paired, d)
    # contigs = graph.enumerate_contigs()
    # IOHandler.print_FASTA(contigs)


if __name__ == "__main__":
    main()
