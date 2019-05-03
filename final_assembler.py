#python3

import sys
from collections import deque, defaultdict
from itertools import tee, combinations

sys.setrecursionlimit(10**6)

""" Implementing an Assembler

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. Instead of outputting an Eulerian path--unique
    construction can be rare--this assembler outputs contigs (non-branching segments).

    Data Sets Tested Against:
    N. deltocephalinicola                - t = 34000,   len = 100,    coverage: 30
    N. deltocephalinicola (read-pairs)   - t = 19679,   len = varies, coverage: unknown (likely around 30)
    E. Coli O104                         - t = 1400000, len = 100,    coverage: 25
    E. Coli O104          (read-pairs)   - t = 700000,  len = 100,    coverage: 25

    -- Approach --
    Filter kmers with low coverage (hamming distance approach)
    We can ignore bubbles and tips after this.

    TO DO:
    Modify DeBruijnGraph to return contigs instead of whole genome.
        Traverse backward until find node with branching, then take previous (non branching start)
        Where to start to ensure long contigs?
    Figure out a filter threshold (likely still 5, maybe 4 for E. Coli, less for read pairs or not?).
        How to account for multiplicity? Can we ignore or not?
        Query statistics?
    Modify Node to optionally take read-pair and distance d
    Modify readInput to check for read-pair or regular read
    Figure how to merge unnecessary edges. (do edges need a new representation? or nodes...?)
"""

# ... Do I actually just scrap most of this class?


class DeBruijnGraph:
    KMER_LEN = 20
    HAMMING_DIST = 5

    def __init__(self, reads):
        self.num_edges = 0
        # Indexed by ID or prefix/suffix
        self.nodes = dict()
        self.tours_with_available_edges = list()
        self.unique_edge_counts = defaultdict(int)
        self.tours_in_order = deque()
        self._build_graph(reads)

    # Use following methods as reference to build contigs, since it a similar process
    # In fact, might be easier since we don't have to join the contigs.
    # Just might want to travel upstream...

    def build_eulerian_path(self):
        ''' Builds tours to save in class, then reconstructs based off build order. '''

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
                next_v_id = self.nodes[cur_v.id].get_edge()
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
        def pairwise(iterable):
            "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

        edge_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = DeBruijnGraph._break_read_into_k_minus_one_mers(
                DeBruijnGraph.KMER_LEN, cur_read)

            for prefix, suffix in pairwise(k_minus_one_mers):
                edge_counts[(prefix, suffix)] += 1

        for edge, count in edge_counts.items():
            prefix, suffix = edge[0], edge[1]
            if count > DeBruijnGraph.HAMMING_DISTANCE:
                if prefix not in self.nodes:
                    self.nodes[prefix] = Node(prefix)
                if suffix not in self.nodes:
                    self.nodes[suffix] = Node(suffix)

                self.nodes[prefix].add_edge(suffix)
                self.nodes[suffix].in_degree += 1
                self.num_edges += 1

    @staticmethod
    def _break_read_into_k_minus_one_mers(k, read):
        ''' Breaks into list of prefixes and suffixes. '''
        return [read[i:i+k-1] for i in range(len(read)-(k-2))]


class Node:
    def __init__(self, id_):
        # ID contains the prefix/suffix string
        self.id = id_
        # Edges are the string, or the key into the constructor
        self.edges = list()
        self.tours = list()
        self.in_degree = 0

    @property
    def out_degree(self):
        return len(self.edges)

    def get_edge(self):
        ''' Gets and removes an edge. '''
        edge = self.edges.pop()
        if not self.edges:
            self.notify_tours_no_available_edges()
        return edge

    def has_available_edges(self):
        return bool(self.edges)

    def notify_tours_no_available_edges(self):
        ''' Used to tell tours not to keep the node in its avaible_edges dict '''
        for tour in self.tours:
            tour.remove_node_from_available_edges(self)

    def subscribe_tour(self, tour):
        self.tours.append(tour)

    def add_edge(self, edge):
        self.edges.append(edge)


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


def read_input():
    ''' Returns a list of reads, not yet broken down. '''
    reads = list()
    for read in sys.stdin:
        reads.append(read.split()[0])
    return reads


def format_output(contigs):
    pass


def main():
    reads = read_input()
    graph = DeBruijnGraph(reads)
    # found, path = graph.build_eulerian_path()
    # if found:
    #     print("".join(map(str, path)))
    # else:
    #     print("NOT FOUND")


if __name__ == "__main__":
    main()
