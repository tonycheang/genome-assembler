#python3

import sys
from collections import deque, defaultdict
from itertools import tee, combinations

sys.setrecursionlimit(10**6)

""" Implementing an Assembler

    The assembler must handle both regular reads and read-pairs. For the latter problem,
    a paired de Bruijn graph can be defined. 

    Data Sets Tested Against:
    N. deltocephalinicola                - t = 34000,   len = 100,    coverage: 30
    N. deltocephalinicola (read-pairs)   - t = 19679,   len = varies, coverage: unknown (likely around 30)
    E. Coli O104                         - t = 1400000, len = 100,    coverage: 25
    E. Coli O104          (read-pairs)   - t = 700000,  len = 100,    coverage: 25

    -- Approach --
    Filter kmers with low coverage (hamming distance approach)
    We can ignore bubbles and tips after this.

    TO DO:
    Modify EulerianPathConstructor to return contigs instead of whole genome. (Contig Enumerator?)
        Contigs are continuous, non-branching segments.
        Where to start to ensure long contigs?
    Figure out a good value for k to break reads down into.
        15 - 20 is used by msot assemblers
    Figure out a filter threshold (likely still 5, maybe 4 for E. Coli, less for read pairs or not?).
        How to account for multiplicity? Can we ignore or not?
        Query statistics?
    Modify Node to optionally take read-pair and distance d
    Figure how to merge unnecessary edges. (do edges need a new representation? or nodes...?)
"""


class EulerianPathConstructor:
    ''' Single use class since we don't save the graph permanently '''

    def __init__(self):
        self.num_edges = 0
        # Indexed by ID or prefix/suffix
        self.nodes = dict()
        self.tours_with_available_edges = list()
        self.unique_edge_counts = defaultdict(int)
        self.tours_in_order = deque()

    def has_node(self, id_):
        return id_ in self.nodes

    def add_node(self, node):
        self.nodes[node.id] = node

    def add_edges_from_k_minus_one_mers(self, k_minus_one_mers, reverse=False):
        ''' Add unique edges for the bubble counting problem. Ignores multiplicity. 
            Reverse creates edges from suffix to prefix to build a reverse graph.
        '''
        def pairwise(iterable):
            "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

        for prefix, suffix in pairwise(k_minus_one_mers):
            self.unique_edge_counts[(prefix, suffix)] += 1
            if (prefix, suffix) in self.unique_edge_counts:
                continue

            if not self.has_node(prefix):
                self.add_node(Node(prefix))
            if not self.has_node(suffix):
                self.add_node(Node(suffix))

            if not reverse:
                self.nodes[prefix].add_edge(suffix)
                self.nodes[suffix].in_degree += 1
            else:
                self.nodes[suffix].add_edge(prefix)
                self.nodes[prefix].in_degree += 1
            self.num_edges += 1

    def dfs(self, cur_node_id, visited, post_visit=None):
        if cur_node_id in visited and visited[cur_node_id]:
            return
        visited[cur_node_id] = True
        for next_node_id in self.nodes[cur_node_id].edges:
            self.dfs(next_node_id, visited, post_visit)
        if post_visit is not None:
            post_visit(self, cur_node_id)

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


class Node:
    def __init__(self, id_):
        # ID contains the prefix/suffix string
        self.id = id_
        # Edges are the string, or the key into the constructor
        self.edges = list()
        self.tours = list()
        self.in_degree = 0
        # Values for tip removal
        self.clock_val = -1
        self.scc_num = -1

    def get_edge(self):
        ''' Gets and removes an edge. '''
        edge = self.edges.pop()
        if not self.edges:
            self.notify_tours_no_available_edges()
        return edge

    @property
    def out_degree(self):
        return len(self.edges)

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


class ReadBreaker:
    ''' Takes input, breaks down into k-1-mers, builds a EulerianPathConstructor as a graph'''

    @staticmethod
    def read_input():
        ''' Returns a list of reads, not yet broken down. '''
        reads = list()
        kmer_len = 20  # 15-20
        for read in sys.stdin:
            reads.append(read.split()[0])
        return (kmer_len, reads)

    @staticmethod
    def _break_read_into_k_minus_one_mers(k, read):
        ''' Breaks into list of prefixes and suffixes. '''
        return [read[i:i+k-1] for i in range(len(read)-(k-2))]

    @staticmethod
    def build_path_constructor(reads, kmer_len):
        ''' Builds the path constructor with coverage considerations
            i.e. Filter out edges with coverage under our set Hamming distance
        '''
        def pairwise(iterable):
            "s -> (s0,s1), (s1,s2), (s2, s3), ... (FROM ITERTOOLS)"
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

        HAMMING_DISTANCE = 5
        path_constructor = EulerianPathConstructor()
        edge_counts = defaultdict(int)
        for cur_read in reads:
            k_minus_one_mers = ReadBreaker._break_read_into_k_minus_one_mers(
                kmer_len, cur_read)

            for prefix, suffix in pairwise(k_minus_one_mers):
                edge_counts[(prefix, suffix)] += 1

        for edge, count in edge_counts.items():
            prefix, suffix = edge[0], edge[1]
            if count > HAMMING_DISTANCE:
                if not path_constructor.has_node(prefix):
                    path_constructor.add_node(Node(prefix))
                if not path_constructor.has_node(suffix):
                    path_constructor.add_node(Node(suffix))

                path_constructor.nodes[prefix].add_edge(suffix)
                path_constructor.nodes[suffix].in_degree += 1
                path_constructor.num_edges += 1

        return path_constructor


def main():
    k, reads = ReadBreaker.read_input()
    path_constructor = ReadBreaker.build_path_constructor(reads, k)
    found, path = path_constructor.build_eulerian_path()
    if found:
        print("".join(map(str, path)))
    else:
        print("NOT FOUND")


if __name__ == "__main__":
    main()
