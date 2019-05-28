import sys


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
