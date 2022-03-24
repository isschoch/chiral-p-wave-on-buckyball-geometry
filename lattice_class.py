import matplotlib.pyplot as plt

class Lattice:
    def __init__(self):
        self._num_nodes = 0
        self._num_edges = 0
        self._nodes = []
        self._edges = []

    def add_node(self, coords):
        self._nodes.append(coords)
        self._edges.append([])
        self._num_nodes = self._num_nodes + 1

    def add_edge(self, nodes):
        assert nodes[0] < len(self._nodes) and nodes[1] < len(self._nodes)
        self._edges[nodes[0]].append(nodes[1])
        self._edges[nodes[1]].append(nodes[0])
        self._num_edges = self._num_edges + 1

    def plot(self):
        x_vals = [self._nodes[i][0] for i in range(self._num_nodes)]
        y_vals = [self._nodes[i][1] for i in range(self._num_nodes)]

        fig, ax = plt.subplots()

        ax.set_box_aspect(1)

        # plot edges
        for idx, node in enumerate(self._nodes):
            for neighbour in self._edges[idx]:
                ax.plot([node[0], self._nodes[neighbour][0]], [node[1], self._nodes[neighbour][1]], c="lightsteelblue",zorder=0)

        # plot nodes
        ax.scatter(x_vals, y_vals, c="red", zorder=1)

        plt.show()

    def get_num_nodes(self):
        return self._num_nodes
    
    def get_num_edges(self):
        return self._num_edges

class SquareLattice(Lattice): 
    def __init__(self, N):
        super().__init__()
        for i in range(N):
            for j in range(N):
                self.add_node((float(j/(N-1)), float((N-1-i)/(N-1))))
        
        for i in range(N):
            for j in range(N-1):
                self.add_edge((i*N+j, i*N+(j+1)%N))
                self.add_edge((((i+1)*N+j)%(N**2), i*N+j))

class ChainLattice(Lattice):
    def __init__(self, N):
        super().__init__()
        for i in range(N):
            self.add_node((float(i)/(N-1), 0.0))
        
        for i in range(N):
            self.add_edge((i, (i+1)%N))