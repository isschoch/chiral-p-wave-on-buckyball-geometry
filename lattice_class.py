import matplotlib.pyplot as plt
import math

class Lattice:
    def __init__(self):
        self._num_sites = 0
        self._num_edges = 0
        self._sites = []
        self._edges = []

    def add_site(self, coords):
        self._sites.append(coords)
        self._edges.append([])
        self._num_sites = self._num_sites + 1

    def add_edge(self, sites):
        assert sites[0] < len(self._sites) and sites[1] < len(self._sites)
        self._edges[sites[0]].append(sites[1])
        self._edges[sites[1]].append(sites[0])
        self._num_edges = self._num_edges + 1
    
    def remove_edge(self, sites):
        assert sites[0] < len(self._sites) and sites[1] < len(self._sites)
        self._edges[sites[0]].remove(sites[1])
        self._edges[sites[1]].remove(sites[0])
        self._num_edges = self._num_edges - 1

    def plot(self, show_idx_bool=False):
        x_vals = [self._sites[i][0] for i in range(self._num_sites)]
        y_vals = [self._sites[i][1] for i in range(self._num_sites)]

        fig, ax = plt.subplots()

        ax.set_box_aspect(1)

        # plot edges
        for idx, site in enumerate(self._sites):
            if show_idx_bool:
                ax.text(site[0], site[1], "{}".format(idx))
            for neighbour in self._edges[idx]:
                ax.plot([site[0], self._sites[neighbour][0]], [site[1], self._sites[neighbour][1]], c="lightsteelblue",zorder=0)

        # plot sites
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        plt.grid(linestyle=":")
        ax.scatter(x_vals, y_vals, c="red", zorder=1)

        plt.show()

    def get_num_sites(self):
        return self._num_sites
    
    def get_num_edges(self):
        return self._num_edges
    
    def get_sites(self):
        return self._sites

    def get_edges(self):
        return self._edges

class BravaisLattice(Lattice):
    def __init__(self, a1, a2, N, M, BC="open"):
        super().__init__()

        self._lattice_vec1 = a1
        self._lattice_vec2 = a2
        
        self._num_x_repititions = M
        self._num_y_repititions = N

        for i in range(self._num_y_repititions-1, -1, -1):
            for j in range(self._num_x_repititions):
                self.add_site((i*self._lattice_vec2[0] + j*self._lattice_vec1[0], i*self._lattice_vec2[1] + j*self._lattice_vec1[1]))

        assert BC == "open" or BC == "periodic"

        if BC == "open":
            for i in range(self._num_y_repititions-1):
                self.add_edge((i*self._num_x_repititions+self._num_x_repititions-1, (i+1)*self._num_x_repititions+self._num_x_repititions-1))
                for j in range(self._num_x_repititions-1):
                    self.add_edge((i*self._num_x_repititions+j, i*self._num_x_repititions+(j+1)))
                    self.add_edge(((i+1)*self._num_x_repititions+j, i*self._num_x_repititions+j))
            for j in range(self._num_x_repititions-1):
                self.add_edge(((self._num_y_repititions-1)*self._num_x_repititions+j, (self._num_y_repititions-1)*self._num_x_repititions+j+1))
        elif BC == "periodic":
            print("hello there")
            for i in range(self._num_y_repititions):
                for j in range(self._num_x_repititions):
                    self.add_edge((i*self._num_x_repititions+j, i*self._num_x_repititions+(j+1)%self._num_x_repititions))
                    self.add_edge((((i+1)*self._num_x_repititions+j)%(self._num_x_repititions*self._num_y_repititions), i*self._num_x_repititions+j))
        
    def get_reciprocal_sites(self):
        a1 = self._lattice_vec1
        a2 = self._lattice_vec2
        reciprocal_a1 = [+2.0 * math.pi / (a1[0]*a2[1]-a1[1]*a2[0]) * a2[1], -2.0 * math.pi / (a1[0]*a2[1]-a1[1]*a2[0]) * a2[0]]
        reciprocal_a2 = [-2.0 * math.pi / (a1[0]*a2[1]-a1[1]*a2[0]) * a1[1], +2.0 * math.pi / (a1[0]*a2[1]-a1[1]*a2[0]) * a1[0]]

        assert math.isclose((a1[0]*reciprocal_a1[0] + a1[1]*reciprocal_a1[1]) / (2.0*math.pi), 1.0), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose((a2[0]*reciprocal_a2[0] + a2[1]*reciprocal_a2[1]) / (2.0*math.pi), 1.0), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose((a1[0]*reciprocal_a2[0] + a1[1]*reciprocal_a2[1]) / (2.0*math.pi), 0.0), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose((a2[0]*reciprocal_a1[0] + a2[1]*reciprocal_a1[1]) / (2.0*math.pi), 0.0), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"

        reciprocal_sites = []
        for i in range(self._num_y_repititions-1, -1, -1):
            for j in range(self._num_x_repititions):
                reciprocal_sites.append((i*reciprocal_a2[0] + j*reciprocal_a1[0], i*reciprocal_a2[1] + j*reciprocal_a1[1]))
        
        assert len(reciprocal_sites) == len(self._sites), "reciprocal lattice has different number of points than direct lattice"

        return reciprocal_sites

class SquareLatticeAlt(BravaisLattice):
    def __init__(self, N, BC="open"):
        a1 = (1.0, 0.0)
        a2 = (0.0, 1.0)
        super().__init__(a1, a2, N, N, BC=BC)

class SquareLattice(Lattice): 
    def __init__(self, N, BC="open"):
        super().__init__()
        for i in range(N):
            for j in range(N):
                self.add_site((float(j/(N-1)), float((N-1-i)/(N-1))))
        
        assert BC == "open" or BC == "periodic"

        if BC == "open":
            for i in range(N-1):
                self.add_edge((i*N+N-1, (i+1)*N+N-1))
                self.add_edge(((N-1)*N+i, (N-1)*N+i+1))
                for j in range(N-1):
                    self.add_edge((i*N+j, i*N+(j+1)))
                    self.add_edge(((i+1)*N+j, i*N+j))
        elif BC == "periodic":
            for i in range(N):
                for j in range(N):
                    self.add_edge((i*N+j, i*N+(j+1)%N))
                    self.add_edge((((i+1)*N+j)%(N**2), i*N+j))

        
    def get_reciprocal_sites(self):
        return [(2.0 * math.pi * (math.sqrt(self.get_num_sites())-1)**2 * r[0], 2.0 * math.pi * (math.sqrt(self.get_num_sites())-1)**2 * r[1]) for r in self._sites]

class ChainLatticeAlt(BravaisLattice):
    def __init__(self, N, BC="open"):
        a1 = (1.0, 0.0)
        a2 = (0.0, 1.0)
        super().__init__(a1, a2, 1, N, BC=BC)
        
        if BC == "periodic":
            for i in range(self.get_num_sites()):
                self.remove_edge((i, i))

class ChainLattice(Lattice):
    def __init__(self, N, BC="open"):
        super().__init__()
        for i in range(N):
            self.add_site((float(i)/(N-1), 0.0))
        
        assert BC == "open" or "periodic"

        bool = -1
        if BC == "open":
            bool = 1
        elif BC == "periodic":
            bool = 0

        assert bool == 0 or bool == 1

        for i in range(N-bool):
            self.add_edge((i, (i+1)%N))
    
    def get_reciprocal_sites(self):
        return [(2.0 * math.pi * r[0], 2.0 * math.pi * r[1]) for r in self._sites]