import site
import matplotlib.pyplot as plt
import math


class Lattice:
    def __init__(self):
        self._num_sites = 0
        self._num_edges = 0
        self._sites = []
        self._edges = []
        self._local_edge_idx = []

    def add_site(self, coords):
        self._sites.append(coords)
        self._edges.append([])
        self._local_edge_idx.append([])
        self._num_sites = self._num_sites + 1

    def move_site(self, site_idx, coords):
        assert site_idx < self._num_sites
        self._sites[site_idx] = coords

    def remove_site(self, site_idx):
        # remove all incoming and outgoing edges
        assert site_idx < len(self._sites), "Site exists assertion failed"
        edge_cnt = 0
        for (edge_idx, edge_lst) in enumerate(self._edges):
            for neigh in edge_lst:
                if neigh == site_idx:
                    self.remove_edge((site_idx, edge_idx))
                    edge_cnt = edge_cnt + 1

        # remove site from self._sites and self._edges
        self._sites.pop(site_idx)
        self._edges.pop(site_idx)
        self._local_edge_idx.pop(site_idx)

        # shift inidices
        for (edge_idx, edge_lst) in enumerate(self._edges):
            for (neigh_idx, neigh) in enumerate(edge_lst):
                if neigh > site_idx:
                    self._edges[edge_idx][neigh_idx] = (
                        self._edges[edge_idx][neigh_idx] - 1
                    )

        # adjust node and edge count
        self._num_sites = len(self._sites)

        assert (
            self.get_num_sites() >= 0 and self.get_num_edges() >= 0
        ), "Postivie number of sites and edges assertion failed"

    def add_edge(self, sites, local_edge_indicies=(None, None)):
        assert sites[0] < len(self._sites) and sites[1] < len(
            self._sites
        ), "Both sites exist assertion failed while trying to add edge."
        self._edges[sites[0]].append(sites[1])
        self._edges[sites[1]].append(sites[0])

        self._local_edge_idx[sites[0]].append(local_edge_indicies[0])
        self._local_edge_idx[sites[1]].append(local_edge_indicies[1])

        self._num_edges = self._num_edges + 1

    def remove_edge(self, sites):
        assert sites[0] < len(self._sites) and sites[1] < len(self._sites)

        site_0_list_idx = self._edges[sites[0]].index(sites[1])
        site_1_list_idx = self._edges[sites[1]].index(sites[0])

        self._edges[sites[0]].pop(site_0_list_idx)
        self._edges[sites[1]].pop(site_1_list_idx)

        self._local_edge_idx[sites[0]].pop(site_0_list_idx)
        self._local_edge_idx[sites[1]].pop(site_1_list_idx)

        self._num_edges = self._num_edges - 1

    def plot(self, show_idx_bool=False):
        x_vals = [self._sites[i][0] for i in range(self.get_num_sites())]
        y_vals = [self._sites[i][1] for i in range(self.get_num_sites())]

        fig, ax = plt.subplots()

        x_width = max(x_vals) - min(x_vals)
        y_width = max(y_vals) - min(y_vals)
        min_ratio_param = 0.1
        y_x_ratio = max(y_width, min_ratio_param * x_width) / max(
            x_width, min_ratio_param * y_width
        )
        ax.set_box_aspect(y_x_ratio)

        site_size = max(500 / self.get_num_sites(), 10)

        # plot edges
        for idx, site in enumerate(self._sites):
            if show_idx_bool:
                ax.text(
                    site[0],
                    site[1],
                    "{}".format(idx),
                    zorder=100,
                    c="teal",
                    size=8 * (1 / (1 + math.exp(-math.sqrt(site_size))) ** 2),
                    horizontalalignment="right",
                    verticalalignment="bottom",
                )
            for neighbour in self._edges[idx]:
                ax.plot(
                    [site[0], self._sites[neighbour][0]],
                    [site[1], self._sites[neighbour][1]],
                    c="lightsteelblue",
                    zorder=5,
                )

        # plot sites
        plt.xticks(list(dict.fromkeys(x_vals)))
        plt.yticks(list(dict.fromkeys(y_vals)))

        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        plt.grid(linestyle=":", zorder=-10)

        ax.scatter(
            x_vals,
            y_vals,
            c="red",
            zorder=10,
            s=site_size,
        )

        plt.show()

    def get_num_sites(self):
        return self._num_sites

    def get_num_edges(self):
        return self._num_edges

    def get_sites(self):
        return self._sites

    def get_edges(self):
        return self._edges

    def get_local_edge_indices(self):
        return self._local_edge_idx

    def get_local_edge_index(self, site_idx, neighbour_idx):
        assert (
            neighbour_idx in self._edges[site_idx]
        ), "Site {} and {} are not adjacent".format(site_idx, neighbour_idx)
        return self._local_edge_idx[site_idx][
            self._edges[site_idx].index(neighbour_idx)
        ]

    def change_local_edge_index(self, sites, new_local_idx):
        assert (sites[1] in self._edges[sites[0]]) and (
            sites[0] in self._edges[sites[1]]
        ), "Site {} and {} are not adjacent".format(sites[0], sites[1])
        self._local_edge_idx[sites[0]][
            self._edges[sites[0]].index(sites[1])
        ] = new_local_idx


class BravaisLattice(Lattice):
    def __init__(self, a1, a2, N, M, BC="open"):
        super().__init__()

        self._lattice_vec1 = a1
        self._lattice_vec2 = a2

        self._num_vec1_repititions = M
        self._num_vec2_repititions = N

        for i in range(self._num_vec2_repititions - 1, -1, -1):
            for j in range(self._num_vec1_repititions):
                self.add_site(
                    (
                        i * self._lattice_vec2[0] + j * self._lattice_vec1[0],
                        i * self._lattice_vec2[1] + j * self._lattice_vec1[1],
                    )
                )

        assert BC == "open" or BC == "periodic"

        if BC == "open":
            for i in range(self._num_vec2_repititions - 1):
                self.add_edge(
                    (
                        i * self._num_vec1_repititions + self._num_vec1_repititions - 1,
                        (i + 1) * self._num_vec1_repititions
                        + self._num_vec1_repititions
                        - 1,
                    ),
                    (
                        2,
                        0,
                    ),
                )
                for j in range(self._num_vec1_repititions - 1):
                    self.add_edge(
                        (
                            i * self._num_vec1_repititions + j,
                            i * self._num_vec1_repititions + (j + 1),
                        ),
                        (
                            1,
                            3,
                        ),
                    )
                    self.add_edge(
                        (
                            i * self._num_vec1_repititions + j,
                            (i + 1) * self._num_vec1_repititions + j,
                        ),
                        (
                            2,
                            0,
                        ),
                    )
            for j in range(self._num_vec1_repititions - 1):
                self.add_edge(
                    (
                        (self._num_vec2_repititions - 1) * self._num_vec1_repititions
                        + j,
                        (self._num_vec2_repititions - 1) * self._num_vec1_repititions
                        + j
                        + 1,
                    ),
                    (
                        1,
                        3,
                    ),
                )
        elif BC == "periodic":
            print("hello there")
            for i in range(self._num_vec2_repititions):
                for j in range(self._num_vec1_repititions):
                    self.add_edge(
                        (
                            i * self._num_vec1_repititions + j,
                            i * self._num_vec1_repititions
                            + (j + 1) % self._num_vec1_repititions,
                        ),
                        (
                            1,
                            3,
                        ),
                    )
                    self.add_edge(
                        (
                            ((i + 1) * self._num_vec1_repititions + j)
                            % (self._num_vec1_repititions * self._num_vec2_repititions),
                            i * self._num_vec1_repititions + j,
                        ),
                        (
                            0,
                            2,
                        ),
                    )

    def get_reciprocal_sites(self):
        a1 = self._lattice_vec1
        a2 = self._lattice_vec2
        reciprocal_a1 = [
            +2.0 * math.pi / (a1[0] * a2[1] - a1[1] * a2[0]) * a2[1],
            -2.0 * math.pi / (a1[0] * a2[1] - a1[1] * a2[0]) * a2[0],
        ]
        reciprocal_a2 = [
            -2.0 * math.pi / (a1[0] * a2[1] - a1[1] * a2[0]) * a1[1],
            +2.0 * math.pi / (a1[0] * a2[1] - a1[1] * a2[0]) * a1[0],
        ]

        assert math.isclose(
            (a1[0] * reciprocal_a1[0] + a1[1] * reciprocal_a1[1]) / (2.0 * math.pi), 1.0
        ), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose(
            (a2[0] * reciprocal_a2[0] + a2[1] * reciprocal_a2[1]) / (2.0 * math.pi), 1.0
        ), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose(
            (a1[0] * reciprocal_a2[0] + a1[1] * reciprocal_a2[1]) / (2.0 * math.pi), 0.0
        ), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"
        assert math.isclose(
            (a2[0] * reciprocal_a1[0] + a2[1] * reciprocal_a1[1]) / (2.0 * math.pi), 0.0
        ), r"assertion $a_i \cdot b_j = 2*\pi\delta_{i,j}$ failed"

        reciprocal_sites = []
        for i in range(self._num_vec2_repititions - 1, -1, -1):
            for j in range(self._num_vec1_repititions):
                reciprocal_sites.append(
                    (
                        i * reciprocal_a2[0] + j * reciprocal_a1[0],
                        i * reciprocal_a2[1] + j * reciprocal_a1[1],
                    )
                )

        assert len(reciprocal_sites) == len(
            self._sites
        ), "reciprocal lattice has different number of points than direct lattice"

        return reciprocal_sites


class SquareLattice(BravaisLattice):
    def __init__(self, N, BC="open"):
        a1 = (1.0, 0.0)
        a2 = (0.0, 1.0)
        super().__init__(a1, a2, N, N, BC=BC)


class ChainLattice(BravaisLattice):
    def __init__(self, N, BC="open"):
        a1 = (1.0, 0.0)
        a2 = (0.0, 1.0)
        super().__init__(a1, a2, 1, N, BC=BC)

        if BC == "periodic":
            for i in range(self.get_num_sites()):
                self.remove_edge((i, i))


class CrystalLattice(Lattice):
    def __init__(
        self,
        a1,
        a2,
        N,
        M,
        unit_cell_sites,
        interrior_unit_cell_edges=[],
        exterior_unit_cell_edges=[],
        BC="open",
    ):
        super().__init__()

        self._lattice_vec1 = a1
        self._lattice_vec2 = a2

        self._num_vec1_repititions = M
        self._num_vec2_repititions = N

        self._unit_cell_sites = unit_cell_sites
        self._num_unit_cell_sites = len(unit_cell_sites)

        for i in range(self._num_vec2_repititions - 1, -1, -1):
            for j in range(self._num_vec1_repititions):
                for b in unit_cell_sites:
                    self.add_site(
                        (
                            b[0]
                            + i * self._lattice_vec2[0]
                            + j * self._lattice_vec1[0],
                            b[1]
                            + i * self._lattice_vec2[1]
                            + j * self._lattice_vec1[1],
                        )
                    )

        num_unit_cell_sites = self.get_num_unit_cell_sites()
        num_sites = self.get_num_sites()
        for edge in interrior_unit_cell_edges:
            for cell_idx in range(self.get_num_sites() // num_unit_cell_sites):
                self.add_edge(
                    (
                        edge[0] + cell_idx * num_unit_cell_sites,
                        edge[1] + cell_idx * num_unit_cell_sites,
                    )
                )

        cell_neigh_indices = [
            ((i - 1) * self._num_vec1_repititions + (j - 1)) * num_unit_cell_sites
            for i in range(3)
            for j in range(3)
        ]

        for i in range(self._num_vec2_repititions - 1):
            for j in range(self._num_vec1_repititions - 1):
                cell_idx = i * self._num_vec1_repititions + j
                for cell_neigh_idx in exterior_unit_cell_edges:
                    for edge in exterior_unit_cell_edges[cell_neigh_idx]:
                        self.add_edge(
                            (
                                (cell_idx * num_unit_cell_sites + edge[0]) % num_sites,
                                (
                                    cell_idx * num_unit_cell_sites
                                    + edge[1]
                                    + cell_neigh_indices[cell_neigh_idx]
                                )
                                % num_sites,
                            )
                        )

    def get_num_unit_cell_sites(self):
        return self._num_unit_cell_sites

    def get_num_cells(self):
        return self.get_num_sites() // self.get_num_unit_cell_sites()

    def get_local_site_index(self, global_idx):
        return global_idx % self._num_unit_cell_sites


class HoneycombLattice(CrystalLattice):
    def __init__(self, N, M, BC="open"):
        a1 = (1.0, 0.0)
        a2 = (1.0 / 2.0, math.sqrt(3.0) / 2.0)
        unit_cell_sites = [(0.0, 0.25), (0.0, 0.75)]
        super().__init__(
            a1,
            a2,
            N,
            M,
            unit_cell_sites,
            interrior_unit_cell_edges=[(0, 1)],
            exterior_unit_cell_edges={7: [(0, 1)], 8: [(0, 1)]},
            BC=BC,
        )

        # self._lattice_vec1 = a1
        # self._lattice_vec2 = a2

        # self._num_vec1_repititions = M
        # self._num_vec2_repititions = N

        # for i in range(self._num_vec2_repititions - 1, -1, -1):
        #     for j in range(self._num_vec1_repititions):
        #         self.add_site(
        #             (
        #                 i * self._lattice_vec2[0] + j * self._lattice_vec1[0],
        #                 i * self._lattice_vec2[1] + j * self._lattice_vec1[1],
        #             )
        #         )
        #         self.add_site(
        #             (
        #                 b[0] + i * self._lattice_vec2[0] + j * self._lattice_vec1[0],
        #                 b[1] + i * self._lattice_vec2[1] + j * self._lattice_vec1[1],
        #             )
        #         )
        #         self.add_edge((self.get_num_sites()-2, self.get_num_sites()-1))

        # # for i in range(1, self.get_num_sites(), 2):
        # #     self.add_edge((i, (i+2*self._num_vec1_repititions-1)%self.get_num_sites()))
        # #     self.add_edge((i, (i+2*self._num_vec1_repititions+1)%self.get_num_sites()))

        # for i in range(self._num_vec2_repititions-1):
        #     for j in range(self._num_vec1_repititions-1):
        #         self.add_edge((2*(j+i*self._num_vec1_repititions)+1, 2*(j + self._num_vec1_repititions+i*self._num_vec1_repititions)))
        #         self.add_edge((2*(j+i*self._num_vec1_repititions)+1, 2*(j + 1 + self._num_vec1_repititions+i*self._num_vec1_repititions)))
