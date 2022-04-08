import site
import matplotlib.pyplot as plt
import math


class Lattice:
    def __init__(self):
        self._num_sites = 0
        self._num_bonds = 0
        self._sites = []
        self._bonds = []
        self._local_bond_idx = []

    def add_site(self, coords):
        self._sites.append(coords)
        self._bonds.append([])
        self._local_bond_idx.append([])
        self._num_sites = self._num_sites + 1

    def add_sites(self, coords_list):
        for coords in coords_list:
            self.add_site(coords)

    def move_site(self, site_idx, coords):
        assert site_idx < self._num_sites
        self._sites[site_idx] = coords

    def remove_site(self, site_idx):
        # remove all incoming and outgoing bonds
        assert site_idx < len(self._sites), "Site exists assertion failed"
        bond_cnt = 0
        for (bond_idx, bond_lst) in enumerate(self._bonds):
            for neigh in bond_lst:
                if neigh == site_idx:
                    self.remove_bond((site_idx, bond_idx))
                    bond_cnt = bond_cnt + 1

        # remove site from self._sites and self._bonds
        self._sites.pop(site_idx)
        self._bonds.pop(site_idx)
        self._local_bond_idx.pop(site_idx)

        # shift inidices
        for (bond_idx, bond_lst) in enumerate(self._bonds):
            for (neigh_idx, neigh) in enumerate(bond_lst):
                if neigh > site_idx:
                    self._bonds[bond_idx][neigh_idx] = (
                        self._bonds[bond_idx][neigh_idx] - 1
                    )

        # adjust node and bond count
        self._num_sites = len(self._sites)

        assert (
            self.get_num_sites() >= 0 and self.get_num_bonds() >= 0
        ), "Postivie number of sites and bonds assertion failed"

    def remove_sites(self, site_indices):
        site_indices.sort(reverse=True)
        for site_idx in site_indices:
            self.remove_site(site_idx)

    def add_bond(self, sites, local_bond_indicies=(None, None)):
        assert sites[0] < len(self._sites) and sites[1] < len(
            self._sites
        ), "Both sites exist assertion failed while trying to add bond."
        self._bonds[sites[0]].append(sites[1])
        self._bonds[sites[1]].append(sites[0])

        self._local_bond_idx[sites[0]].append(local_bond_indicies[0])
        self._local_bond_idx[sites[1]].append(local_bond_indicies[1])

        self._num_bonds = self._num_bonds + 1

    def add_bonds(self, bonds, local_bond_indices_list=None):
        if local_bond_indices_list == None:
            local_bond_indices_list = [(None, None) for bond in bonds]

        for bond_lcl_idx in zip(bonds, local_bond_indices_list):
            self.add_bond(bond_lcl_idx[0], bond_lcl_idx[1])

    def remove_bond(self, sites):
        assert sites[0] < len(self._sites) and sites[1] < len(self._sites)

        site_0_list_idx = self._bonds[sites[0]].index(sites[1])
        site_1_list_idx = self._bonds[sites[1]].index(sites[0])

        self._bonds[sites[0]].pop(site_0_list_idx)
        self._bonds[sites[1]].pop(site_1_list_idx)

        self._local_bond_idx[sites[0]].pop(site_0_list_idx)
        self._local_bond_idx[sites[1]].pop(site_1_list_idx)

        self._num_bonds = self._num_bonds - 1

    def remove_bonds(self, bonds):
        for bond in bonds:
            self.remove_bond(bond)

    def glue_bond(self, bonds1, bonds2):
        assert len(bonds1) == len(
            bonds2
        ), "Assertion that bond lists have same length failed"
        sites1 = []

    def plot(self, show_idx_bool=False, **kwargs):
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

        # plot bonds
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
            for neighbour in self._bonds[idx]:
                ax.plot(
                    [site[0], self._sites[neighbour][0]],
                    [site[1], self._sites[neighbour][1]],
                    c="black",
                    alpha=0.2,
                    zorder=5,
                )

        # plot sites
        plt.xticks(list(dict.fromkeys(x_vals)))
        plt.yticks(list(dict.fromkeys(y_vals)))

        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        plt.grid(linestyle=":", zorder=-10)

        return ax.scatter(
            x_vals,
            y_vals,
            zorder=10,
            s=site_size,
            **kwargs,
        )

    def get_num_sites(self):
        return self._num_sites

    def get_num_bonds(self):
        return self._num_bonds

    def get_sites(self):
        return self._sites

    def get_bonds(self):
        return self._bonds

    def get_local_bond_indices(self):
        return self._local_bond_idx

    def get_local_bond_index(self, site_idx, neighbour_idx):
        assert (
            neighbour_idx in self._bonds[site_idx]
        ), "Site {} and {} are not adjacent".format(site_idx, neighbour_idx)
        return self._local_bond_idx[site_idx][
            self._bonds[site_idx].index(neighbour_idx)
        ]

    def change_local_bond_index(self, sites, new_local_idx):
        assert (sites[1] in self._bonds[sites[0]]) and (
            sites[0] in self._bonds[sites[1]]
        ), "Site {} and {} are not adjacent".format(sites[0], sites[1])
        self._local_bond_idx[sites[0]][
            self._bonds[sites[0]].index(sites[1])
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
                self.add_bond(
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
                    self.add_bond(
                        (
                            i * self._num_vec1_repititions + j,
                            i * self._num_vec1_repititions + (j + 1),
                        ),
                        (
                            1,
                            3,
                        ),
                    )
                    self.add_bond(
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
                self.add_bond(
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
            for i in range(self._num_vec2_repititions):
                for j in range(self._num_vec1_repititions):
                    self.add_bond(
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
                    self.add_bond(
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
                self.remove_bond((i, i))


class CrystalLattice(Lattice):
    def __init__(
        self,
        a1,
        a2,
        N,
        M,
        unit_cell_sites,
        interrior_unit_cell_bonds=[],
        exterior_unit_cell_bonds=[],
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
        for bond in interrior_unit_cell_bonds:
            for cell_idx in range(self.get_num_sites() // num_unit_cell_sites):
                self.add_bond(
                    (
                        bond[0] + cell_idx * num_unit_cell_sites,
                        bond[1] + cell_idx * num_unit_cell_sites,
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
                for cell_neigh_idx in exterior_unit_cell_bonds:
                    for bond in exterior_unit_cell_bonds[cell_neigh_idx]:
                        self.add_bond(
                            (
                                (cell_idx * num_unit_cell_sites + bond[0]) % num_sites,
                                (
                                    cell_idx * num_unit_cell_sites
                                    + bond[1]
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
            interrior_unit_cell_bonds=[(0, 1)],
            exterior_unit_cell_bonds={7: [(0, 1)], 8: [(0, 1)]},
            BC=BC,
        )

        num_cells = self.get_num_cells()
        num_unit_cell_sites = self.get_num_unit_cell_sites()
        for i in range(self._num_vec2_repititions - 1):
            self.add_bond(
                (
                    i * self._num_vec1_repititions * num_unit_cell_sites
                    + self._num_vec1_repititions * num_unit_cell_sites
                    - 2,
                    (i + 1) * self._num_vec1_repititions * num_unit_cell_sites
                    + self._num_vec1_repititions * num_unit_cell_sites
                    - 1,
                )
            )
