import numpy as np


class HamiltonianConstructor:
    def __init__(self, mu, tx, ty, delta):

        self.H_tunnel_x = 0.5 * np.diag([tx, tx, -tx, -tx])
        self.H_tunnel_y = 0.5 * np.diag([ty, ty, -ty, -ty])

        H_sc_x = 0.5 * np.array(
            [
                [0.0, 0.0, 0.0, complex(0.0, 1.0) * delta],
                [0.0, 0.0, complex(0.0, 1.0) * delta, 0.0],
                [0.0, complex(0.0, 1.0) * delta.conjugate(), 0.0, 0.0],
                [complex(0.0, 1.0) * delta.conjugate(), 0.0, 0.0, 0.0],
            ],
            dtype=complex,
        )

        H_sc_y = 0.5 * np.array(
            [
                [0.0, 0.0, 0.0, -delta],
                [0.0, 0.0, -delta, 0.0],
                [0.0, delta.conjugate(), 0.0, 0.0],
                [delta.conjugate(), 0.0, 0.0, 0.0],
            ],
            dtype=complex,
        )

        H_right = self.H_tunnel_x + H_sc_x
        H_left = self.H_tunnel_x - H_sc_x

        H_up = self.H_tunnel_y + H_sc_y
        H_down = self.H_tunnel_y - H_sc_y

        self.hop_hamiltonians = [H_up, H_right, H_down, H_left]
        self.site_hamiltonian = np.diag([-mu, -mu, mu, mu])
        self.dim_H_BdG = 4

    def get_hop_hamiltonian(self, local_edge_idx):
        return self.hop_hamiltonians[local_edge_idx]

    def get_block_indices(self, num_sites):
        return [
            slice(block_idx * self.dim_H_BdG, (block_idx + 1) * self.dim_H_BdG)
            for block_idx in range(num_sites)
        ]

    def construct_direct_lattice_hamiltonian(self, lattice):
        num_sites = lattice.get_num_sites()
        lattice_edges = lattice.get_edges()

        H_direct_lattice = np.zeros(
            shape=(num_sites * self.dim_H_BdG, num_sites * self.dim_H_BdG),
            dtype=complex,
        )

        block_indices = self.get_block_indices(num_sites)

        for (site_idx, site) in enumerate(lattice.get_sites()):
            H_direct_lattice[
                block_indices[site_idx], block_indices[site_idx]
            ] = self.site_hamiltonian

            for neighbour_idx in lattice_edges[site_idx]:
                H_direct_lattice[
                    block_indices[site_idx], block_indices[neighbour_idx]
                ] = self.hop_hamiltonians[
                    lattice.get_local_edge_index(site_idx, neighbour_idx)
                ]

        assert self._hermitian_check(
            H_direct_lattice
        ), "Hermitian Hamiltonian check assertion failed."

        return H_direct_lattice

    def add_tunneling_phase(self, H_direct_lattice, edge, phase, num_sites):
        block_indices = self.get_block_indices(num_sites)
        np.fill_diagonal(
            H_direct_lattice[block_indices[edge[0]], block_indices[edge[1]]],
            phase
            * np.diag(
                H_direct_lattice[
                    block_indices[edge[0]],
                    block_indices[edge[1]],
                ]
            ),
        )
        np.fill_diagonal(
            H_direct_lattice[block_indices[edge[1]], block_indices[edge[0]]],
            phase.conjugate()
            * np.diag(
                H_direct_lattice[
                    block_indices[edge[1]],
                    block_indices[edge[0]],
                ]
            ),
        )

        assert self._hermitian_check(
            H_direct_lattice
        ), "Hermitian Hamtilonian check assertion failed."

    def _hermitian_check(self, A, rtol=1e-05, atol=1e-08):
        return np.allclose(A, A.conj().T, rtol=rtol, atol=atol)
