from typing import List

import numpy as np
import galois


class ShamirSecretScheme:
    """Incapsulates all functional (creation, encoding, deconding)
      of Shamir secret scheme.
    """
    def __init__(self, m: List[int], n: int, r: int, p: int) -> None:
        """
        Args:
          message (List[int]): message divided into symbols over field Fp
          n (int): number of parties
          r (int): (n - r) - minimal number of parties required to restore
            the secret
          p (int): field order
        """
        if (len(m) % (n - r)) != 0:
            raise ValueError(f"Message length must be divisible by minimal number of parties for decoding: {n - r}")
        
        self.message = m
        self.n_parties = n
        self.degree = n - r

        self.p = p
        self.gf = galois.GF(p)

    def encode(self) -> None:
        """Encoding the message into polynomials,
          distributing shares between parties.
        """
        self._construct_poly()
        self._construct_parties()
    
    def decode(self, indices) -> List[int]:
        """Method that recovers a message from
          requested parties.

        Args:
          indices (List[int]): indices of parties
            to request shares from.
        """

        # choose parties to recover the message
        if len(indices) < self.degree:
            raise ValueError(f"To small amount of parties. At least {self.degree} to recover.")
        
        parties = self._parties[indices[:self.degree]]
        
        message = []
        for idx in range(len(self._polynoms)):
            shares = parties[:, idx]  # (n - r, 2) -> [xi, f(xi)]
            coeffs = galois.lagrange_poly(
                self.gf(shares[:, 0]),
                self.gf(shares[:, 1])
            ).coeffs

            message.append(np.array(coeffs))

        return list(np.concatenate(message))
    
    def _construct_poly(self) -> None:
        # construct polynoms
        n_poly = len(self.message) // self.degree
        polynoms = []
        for idx in range(n_poly):
            coeffs = self.message[idx * self.degree:(idx + 1) * self.degree]
            poly = galois.Poly(coeffs, field=self.gf)

            polynoms.append(poly)
        
        self._polynoms = polynoms
    
    def _construct_parties(self) -> None:
        # generate x ([0; n_parties])
        xs = np.arange(self.n_parties) 

        # distribute (xi, f1(xi), ..., fk(xi)) among parties
        parties = []
        for poly in self._polynoms:
            ys = np.array(poly(xs))

            shares = np.stack([xs, ys], axis=1)
            parties.append(shares)
        
        parties = np.array(parties)  # (n_poly, n_parties, 2)
        self._parties = parties.transpose([1, 0, 2])  # (n_parties, n_poly, 2)
