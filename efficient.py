from typing import List

import numpy as np
import galois


class CommunicationEfficientScheme:
    """Shamir-based communication efficient secret sharing scheme.
      For implementation details refer to: https://arxiv.org/pdf/1505.07515.pdf
    """
    def __init__(
        self, m: List[int], n: int, r: int, p: int, D: int
    ) -> None:
        """
        Args:
          m (List[int]): message (List[int]): message divided into symbols over field Fp
          n (int): number of parties
          r (int): (n - r) - minimal number of parties required to restore
            the secret
          p (int): field order
          D (int): power |D| of set {d1, d2, ..., d|D|}
        """
        # Nessecary conditions
        if D > (n / (n - r)):
            raise ValueError(f"Size of the parties set should be less or equal: {n / (n - r)}, got {D}")
        if len(m) != D * (n - r):
            raise ValueError(f"Message length should be equal to {D * (n - r)}, got {len(m)}")
        
        self.message = m
        self.n_parties = n

        self.Ds = [len(m) // i for i in range(1, D + 1)]

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
        if len(indices) not in self.Ds:
            raise ValueError(f"Number of parties must be one of: {self.Ds}, got: {len(indices)}")
        
        d_idx = [idx for idx, d in enumerate(self.Ds) if len(indices) == d][0]
        parties = self._parties[indices]

        # reconstruct iteratively (start from the smallest polynom based on d_idx)
        ho_message = []
        lo_message = []
        for idx in range(d_idx, 0, -1):
            shares = parties[:, idx].copy()  # (d, 2) -> [xi, f(xi)
            n_ho = self.Ds[idx - 1] - self.Ds[idx]

            xs = self.gf(shares[:, 0])
            ys = self.gf(shares[:, 1])

            if len(ho_message):
                max_deg = self.Ds[idx] - 1
                for i, m_ in enumerate(ho_message):
                    ys = ys - m_ * (xs ** (max_deg - i))    
            
            coeffs = galois.lagrange_poly(
                xs, ys
            ).coeffs
            coeffs = list(np.array(coeffs))

            ho_message += coeffs[:n_ho]
            lo_message += coeffs[n_ho:][::-1]

        # recover final message from biggest polynomial
        message = lo_message + ho_message[::-1]
        message = message[::-1]

        shares = parties[:, 0].copy()
        xs = self.gf(shares[:, 0])
        ys = self.gf(shares[:, 1])

        max_deg = self.Ds[0] - 1
        for i, m_ in enumerate(message):
            ys = ys - m_ * (xs ** (max_deg - i))   

        coeffs = galois.lagrange_poly(
            xs, ys
        ).coeffs
        coeffs = list(np.array(coeffs))

        message = coeffs[::-1] + message[::-1]    

        return message
    
    def _construct_poly(self) -> None:
        # construct polynoms (from big to small)

        n_ho = []
        for i in range(1, len(self.Ds)):
            n_ho.append(min(self.Ds[i - 1] - self.Ds[-1], self.Ds[i]))

        n_lo = []
        for i in range(1, len(self.Ds)):
            n_lo.append(self.Ds[i] - n_ho[i - 1])

        polynoms = []
        ilo_idx = -n_ho[0]

        for n_iho, n_ilo in zip(n_ho, n_lo):
            lho = self.message[ilo_idx - n_ilo:ilo_idx]
            iho = self.message[-n_iho:]
            
            coeffs = (lho + iho)[::-1]  # reverse so last symbols of high order
            poly = galois.Poly(coeffs, field=self.gf)
            polynoms.append(poly)

            ilo_idx -= n_ilo
        
        # construct first (biggest) polynom
        coeffs = self.message[::-1]
        poly = galois.Poly(coeffs, field=self.gf)
        polynoms.insert(0, poly)

        self._polynoms = polynoms
    
    def _construct_parties(self) -> None:
        # generate x ([0; n_parties])
        xs = np.arange(self.n_parties)

        # distribute (xi, pd_1(xi), ..., pd_D(xi)) among parties
        parties = []
        for poly in self._polynoms:
            ys = np.array(poly(xs))

            shares = np.stack([xs, ys], axis=1)
            parties.append(shares)
        
        parties = np.array(parties)  # (n_poly, n_parties, 2)
        self._parties = parties.transpose([1, 0, 2])  # (n_parties, n_poly, 2)
