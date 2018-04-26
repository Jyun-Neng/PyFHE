from NTT import NTT
from sympy.ntheory.modular import crt
import numpy as np


class CRTPoly:
    """
    Data structure:
        crt_poly, prime_set 
    """

    def __init__(self, poly=None, primes=None, fft=True, crt=False, N=None):
        self.do_fft = fft
        if crt:
            self.N = N
            self.initial_w_crt(poly, primes)
        else:
            self.N = len(poly)
            self.initial_wo_crt(poly, primes)

    def initial_wo_crt(self, poly, primes):
        self.prime_set = list(primes)
        self.crt_poly = self.crtPoly(poly)

    def initial_w_crt(self, poly, primes):
        self.prime_set = list(primes)
        self.crt_poly = poly

    def __name__(self):
        return "CRTPoly"

    def __str__(self):
        return str(self.crt_poly)

    def __add__(self, other):

        assert self.do_fft == other.do_fft
        assert self.N == other.N

        add_result = []
        same_size = False

        # check if two poly prime sets are equal
        if len(self.crt_poly) > len(other.crt_poly):
            small_obj = other
            large_obj = self
        elif len(self.crt_poly) < len(other.crt_poly):
            small_obj = self
            large_obj = other
        else:
            small_obj = self
            same_size = True

        for i, poly in enumerate(small_obj.crt_poly):
            if self.do_fft:
                add_result.append(poly + other.crt_poly[i])

            else:
                _result = np.asarray(poly) + np.asarray(other.crt_poly[i])
                result = np.fmod(_result, self.prime_set[i])
                add_result.append(result.tolist())

        if not same_size:
            prime_set = list(small_obj.prime_set)
            N = small_obj.N
        else:
            prime_set = list(self.prime_set)
            N = self.N

        return CRTPoly(add_result, prime_set, self.do_fft, True, N=N)

    def __sub__(self, other):

        assert self.do_fft == other.do_fft
        assert self.N == other.N
        assert len(self.prime_set) >= len(other.prime_set)

        sub_result = []

        for i, poly in enumerate(other.crt_poly):
            if self.do_fft:
                sub_result.append(self.crt_poly[i] - poly)
            else:
                _result = np.asarray(self.crt_poly[i]) - np.asarray(poly)
                result = np.fmod(_result, self.prime_set[i])
                sub_result.append(result.tolist())

        prime_set = list(other.prime_set)
        return CRTPoly(sub_result, prime_set, self.do_fft, True, N=other.N)

    def __mul__(self, other):

        mul_result = []

        if type(other).__name__ == 'int':
            for poly in self.crt_poly:
                if self.do_fft:
                    result = poly * other
                else:
                    result = (np.asarray(poly) * other).tolist()
                mul_result.append(result)
            prime_set = self.prime_set
        else:
            assert self.do_fft and other.do_fft

            for i, poly in enumerate(self.crt_poly):
                result = poly * other.crt_poly[i]
                mul_result.append(result)
            if len(self.prime_set) <= len(other.prime_set):
                prime_set = self.prime_set
            else:
                prime_set = self.prime_set

        return CRTPoly(mul_result, prime_set, self.do_fft, True, N=self.N)

    def crtPoly(self, poly, primes=None):
        """
        Transform poly to CRT form, then transform each CRT poly to frequency domain,
        if do_fft is true. 
        """
        if primes is None:
            primes = self.prime_set
        result = []
        for prime in primes:
            crt_poly = np.remainder(poly, prime).tolist()
            if self.do_fft:
                fft_crt_poly = NTT(crt_poly, prime, self.N)
                result.append(fft_crt_poly)
            else:
                result.append(crt_poly)
        return result

    def toPoly(self):
        """
        CRT poly. 
        """
        if self.do_fft:
            polys = []
            for fft_poly in self.crt_poly:
                polys.append(fft_poly.intt())

        else:
            polys = self.crt_poly

        poly = []
        residue_array = np.asarray(polys).T
        for residues in residue_array:
            coeff, M = crt(self.prime_set, residues)
            poly.append(coeff)

        return poly
