from FHE import FHE
from CRTPoly import CRTPoly
import numpy as np


class Ctxt():
    """
    """

    def __init__(self, d, stdev, primes, P, L, cur_level=0, c=None):
        self.prime_set = list(primes)
        self.prime_set.sort(reverse=True)
        del self.prime_set[:cur_level]
        self.L = cur_level
        self.P = P
        self.f = FHE(d, stdev, primes, P, L, cur_level=cur_level)
        self.c = c

    def __name__(self):
        return 'Ctxt'

    def __add__(self, other):
        """
        Homomorphic addition.
        Return a new Ctxt object with its ciphertext is c + c' = (c0 + c'0, c1 + c'1)
        """

        assert other.__name__() == 'Ctxt'

        while self.L != other.L:
            if self.L < other.L:
                self.scaleDown()
            else:
                other.scaleDown()

        modulus = 1
        for prime in self.prime_set:
            modulus *= prime

        result0 = (np.asarray(self.c[0]) + np.asarray(other.c[0])) % modulus
        result1 = (np.asarray(self.c[1]) + np.asarray(other.c[1])) % modulus

        result = []
        result.append(result0.tolist())
        result.append(result1.tolist())

        return Ctxt(self.f.d, self.f.stdev, self.f.prime_set, self.P, self.f.L, self.L, result)

    def __mul__(self, other):
        """
        Homomorphic multiplication.
        Return a new Ctxt object with its ciphertext is c * c' = (c0 * c'0, c0 * c'1 + c'0 * c1, c1 * c'1)
        """
        assert other.__name__() == 'Ctxt'

        while self.L != other.L:
            if self.L < other.L:
                self.scaleDown()
            else:
                other.scaleDown()

        crt_c10 = CRTPoly(self.c[0], self.prime_set)
        crt_c11 = CRTPoly(self.c[1], self.prime_set)
        crt_c20 = CRTPoly(other.c[0], other.prime_set)
        crt_c21 = CRTPoly(other.c[1], other.prime_set)

        crt_result0 = crt_c10 * crt_c20
        crt_result1 = crt_c10 * crt_c21 + crt_c11 * crt_c20
        crt_result2 = crt_c11 * crt_c21

        result = []
        result.append(crt_result0.toPoly())
        result.append(crt_result1.toPoly())
        result.append(crt_result2.toPoly())

        return Ctxt(self.f.d, self.f.stdev, self.f.prime_set, self.P, self.f.L, self.L, result)

    def enc(self, m, pk):
        """
        FHE encryption. 
        """
        self.c = self.f.homoEnc(m, pk)
        return self.c

    def dec(self, sk):
        """
        FHE decryption. 
        """
        m = self.f.homoDec(self.c, sk)
        return m

    def scaleDown(self):
        """
        Modulus switching down to next level. 
        """
        self.c = self.f.modSwitch(self.c, self.L)
        self.L += 1
        del self.prime_set[0]

    def relinearize(self, switch_key):
        """
        Do key switching and modulus switching. 
        """
        self.c = self.f.keySwitch(self.c, switch_key)
        self.scaleDown()
