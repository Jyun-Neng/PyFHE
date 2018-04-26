############################################################
#                                                          #
# File Name: numTh.py                                      #
#                                                          #
# Author: Jyun-Neng Ji (jyunnengji@gmail.com)              #
#                                                          #
# Creation Date: 2018/02/15                                #
#                                                          #
# Last Modified: 2018/04/12                                #
#                                                          #
# Description: Number theory library.                      #
#                                                          #
############################################################
import math
import random
import numpy as np
from sympy import isprime, nextprime
from sympy.ntheory.residue_ntheory import nthroot_mod


def findPrimes(prime_bit, N, num):
    """
    Give prime bits and number of primes,
    then generate primes, which congruence to 1 modular 2N.

    Examples
    ========

    >>> from numTh import findPrimes
    >>> findPrimes(12, 4, 5)
    [2081, 2089, 2113, 2129, 2137], 60
    """
    primes = []
    total_bits = 0
    prime = pow(2, prime_bit - 1)
    while len(primes) != num:
        prime = nextprime(prime)
        if prime % (2 * N) == 1:
            primes.append(prime)
            total_bits += prime.bit_length()

    return primes, total_bits


def findPrimitiveNthRoot(M, N):
    """
    Generate the smallest primitive Nth root of unity (expect 1).
    Find w s.t w^N = 1 mod M and there are not other numbers k (k < N) 
    s.t w^k = 1 mod M 

    """
    roots = nthroot_mod(1, N, M, True)[1:]  # find Nth root of unity
    for root in roots:  # find primitive Nth root of unity
        is_primitive = True
        for k in range(1, N):
            if pow(root, k, M) == 1:
                is_primitive = False
        if is_primitive:
            return root
    return None


def isPrimitiveNthRoot(M, N, beta):
    """
    verify B^N = 1 (mod M)
    """
    return pow(beta, N, M) == 1     # modular(M).modExponent(beta, N) == 1


def uniform_sample(upper, num):
    """
    Sample num values uniformly between [0,upper).
    """
    sample = []
    for i in range(num):
        value = random.randint(0, upper - 1)
        sample.append(value)
    return sample


def gauss_sample(num, stdev):
    """
    Sample num values from gaussian distribution
    mean = 0, standard deviation = stdev 
    """
    sample = np.random.normal(0, stdev, num)
    sample = sample.round().astype(int)
    return sample


def hamming_sample(num, hwt):
    """
    Sample a vector uniformly at random from -1, 0, +1,
    subject to the condition that it has exactly hwt nonzero entries. 
    """
    i = 0
    sample = [0] * num
    while i < hwt:
        degree = random.randint(0, num - 1)
        if sample[degree] == 0:
            coeff = random.randint(0, 1)
            if coeff == 0:
                coeff = -1
            sample[degree] = coeff
            i += 1
    return sample


def small_sample(num):
    """
    Sample vectors with entires -1, 0, +1.
    Each element is 0 with probabilty 0.5 and +-1 with probabilty 0.25. 
    """
    sample = [0] * num
    for i in range(num):
        u = random.randint(0, 3)
        if u == 3:
            sample[i] = -1
        if u == 2:
            sample[i] = 1
    return sample

# Not used


class modular:
    def __init__(self, M):
        self.mod = M
        self.M_bit = M.bit_length()
        self.u = (1 << (2 * self.M_bit)) // M

    def modReduce(self, x):
        """
        Barrett modular reduction algorithm.
        Compute x mod M. M is initialized by modular class.

        Examples
        ========

        >>> from numTh import numTh
        >>> modular(11).modReduce(12)
        1
        """

        assert 0 <= x < pow(self.mod, 2), 'out of range.'
        q = (x * self.u) >> (2 * self.M_bit)
        r = x - q * self.mod
        while r >= self.mod:
            r -= self.mod
        return r

    def modReducem(self, x, M):
        """
        Barrett modular reduction algorithm.
        Compute x mod M. M can be redefined.

        Examples
        ========

        >>> from numTh import numTh
        >>> modular(5).modReducem(12, 11)
        1
        """
        tmp_mod, tmp_M_bit, tmp_u = self.mod, self.M_bit, self.u
        self.mod = M
        self.M_bit = M.bit_length()
        self.u = (1 << (2 * self.M_bit)) // M
        r = self.modReduce(x)
        # return initial modular, bit size of modular and precompute u
        self.mod, self.M_bit, self.u = tmp_mod, tmp_M_bit, tmp_u
        return r

    def modInv(self, x):
        """
        Calculate modular inverse.

        Examples
        ========

        >>> from numTh import modular  
        >>> modular(5).modInv(3)
        2
        """
        t, new_t, r, new_r = 0, 1, self.mod, x

        while new_r != 0:
            q = r // new_r
            r, new_r = new_r, (r % new_r)
            t, new_t = new_t, (t - q * new_t)
        assert r <= 1, 'x is not invertible'
        return t if t > 0 else t + self.mod

    # Slower than pow(base, power, modulus)
    def modExponent(self, base, power):
        """
        Modular exponentiation algorithm.
        It's a fast method to compute a^b mod p

        Examples
        ========

        >>> from numTh import modular 
        >>> modular(2013265921).modExponent(1003203377, 2048)
        1
        """
        result = 1
        power = int(power)
        base = base % self.mod
        while power > 0:
            if power & 1:
                # self.modReduce(result * base)
                result = result * base % self.mod
            base = base * base % self.mod  # self.modReduce(base * base)
            power = power >> 1
        return result
