import numpy as np
from CRTPoly import CRTPoly
from numTh import uniform_sample, gauss_sample, hamming_sample, small_sample


class FHE:
    """
    Implementation of BGV-FHE scheme. 
    """

    def __init__(self, d, stdev, primes, P, L, cur_level=0):
        """
        Initialize parameters.
        L : limitation of homomorphic multiplications
        cur_level : homomorphic multiplication times
        d : polynomial degree
        stdev : standard deviation of gaussian distribution
        prime_set : total primes
        modulus : product of primes 
        """
        self.L = L          # levels
        self.cur_level = cur_level  # current level
        self.d = d          # polynomial degree
        self.stdev = stdev  # standard deviation
        self.prime_set = list(primes)   # primes
        self.prime_set.sort(reverse=True)
        self.special_prime = P
        self.modulus = 1    # modulus
        for i in range(cur_level, L):
            self.modulus *= primes[i]

    def setCoeffs(self, poly, q=None):
        """
        Let each coefficient in the polynomial in the range of [-q/2,q/2].
        """
        if q is None:
            q = self.modulus
        for i, coeff in enumerate(poly):
            if coeff > q // 2:
                poly[i] -= q

    def secretKeyGen(self, h):
        """
        Generate secret key. 
        sk = (1, s')
        """
        secret_key = []
        sk0 = [0] * self.d
        sk0[0] = 1
        sk1 = hamming_sample(self.d, h)
        secret_key.append(sk0)
        secret_key.append(sk1)
        return secret_key

    def publicKeyGen(self, sk, modulus=None):
        """
        Generate public key. 
        pk = (b, -A'), b = A's'+2e. 
        """
        prime_set = list(self.prime_set)
        if modulus is None:
            modulus = self.modulus
        else:
            prime_set.append(self.special_prime)

        public_key = []
        e = gauss_sample(self.d, self.stdev)
        A = uniform_sample(modulus, self.d)
        # set coefficients range [-q/2,q/2]
        self.setCoeffs(A, modulus)
        # CRT-FFT representation
        fft_sk1 = CRTPoly(sk[1], prime_set)
        fft_A = CRTPoly(A, prime_set)
        fft_2e = CRTPoly((2 * np.asarray(e)).tolist(), prime_set)
        # b = A's'+2e
        fft_b = fft_A * fft_sk1 + fft_2e
        # polynomial representation
        b = fft_b.toPoly()
        # set coefficients range [-q/2,q/2]
        self.setCoeffs(b, modulus)
        neg_A = (-(np.asarray(A))).tolist()
        public_key.append(b)
        public_key.append(neg_A)
        return public_key

    def switchKeyGen(self, sk):
        """
        Generate L-1 switch keys. 
        Each switch key is in R_Qi, where Qi = P * modulus_i and i is level.
        And each switch key is (b + P * s^2, -a), 
        where b = a * s + 2e, and a is sampled uniformly in [-Qi/2,Qi/2].
        """
        modulus = self.modulus * self.special_prime
        prime_set = list(self.prime_set)
        prime_set.append(self.special_prime)
        switch_keys = []
        switch_key = []
        for i in range(0, self.L - 1):
            switch_key = []
            if i != 0:
                modulus //= self.prime_set[i - 1]
            pk = self.publicKeyGen(sk, modulus)  # pk = (a * s + 2e, -a)
            # CRT-FFT representation
            crt_b = CRTPoly(pk[0], prime_set[i:])
            crt_sk1 = CRTPoly(sk[1], prime_set[i:])
            # key0 = a * s + 2e + P * s^2
            crt_switch_key0 = crt_b + crt_sk1 * crt_sk1 * self.special_prime
            # polynomial representation
            key0 = crt_switch_key0.toPoly()
            # set coefficients in [-Q/2, Q/2]
            self.setCoeffs(key0, modulus)
            # switch key = (b + P * s^2, -a)
            switch_key.append(key0)
            switch_key.append(pk[1])
            switch_keys.append(switch_key)
        return switch_keys

    def homoEnc(self, m, pk):
        """
        FHE encryption:
        c = (c0, c1)
        c0 = pk0 * r + 2e0 + m 
        c1 = pk1 * r + 2e1
        """
        r = small_sample(self.d)
        e0 = gauss_sample(self.d, self.stdev)
        e1 = gauss_sample(self.d, self.stdev)
        if len(m) < self.d:
            m += [0] * (self.d - len(m))
        # CRT-FFT representation
        crt_m = CRTPoly(m, self.prime_set)
        crt_pk0 = CRTPoly(pk[0], self.prime_set)
        crt_pk1 = CRTPoly(pk[1], self.prime_set)
        crt_r = CRTPoly(r, self.prime_set)
        crt_2e0 = CRTPoly((2 * np.asarray(e0)).tolist(), self.prime_set)
        crt_2e1 = CRTPoly((2 * np.asarray(e1)).tolist(), self.prime_set)
        # c0 = pk0 * r + 2e0 + m, c1 = pk1 * r + 2e1
        crt_c0 = crt_m + crt_2e0
        crt_c1 = crt_2e1
        crt_c0 += crt_pk0 * crt_r
        crt_c1 += crt_pk1 * crt_r
        # polynomial representation
        c0 = crt_c0.toPoly()
        c1 = crt_c1.toPoly()
        # set coefficients in [-q/2,q/2]
        self.setCoeffs(c0)
        self.setCoeffs(c1)
        c = []
        c.append(c0)
        c.append(c1)
        return c

    def homoDec(self, c, sk):
        """
        FHE decryption:
        m = (c0 + c1 * s') mod 2
        """
        # CRT-FFT representation
        crt_c0 = CRTPoly(c[0], self.prime_set[self.cur_level:])
        crt_c1 = CRTPoly(c[1], self.prime_set[self.cur_level:])
        crt_sk1 = CRTPoly(sk[1], self.prime_set[self.cur_level:])
        # m = [[<c,s>] mod p] mod 2
        crt_m = crt_c0 + crt_c1 * crt_sk1
        # polynomial represenatation
        m = crt_m.toPoly()
        # set coefficients range [-q/2,q/2]
        self.setCoeffs(m)
        return np.remainder(m, 2).tolist()

    def scale(self, c, from_q, to_q):
        """
        Change the modulus. 
        c = p_t * qoutient + rem
        odd number coefficients in rem +- p_t (i.e +- 1 becomes even)
        _c = p_t * (qoutient +- 1) + rem
        the coefficient is even or odd is effected by rem
        result = (c - _c) / p_t = quotient +- 1
        """
        p_t = from_q // to_q
        _c = np.asarray(c) % p_t
        for i, _c_i in enumerate(_c):
            self.setCoeffs(_c_i, p_t)
            for j, coeff in enumerate(_c_i):
                if coeff % 2 == 1:
                    if coeff > 0:
                        _c[i][j] -= p_t
                    else:
                        _c[i][j] += p_t
        c_dagger = np.asarray(c) - _c
        result = c_dagger // p_t
        return np.remainder(result, to_q,).tolist()

    def modSwitch(self, c, level):
        """
        Scale modulus down. c' is closest to c/p. 
        And new c' must satisfy c' = c mod 2. 
        """
        assert level < self.L - 1, "cannot reduce noise"
        # new modulus
        to_modulus = self.modulus // self.prime_set[level]
        # scale(c, q, q')
        result = self.scale(c, self.modulus, to_modulus)
        # down to new modulus
        self.modulus = to_modulus
        # increase level
        self.cur_level += 1
        return np.remainder(result, self.modulus).tolist()

    def keySwitch(self, c, switch_key):
        """
        Key switching. 
        """
        modulus = self.modulus * self.special_prime
        prime_set = list(self.prime_set[self.cur_level:])
        prime_set.append(self.special_prime)
        # CRT-FFT representation
        crt_c0 = CRTPoly(c[0], prime_set)
        crt_c1 = CRTPoly(c[1], prime_set)
        crt_c2 = CRTPoly(c[2], prime_set)
        crt_b = CRTPoly(switch_key[0], prime_set)
        crt_a = CRTPoly(switch_key[1], prime_set)
        # c'0 = P * c0 + b * c2, c'1 = P * c1 + a * c2
        crt_result0 = crt_c0 * self.special_prime + crt_b * crt_c2
        crt_result1 = crt_c1 * self.special_prime + crt_a * crt_c2
        # polynomial representation
        result0 = crt_result0.toPoly()
        result1 = crt_result1.toPoly()
        # set coefficients in [-Q/2, Q/2]
        self.setCoeffs(result0, modulus)
        self.setCoeffs(result1, modulus)

        result = []
        result.append(result0)
        result.append(result1)
        # scale(c', Q, q)
        result = self.scale(result, modulus, self.modulus)
        return result
