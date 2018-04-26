############################################################
#                                                          #
# File Name: NTT.py                                        #
#                                                          #
# Author: Jyun-Neng Ji (jyunnengji@gmail.com)              #
#                                                          #
# Creation Date: 2018/02/15                                #
#                                                          #
# Last Modified: 2018/03/25                                #
#                                                          #
# Description: Implement Cooley-Tukey FFT in finite field. #
#                                                          #
############################################################
from numTh import *
from sympy.core.numbers import mod_inverse
from sympy.ntheory import sqrt_mod


class NTT:
    def __init__(self, poly, M, N, ideal=True, ntt=False, w=None, phi=None):
        """
        Initialize the parameters for NTT 
        and transfer poly into frequency domain if ntt is True.
        Parameters: mod(modulus), N(NTT points), w(root of unity), phi(square root of unity),
                    ideal(ideal ring), fft_poly(poly in freq domain). 
        """
        if ntt:  # poly is in frequency domain
            self.initial_w_ntt(poly, M, N, ideal, w, phi)
        else:   # poly is not in frequency domain
            self.initial_wo_ntt(poly, M, N, ideal, w, phi)

    def initial_wo_ntt(self, poly, M, N, ideal, w, phi):
        self.mod = M    # modulus
        self.N = N      # N points NTT
        if w is None and phi is None:   # generate root of unity and square root of unity
            self.w = findPrimitiveNthRoot(M, N)  # primitive Nth root of unity
            self.phi = sqrt_mod(self.w, M)  # phi: w^2 = 1 (mod M)
        else:
            self.w = w
            self.phi = phi
        self.ideal = ideal
        # If it is computed over Z[x]/<x^n+1>, it should multiply phi
        if ideal:
            poly_bar = self.mulPhi(poly)
        else:
            poly_bar = poly
        self.fft_poly = self.ntt(poly_bar)

    def initial_w_ntt(self, poly, M, N, ideal, w, phi):
        self.mod = M    # modulus
        self.N = N      # N points NTT
        self.w = w      # primitive Nth root of unity
        self.phi = phi  # phi: w^2 = 1 (mod M)
        self.ideal = ideal
        self.fft_poly = poly

    def __name__(self):
        return "NTT"

    def __str__(self):
        poly = " ".join(str(coeff) for coeff in self.fft_poly)
        return "NTT points [" + poly + "] modulus " + str(self.mod)

    def __mul__(self, other):
        """
        Multiply in frequency domain. 
        """
        assert type(other).__name__ == 'NTT' or type(other).__name__ == 'int', 'type error'
            

        if type(other).__name__ == 'int':
            mul_result = self.mulConstant(other)

        else:
            assert self.N == other.N, "points different"
            assert self.mod == other.mod, "modulus different"
            assert self.ideal == other.ideal

            mul_result = []
            for i, point in enumerate(self.fft_poly):
                mul_result.append((point * other.fft_poly[i]) % self.mod)

        return NTT(mul_result, self.mod, self.N, self.ideal, True, self.w, self.phi)

    def __add__(self, other):
        """
        Addition in  frequency domain
        """
        assert self.N == other.N, 'points different'
        assert self.mod == other.mod, 'modulus different'
        assert self.ideal == other.ideal

        add_result = []
        for i, point in enumerate(self.fft_poly):
            add_result.append((point + other.fft_poly[i]) % self.mod)

        return NTT(add_result, self.mod, self.N, self.ideal, True, self.w, self.phi)

    def __sub__(self, other):
        """
        Substraction in frequency domain
        """
        assert self.N == other.N, 'points different'
        assert self.mod == other.mod, 'modulus different'
        assert self.ideal == other.ideal

        sub_result = []
        for i, point in enumerate(self.fft_poly):
            sub_result.append((point - other.fft_poly[i]) % self.mod)

        return NTT(sub_result, self.mod, self.N, self.ideal, True, self.w, self.phi)

    def bitReverse(self, num, len):
        """
        Reverse bits of a number.

        Examples
        ========

        >>> NTT.bitReverse(3, 3) # 3(011)
        6                        # 6(110)   
        """
        rev_num = 0

        for i in range(0, len):

            if (num >> i) & 1:
                rev_num |= 1 << (len - 1 - i)

        return rev_num

    def orderReverse(self, poly, N_bit):
        """
        Change the order of coefficients of polynomial to fit DIT FFT input.

        Examples
        ========

        >>> NTT.orderReverse([1, 2, 3, 4], 2)
        [1, 3, 2, 4]
        """
        _poly = list(poly)
        for i, coeff in enumerate(_poly):

            rev_i = self.bitReverse(i, N_bit)

            if rev_i > i:
                coeff ^= _poly[rev_i]
                _poly[rev_i] ^= coeff
                coeff ^= _poly[rev_i]
                _poly[i] = coeff

        return _poly

    def ntt(self, poly, w=None):
        """
        Compute FFT in finite feild.
        Use Cooley-Tukey DIT FFT algorithm.
        input: polynomial, primitive nth root of unity
        output FFT(poly)
        complexity: O(N/2 log N)
        """

        if w is None:
            w = self.w

        N_bit = self.N.bit_length() - 1
        rev_poly = self.orderReverse(poly, N_bit)

        for i in range(0, N_bit):

            points1, points2 = [], []

            for j in range(0, int(self.N / 2)):
                shift_bits = N_bit - 1 - i
                P = (j >> shift_bits) << shift_bits
                w_P = pow(w, P, self.mod)
                odd = rev_poly[2 * j + 1] * w_P
                even = rev_poly[2 * j]
                points1.append((even + odd) % self.mod)
                points2.append((even - odd) % self.mod)
                points = points1 + points2

            if i != N_bit:
                rev_poly = points

        return points

    def intt(self):
        """
        Compute IFFT in finite feild.
        The algorithm is the same as NTT, but it need to change w to w^(-1) 
        and multiply N^(-1) for each coefficient of polynomial.
        input: FFT(poly), primitive nth root of unity
        output: polynomial
        complexity: (N/2 log N)
        """

        inv_w = mod_inverse(self.w, self.mod)
        inv_N = mod_inverse(self.N, self.mod)

        poly = self.ntt(self.fft_poly, inv_w)

        for i in range(0, self.N):
            poly[i] = poly[i] * inv_N % self.mod

        # if it is computed over Z[x]/<x^n+1>, it should multiply phi^(-1)
        if self.ideal:
            inv_phi = mod_inverse(self.phi, self.mod)
            poly = self.mulPhi(poly, inv_phi)

        return poly

    def mulPhi(self, poly, phi=None):
        """
        The coefficients in polynomial multiply phi^i mod modulus
        a_i_bar = a_i * phi^i (mod m)

        It's used before NTT multiplication over Zp/<x^n+1>
        """

        if phi is None:
            phi = self.phi

        poly_bar = list(poly)

        for i, coeff in enumerate(poly):
            poly_bar[i] = (poly[i] * pow(phi, i, self.mod)) % self.mod

        return poly_bar

    def mulConstant(self, constant):
        mul_result = []
        for coeff in self.fft_poly:
            result = coeff * constant % self.mod
            mul_result.append(result)
        return mul_result
