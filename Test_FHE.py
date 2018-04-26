from FHE import FHE
from CRTPoly import CRTPoly
from numTh import findPrimes
import numpy as np


def multiply(c1, c2, primes):
    result = []
    fft_c10 = CRTPoly(c1[0], primes)
    fft_c11 = CRTPoly(c1[1], primes)
    fft_c20 = CRTPoly(c2[0], primes)
    fft_c21 = CRTPoly(c2[1], primes)
    fft_result0 = fft_c10 * fft_c20
    fft_result1 = fft_c10 * fft_c21 + fft_c11 * fft_c20
    fft_result2 = fft_c11 * fft_c21
    result.append(fft_result0.toPoly())
    result.append(fft_result1.toPoly())
    result.append(fft_result2.toPoly())
    return result


def polyMul(p1, p2, primes):
    fft_p1 = CRTPoly(p1, primes)
    fft_p2 = CRTPoly(p2, primes)
    modulus = 1
    for prime in primes:
        modulus *= prime
    fft_result = fft_p1 * fft_p2
    result = fft_result.toPoly()
    for i, coeff in enumerate(result):
        if coeff > modulus // 2:
            result[i] -= modulus
    return np.remainder(result, 2).tolist()


poly_degree = 4096
stdev = 3.2
L = 4
#primes = [549755860993, 549755873281, 549755876353]
primes, bits = findPrimes(22, 4096, 4)
a, bits = findPrimes(10, 4096, 1)
P = a[0]
# primes = [521, 569, 577]
modulus = 1
for prime in primes:
    modulus *= prime
f = FHE(poly_degree, stdev, primes, P, L)
sk = f.secretKeyGen(64)
# sk = [[1, 0, 0, 0], [0, 1, -1, 0]]
pk = f.publicKeyGen(sk)
# pk = [[-24187115, -62847359, 2213875, 53855074], [-13973837, -16187706, -70042772, 76821192]]
switch_keys = f.switchKeyGen(sk)
m = np.random.randint(0, 2, 4096).tolist()
m1 = np.random.randint(0, 2, 4096).tolist()
print('plaintext')
# print(m)
# print(m1)
print('Encryption')
c = f.homoEnc(m, pk)
print('homo Multiply')
mul_result = multiply(c, c, primes)
mul_result = f.keySwitch(mul_result, switch_keys[0])
mul_result = f.modSwitch(mul_result, 0)
dec_mul_result = f.homoDec(mul_result, sk)
print('Decrypt mul result')
#print(polyMul(m, m, primes))
# print(dec_mul_result)
print(dec_mul_result == polyMul(m, m, primes))
dec_mm = f.homoDec(c, sk)
print('Decrypt m')
# print(dec_mm)
c1 = f.homoEnc(m1, pk)
"""
print('Modulus Switching')
c = f.modSwitch(c, 0)
c = f.modSwitch(c, 1)
"""
dec_mm = f.homoDec(c, sk)
print('Decrypt m')
print(dec_mm == m)
new_c = ((np.asarray(c1) + np.asarray(c)) % f.modulus).tolist()
print('Decryption')
#dec_m = f.homoDec(new_c,sk)


# if m == dec_m:
#    print('success')
# else:
#    print('fail')
#m = []
# for bit in dec_m:
#    m.append(int(bit))
# print(m)
