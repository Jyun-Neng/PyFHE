from CRTPoly import CRTPoly
from numTh import *
import numpy as np

poly = [12, 33, 14, 53]
poly2 = [8, 7, 6, 7]
primes = [7, 11, 13]
crt_poly = CRTPoly(poly, primes, fft=False)
crt_poly2 = CRTPoly(poly2, primes, fft=False)
print(crt_poly)
mul_const_result = crt_poly2 * 300
add_result = crt_poly2 + crt_poly
sub_result = crt_poly - crt_poly2
print(add_result.toPoly())
print(sub_result.toPoly())
print(mul_const_result.toPoly())

