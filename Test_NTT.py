from NTT import NTT
from numTh import *
import numpy as np 

poly1 = [1,0,1,1]
poly2 = [1,0,1,1]
modulus = 17
N = 4
fft_poly1 = NTT(poly1, modulus, N)
fft_poly2 = NTT(poly2, modulus, N)
print(type(fft_poly1).__name__)
mult_result = fft_poly1 * fft_poly2
add_result = fft_poly1 + fft_poly2
mult_const_result = fft_poly1 * 2
sub_result = fft_poly1 - fft_poly2
print('poly1=', poly1)
print('poly2=', poly2)
print('multiply:', mult_result.intt())
print('multiply 2:', mult_const_result.intt())
print('addition:', add_result.intt())
print('substraction:', sub_result.intt())

