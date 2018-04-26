from Ctxt import Ctxt
from FHE import FHE
from numTh import *
import numpy as np
import cProfile as cp
import time

def setup(d, stdev, bits_per_level, P_bits, L):
    primes, total_bits = findPrimes(bits_per_level, d, L)
    special_primes, total_bits = findPrimes(P_bits, d, 1)
    P = special_primes[0]
    f = FHE(d, stdev, primes, P, L)
    return f, primes, P

def keyGen(f, h):
    print('-----Secret Key Generation-----')
    sk = f.secretKeyGen(h)
    print('-----Public Key Generation-----')
    pk = f.publicKeyGen(sk)
    print('-----Switch Key Generation-----')
    swk = f.switchKeyGen(sk)
    return sk, pk, swk 

def set_timer_profile():
    pr = cp.Profile()
    pr.enable()
    start = time.time()
    return pr, start  

def end_timer_profile(pr, start, filename):
    cost = time.time() - start 
    pr.disable()
    pr.dump_stats(filename)
    return cost 

def hex_rep(m):
    binary_str = ''.join(str(bit) for bit in m)
    return hex(int(binary_str, 2))

key_profile = './profile/keygeneration.prof'
enc_profile = './profile/encrytpion.prof'
dec_profile = './profile/decryption.prof'
l0_mul_profile = './profile/l0_mul.prof'
l1_mul_profile = './profile/l1_mul.prof'
add_profile = './profile/add.prof'
f, primes, P = setup(4096, 3.2, 22, 49, 4)

# print('-----plaintext-----')
msg1 = np.random.randint(0, 2, 4096).tolist()
print(hex_rep(msg1))
pr, start = set_timer_profile()

sk, pk, swk = keyGen(f, 64)
cost = end_timer_profile(pr, start, key_profile)
print('Done {0: .3f}s'.format(cost))
print('-----Encryption-----')
c1 = Ctxt(4096, 3.2, primes, P, 4)

pr, start = set_timer_profile()
# message encryption
c1.enc(msg1, pk)

cost = end_timer_profile(pr, start, enc_profile)

print('Done {0: .3f}s'.format(cost))
print('-----Homomorphic Operation-----')
pr, start = set_timer_profile()

c_mul = c1 * c1
c_mul.relinearize(swk[0])

cost = end_timer_profile(pr, start, l0_mul_profile)
print('First multiplication done {0: .3f}s'.format(cost))
pr, start = set_timer_profile()

c_mul *= c_mul
c_mul.relinearize(swk[1])

cost = end_timer_profile(pr, start, l1_mul_profile)
print('Second multiplication done {0: .3f}s'.format(cost))

pr, start = set_timer_profile()

c_add = c1 + c1

cost = end_timer_profile(pr, start, add_profile)
print('Addition done {0: .3f}s'.format(cost))

print('-----Decryption-----')
pr, start = set_timer_profile()

m1 = c1.dec(sk)

cost = end_timer_profile(pr, start, dec_profile)
print('Done {0: .3f}s'.format(cost))
mul_result = c_mul.dec(sk)
add_result = c_add.dec(sk)
print(msg1 == m1)
#print('msg1 ', m1)
#print('msg1^4 = ', mul_result)
#print('msg1 + msg1 = ', add_result)
