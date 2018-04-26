## PyFHE
PyFHE is an implementation of the [Brakerski-Gentry-Vaikuntanathan][1] (BGV) scheme, along with some optimizations described in the [Gentry-Halevi-Smart][2] optimizations. This project takes [HElib][3] as reference but is simpler than HElib. It only implements some basic homomorphic operations such as multiplication and addition and does not use ciphertext packing techniques.

This project is written in Python 3.6.5 and uses Python packages [sympy][4] and [numpy][5].

## Demo
Run `Test_Ctxt.py`, which executes a 4 levels 4096 degree polynomial BGV encryption.
```
  > python Test_Ctxt.py
```
## TODO
Use multiprocessing.

  [1]: http://eprint.iacr.org/2011/277       "BGV12"
  [2]: http://eprint.iacr.org/2012/099       "GHS12"
  [3]: https://github.com/shaih/HElib        "HELIB" 
  [4]: http://www.sympy.org/en/index.html    "SYMPY"
  [5]: http://www.numpy.org                  "NUMPY"
