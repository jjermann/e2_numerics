# Numerical Fourier coefficients of F(z)+z
# ----------------------------------------
#
# Usage example:
#
# num_fourier_coeff(1)
#
# Returns the first non-constant Fourier coefficient of f^(-1)(z)-z


import sage.libs.mpmath.all as mpmath
from mpmath import mp,mpmathify

# Fourier expansion of f^(-1)(z)-z...
def num_fourier_coeff(n):
  E2=lambda z: 1-24*mpmath.nsum(lambda n: n/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');
  f=lambda z: lambda w: w-6*i/pi/E2(w)-z;
  F=lambda z: mpmath.findroot(f(z),z+6*i/pi,solver='muller');
  H=lambda z: F(z)-z;
  fourier_fun=lambda n: lambda x: H(x)*mpmath.exp(-2*pi*i*n*x);
  coeff=mpmath.quad(fourier_fun(n), [-1/2,1/2], method='gauss-legendre');
  return coeff;
