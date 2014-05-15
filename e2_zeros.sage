# Zeros of E2
# -----------
#
# Usage example:
#
# getE2Zero(7/100)
#
# Returns the zero of E2 around the cusp 7/100, note that the zero is
# very close to 7/100+pi/6/100^2*i.


import sage.libs.mpmath.all as mpmath
from mpmath import mp,mpmathify

# Returns the matrix that maps the usual fundamental domain to the fundamental
# domain that contains a zero of E2 around the cusp q=a/c, a,c coprime, c>0.
def cuspMatrix(q):
  qrat=QQ(q);
  a=q.numer();
  c=q.denom();
  (g,d,b)=a._xgcd(-c,minimal=True);
  if ((abs(c)-d)<d):
    d-=sign(c)*c;   
    b-=sign(c)*a;   

  A=Gamma0(1)([a,b,c,d]);
  return A;

# Get the zero in the fundamental domain with x value close to -d/c
def getWcdZero(c,d):
  # Since we evaluate in the fundamental domain the usual definition of E2 is ok
  E2=lambda z: 1-24*mpmath.nsum(lambda n: n/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');
  F=lambda w: E2(w)-i*6/pi/(w+d/c);
  west=-d/c+6/pi*i;

  return mpmath.findroot(F,west,solver='muller');

# Get the zero in the fundamental domain corresponding to the cusp q
def getWZero(q):
  A=cuspMatrix(Rational(q))
  return getWcdZero(A.c(),A.d());

# Returns the unique zero of E2 around the cusp q...
def getE2Zero(q):
  A=cuspMatrix(q)
  c=A.c();
  d=A.d();

  # Since we evaluate in the fundamental domain the usual definition of E2 is ok
  E2=lambda z: 1-24*mpmath.nsum(lambda n: n/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');
  F=lambda w: E2(w)-i*6/pi/(w+d/c);
  west=-d/c+6/pi*i;
  with mp.extradps(floor(mp.log((2*(c+1))^2)/mp.log(10))+1):
    wzero=mpmath.findroot(F,west,solver='muller');
    zest=q+pi/6/c/c*i;
    zero=A.acton(wzero);

  return zero;
