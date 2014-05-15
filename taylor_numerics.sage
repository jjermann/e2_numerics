# Taylor expansion of f^(-1)(z)
# -----------------------------
#
# Usage example 1:
#
# taylor_coeff(1)
# 
# Returns the first Taylor coefficient of f^(-1) around 0
# as a rational function over Q in x=E4(w), y=E6(w), z=E2(w), c=2*pi*i
# where w is the zero of E2 on the imaginary axis.
#
#
# Usage example 2:
#
# num_taylor_pol(3)
#
# Returns the Taylor polynomial of order 3 of f^(-1) around 0
# (as a numerical function).


import sage.libs.mpmath.all as mpmath
from mpmath import mp,mpmathify

load derivative.sage
load fast_eisenstein.sage
load e2_zeros.sage

#we need one more variable "w"
R.<x,y,z,t,c,w> = PolynomialRing(QQ,"x,y,z,t,c,w");

# Taylor coefficient of F, where F(z) is the inverse of f(w)=w-6*i/(pi*E2(w)).
# It is given as a rational function over Q in x=E4(w), y=E6(w), z=E2(w)=-12/c/w, c=2*pi*i,
# where w is the zero of E2 on the imaginary axis.
def taylor_coeff(n):
  if (n==0):  
    pfinal=-12/c/z;
  else:
    pderiv=dmod(D^(n-1),1/factorial(n)*((t-w)/(t+12/c/z))^n);
    pdenom=pderiv.denominator();
    step=0;
    while (pdenom.subs(t=-12/c/z,w=-12/c/z)==0):
      step=step+1;
      pnum=dmod(D,pderiv.numerator());
      pdenom=dmod(D,pderiv.denominator());
      pderiv=(pnum/pdenom);
      pderiv.reduce();
    pfinal=pderiv.subs(t=-12/c/z,w=-12/c/z);  
    pfinal.reduce();
  return pfinal


# Substitute numerical values for the exact formula of the (given) Taylor coefficient.
def num_taylor_eval(pfinal):
  ts=lambda z: mpmath.mpmath_to_sage(mpmathify(z),mp.prec)
  w0=getWZero(0);
  z4=FE4(w0);
  z6=FE6(w0);
  return pfinal.subs(x=ts(z4),y=ts(z6),z=ts(6*i/pi/w0),c=2*pi*i).n(prec=mp.prec)

# Returns the Taylor series as a numerical function for the given coefficients.
def fseries(V,z0=0):
  return lambda z: mpmath.nsum(lambda n: V[Integer(mpmath.mpmath_to_sage(n,mp.prec))]*(z-z0)^n,[0,len(V)-1]);

# Taylor polynomial of F(z) around 0 where F(z) is the inverse of f(w)=w-6i/(pi*E2(w))
def num_taylor_pol(m):
  V=[num_taylor_eval(taylor_coeff(n)) for n in (0..m)];  
  #V[1]-=1; 
  return fseries(V); 

