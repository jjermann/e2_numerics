# Fast evaluation of Eisenstein series E2, E4, E6
# -----------------------------------------------
#
# Usage example:
#
# FE2(5.12+0.00013*I)
#
# Returns the numerical value of E2(5.12+0.00013*I) using smart modular evaluation


import sage.libs.mpmath.all as mpmath
from mpmath import mp,mpmathify

# Returns [A,w] such that A.acton(w)=z and such that w lies in the usual
# strict fundamental domain of SL2(Z)
def getFD(z):
  A=Gamma0(1)([1,0,0,1]);
  T=Gamma0(1)([1,1,0,1]);
  S=Gamma0(1)([0,-1,1,0]);
  w=z;
  while (abs(w)<1 or abs(w.real())>1/2):
    if (abs(w)<1):
      w=S.acton(w);
      A=S*A;
    if (w.real()>=1/2):
      w=(T^(-1)).acton(w);
      A=(T^(-1))*A;
    elif (w.real()<1/2):
      w=T.acton(w);
      A=T*A;
  if (w.real()==1/2):
    w=(T^(-1)).acton(w);
    A=(T^(-1))*A;
  if (abs(w)==1 and w.real()>0):
    w=S.acton(w);
    A=S*A;
  return [A^(-1),A.acton(z)];

# Eisenstein series of weight 2 (efficient)
def FE2(z):
  # The usual definition (slow for Im(z) small).
  E2=lambda z: 1-24*mpmath.nsum(lambda n: n/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');

  # A prelimenary step to get c which we need to estimate the needed precision.
  with mp.workdps(100):
    c=abs((getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec))[0])[1][0]);

  # Regarding the precision the result is very roughly of the form 4c^2*E2(w),
  # so we want to do calculations with the appropriate additional precision.
  with mp.extradps(floor(mp.log((2*(c+1))^2)/mp.log(10))+1):
    [A,w]=getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec));
    c=A[1][0];
    d=A[1][1];
    w=mpmathify(w);
    Ew=E2(w);

    # We use the transformation formula for E2 to do calculations in the usual
    # fundamental domain. This is fast since Im(w) > sqrt(3)/2 is "big".
    return Ew*(c*w+d)^2+6*c/pi/i*(w*c+d);

# Eisenstein series of weight 4 (efficient)
def FE4(z):
  # The usual definition (slow for Im(z) small).
  E4=lambda z: 1+240*mpmath.nsum(lambda n: n^3/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');

  # A prelimenary step to get c which we need to estimate the needed precision.
  with mp.workdps(100):
    c=abs((getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec))[0])[1][0]);

  # Regarding the precision the result is very roughly of the form const*c^4*E4(w),
  # so we want to do calculations with the appropriate additional precision.
  with mp.extradps(floor(mp.log((2*(c+1))^4)/mp.log(10))+1):
    [A,w]=getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec));
    c=A[1][0];
    d=A[1][1];
    w=mpmathify(w);
    Ew=E4(w);

    # We use the transformation formula for E4 to do calculations in the usual
    # fundamental domain. This is fast since Im(w) > sqrt(3)/2 is "big".
    return Ew*(c*w+d)^4;

# Eisenstein series of weight 6 (efficient)
def FE6(z):
  # The usual definition (slow for Im(z) small).
  E6=lambda z: 1-504*mpmath.nsum(lambda n: n^5/(mpmath.exp(-2*pi*i*n*z)-1),[1,mpmath.inf],method='d');

  # A prelimenary step to get c which we need to estimate the needed precision.
  with mp.workdps(100):
    c=abs((getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec))[0])[1][0]);

  # Regarding the precision the result is very roughly of the form const*c^6*E6(w),
  # so we want to do calculations with the appropriate additional precision.
  with mp.extradps(floor(mp.log((2*(c+1))^6)/mp.log(10))+1):
    [A,w]=getFD(mpmath.mpmath_to_sage(mpmathify(z),mp.prec));
    c=A[1][0];
    d=A[1][1];
    w=mpmathify(w);
    Ew=E6(w);

    # We use the transformation formula for E6 to do calculations in the usual
    # fundamental domain. This is fast since Im(w) > sqrt(3)/2 is "big".
    return Ew*(c*w+d)^6;
