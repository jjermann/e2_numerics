# Derivatives of modular functions/forms
# --------------------------------------
#
# Usage example:
#
# dmod(D^4,x^3/(x^3-y^2+t))
#
# This gives the 4'th derivative of E4^3(t)/(E4^3(t)-E6^2(t)-t)
# as a function in x=E4(t), y=E6(t), z=E2(t), t, c=2*pi*i.
#
# If "x,y,z,t,c=var('x,y,z,t,c')" instead of using a polynomial ring,
# then arbitrary functions can be given as an argument...


# Classical polynomial algebra
R.<x,y,z,t,c> = PolynomialRing(QQ,"x,y,z,t,c");

# Algebra of differential operators
# X,Y,Z,T correspond to E4(t),E6(t),E2(t),t
# C=2*pi*i is a constant
A.<X,Y,Z,T,dX,dY,dZ,dT,C>=FreeAlgebra(QQ,9);
G=A.g_algebra({dX*X:1+X*dX,dY*Y:1+Y*dY,dZ*Z:1+Z*dZ,dT*T:1+T*dT});
(X,Y,Z,T,dX,dY,dZ,dT,C) = G.gens();

# Specific differential operator corresponding to the derivative
D=C*(1/3*(X*Z-Y)*dX+1/2*(Y*Z-X^2)*dY+1/12*(Z^2-X)*dZ)+dT;

# Scalar multiplication on the D-module
#
# Note that this works with arbitrary differentiable functions.
#
def dmod(A,b):
  L=A.monomials();
  RES=0;
  for s in L:
    TEMPRES=b;
    degs=s.degrees();
    for i in range(degs[4]):
      TEMPRES=TEMPRES.derivative(x);
    for i in range(degs[5]):
      TEMPRES=TEMPRES.derivative(y);
    for i in range(degs[6]):
      TEMPRES=TEMPRES.derivative(z);
    for i in range(degs[7]):
      TEMPRES=TEMPRES.derivative(t);
    TEMPRES*=x^(degs[0]);
    TEMPRES*=y^(degs[1]);
    TEMPRES*=z^(degs[2]);
    TEMPRES*=t^(degs[3]);
    TEMPRES*=c^(degs[8]);
    RES+=A.monomial_coefficient(s)*TEMPRES;
  return RES;
