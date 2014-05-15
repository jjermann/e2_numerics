# The Fourier expansion of f^(-1)
# -------------------------------
#
# Usage example:
#
# fourier_coeff(1)
# 
# This gives the first non-constant Fourier coefficient of F=f^(-1)


# Returns the rational Fourier series H with default_prec=num
def getH(num=100):
  R.<q>  = PowerSeriesRing(QQ,"q",default_prec=num)
  E2     = 1-24*sum([i*q^i/(1-q^i) for i in (1..num+1)])
  E4     = 1+240*sum([i^3*q^i/(1-q^i) for i in (1..num+1)])
  E6     = 1-504*sum([i^5*q^i/(1-q^i) for i in (1..num+1)])

  H      = -12 + log((exp(12/E2-12)*q).reversion()/q)

  return H


FH=getH();

# Returns the n'th Fourier coefficient (up to 100, see above)
def fourier_coeff(n,H=FH):
  return H[n]/(2*pi*I*exp(12*n));

# Returns F=f^(-1) based on its Fourier series (with +-high order) with numerical precision prec
def getF(prec=2000,H=FH):
  q=H.parent().gen();
  return lambda t: t+H.polynomial().subs(q=exp(2*pi*I*(t)-12).n(prec))/(2*pi*I).n(prec);

# Returns the n'th order Fourier series approximation of F(t)-t
def num_fourier(n,prec=2000,H=FH):
  q=H.parent().gen();
  return lambda t: (H+O(q^(n+1))).polynomial().subs(q=exp(2*pi*I*(t)-12).n(prec))/(2*pi*I).n(prec)
