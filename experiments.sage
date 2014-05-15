# Experiments
# -----------
#


load derivative.sage
load e2_zeros.sage 
load fast_eisenstein.sage
load fourier_numerics.sage
load fourier.sage 
load taylor_numerics.sage
load plot_results.sage


F=getF();
Z0=(-1/F(0)).imag();
def boundary_curve(x):
  mx=abs(mpmathify(x)%1);
  return mpmath.findroot(lambda y: (FE4(mx+I*y)*conjugate(FE2(mx+I*y))^2).real, (Z0+2*(sqrt(3)/2-0.000001-Z0)*mx,Z0+2*(1.2-Z0)*mx),solver='anderson')

bfourier_fun=lambda n: lambda x: boundary_curve(x)*mpmath.exp(-2*pi*I*n*x);
bcoeff=lambda n: mpmath.quad(bfourier_fun(n),[-1/2,1/2])

