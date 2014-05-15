# Plotting functions
# ------------------
#
# Usage example:
#
# getFD_plot(): Returns a basic plot of the fundamental domain.


from sage.plot.circle import circle;

def plot_fun(f):
  a=lambda u: f(u).real()
  b=lambda u: f(u).imag()
  return (a,b)
def mpplot_fun(f):
  a=lambda u: mpmath.mpmath_to_sage(f(mpmathify(u)).real,mp.prec)
  b=lambda u: mpmath.mpmath_to_sage(f(mpmathify(u)).imag,mp.prec)
  return (a,b)

def getFD_plot():
  PA = plot((cos,sin),(0,pi),parametric=True,aspect_ratio=1);
  PB = plot((-1/2,lambda t: t),(sqrt(3)/2,2),parametric=True,aspect_ratio=True);
  PC = plot((1/2,lambda t: t),(sqrt(3)/2,2),parametric=True,aspect_ratio=True);
  PD = plot((lambda t: t,0),(-1,1),parametric=True,aspect_ratio=True);
  PE = plot((0,lambda t: t),(0,2),parametric=True,aspect_ratio=True);
  PF = point((-1/2,sqrt(3)/2),color="red",size=30);
  PG = point((1/2,sqrt(3)/2),color="red",size=30);
  PLOT=PA+PB+PC+PD+PE+PF+PG;
  PLOT.set_aspect_ratio(1);
  
  return PLOT;

# Returns a graphics of a (list of) complex point(s)
def getCPoints(list,psize=5,color='red'):
  return point(([((mpmath.mpmathify(P).real,mpmath.mpmathify(P).imag)) for P in list]),rgbcolor=color,size=psize,zorder=10);

# Returns the graphics of all zeros of weight n in the fundamental domain
def get_H_height(n):
  load e2_zeros.sage

  PG=Graphics();
  for a in range(1,n+1):
    if (gcd(a,n)==1):
      PG+=getCPoints([getWZero(Rational((a,n)))]);
  return PG;


# ----------------------------------------------------------------------
# Some summarized results
# ----------------------------------------------------------------------

# Plots the Taylor approximation to F with the given order and
# plots the zeros of height n, together with a basic bounding box
def plot_H_taylor_height(n=100,order=6):
  load taylor_numerics.sage
  err=0.00029;   

  PG  = plot(mpplot_fun(num_taylor_pol(order)), -1/2,1/2,parametric=True);
  PG += line(([( -1/2-err, 6/pi-err ),(  1/2+err, 6/pi-err )]));
  PG += line(([( -1/2-err, 6/pi+err ),(  1/2+err, 6/pi+err )]));
  PG += line(([( -1/2-err, 6/pi-err ),( -1/2-err, 6/pi+err )]));
  PG += line(([(  1/2+err, 6/pi-err ),(  1/2+err, 6/pi+err )]));

  PG += get_H_height(n);
  return PG;

# Plots the "corrected" zeros of height n together with the lower and upper bounds
def plot_corrected_height(n=100):
  load e2_zeros.sage
  load fourier.sage

  PG  = circle((0,6/pi),0.00027);
  PG += circle((0,6/pi),0.00029);

  for a in range(1,n+1):
    if (gcd(a,n)==1):
      PG+=getCPoints([getWcdZero(n,-a)-a/n]);

  PG.set_aspect_ratio(1);
  return PG;

# Plots the "corrected" zeros of height n together with the first order
# Fourier series approximation of f^(-1)(t)-t (which is actually already very accurate) 
def plot_corrected_fourier_height(n=100):
  load e2_zeros.sage
  load fourier.sage

  PG  = plot(plot_fun(num_fourier(1)), -1/2,1/2,parametric=True);

  for a in range(1,n+1):
    if (gcd(a,n)==1):
      PG+=getCPoints([getWcdZero(n,-a)-a/n],20);

  PG.set_aspect_ratio(1);
  return PG;

# Returns a list of the first N numerically calculated Fourier coefficients of f^(-1)(t)-t,
# set mp.dps/mp.prec accordingly for higher precision
def list_num_fourier(N=10):
  load fourier_numerics.sage
  return [num_fourier_coeff(n) for n in range(0,N+1)]

# Returns a list of the first N Taylor coefficients of f^(-1),
# set mp.dps/mp.prec accordingly for higher precision
def list_num_taylor(N=10):
  load taylor_numerics.sage
  return [num_taylor_eval(taylor_coeff(n)) for n in range(0,N+1)]

# Returns a list of the first N Taylor coefficients (starting at n=1) of f^(-1) as polynomials,
# where each coefficient is divided by (-1)^n*(12/c)^(1-n)*x^(1-2n)
def list_modified_taylor_as_pol(N=10):
  load taylor_numerics.sage
  return [taylor_coeff(n)*(-1)^n*(12/c)^(n-1)*x^(2*n-1) for n in range(1,N+1)]

# Fourier expansion over Q of H up to order n            
def fourier_H(n=7):
  load fourier.sage
  print getH(n)


# Plots a picture of (one) domain of definition of f
def getfDomain_plot():
  cplot=lambda f:(lambda t:f(t).real(),lambda t:f(t).imag());
  FD=getFD_plot();
  CPA1=contour_plot(lambda x,y:-1*(x^2+y^2-1),(x,-0.5,0.5),(y,0,2),contours=1,cmap=["lightgrey","white"]);
  CPA2=contour_plot(lambda x,y:-1*((x-1)^2+y^2-1),(x,0.5,1),(y,0,2),contours=1,cmap=["lightgrey","white"]);
  CPA3=contour_plot(lambda x,y:-1*((x+1)^2+y^2-1),(x,-1,-0.5),(y,0,2),contours=1,cmap=["lightgrey","white"]);
  CPA4=contour_plot(lambda x,y:1,(x,-1,1),(y,-1.1,2),contours=1,cmap=["lightgrey","white"]);
  CP1B = plot((lambda t: cos(t),lambda t: sin(t)),(pi/3,2*pi/3),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP2B = plot((lambda t: cos(t)+1,lambda t: sin(t)),(pi/2,2*pi/3),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP3B = plot((lambda t: cos(t)-1,lambda t: sin(t)),(pi/3,pi/2),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP=plot(cplot(lambda t: F(t-I*sqrt(3)/2)),(-1,1),parametric=True,color="red",thickness=1.5,zorder=10)
  CP8=plot(cplot(lambda t: F(t-(1-1/8)*I*sqrt(3)/2)),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP4=plot(cplot(lambda t: F(t-(1-1/4)*I*sqrt(3)/2)),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP2=plot(cplot(lambda t: F(t-(1-1/2)*I*sqrt(3)/2)),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP1=plot(cplot(lambda t: F(t)),(-1,1),parametric=True,color="red",thickness=1,zorder=10)
  PF = point((-1/2,sqrt(3)/2),color="red",size=40,zorder=10);
  PG = point((1/2,sqrt(3)/2),color="red",size=40,zorder=10);

  PLOT=CPA4+CPA1+CPA2+CPA3+FD+CP+CP1+CP2+CP4+CP8+CP1B+CP2B+CP3B+PF+PG;
  PLOT.set_aspect_ratio(1);
  
  return PLOT;

# Plots a picture of the corresponding range of f
def getfRange_plot():
  cplot=lambda f:(lambda t:f(t).real(),lambda t:f(t).imag());
  FD=getFD_plot();
  CPA1=contour_plot(lambda x,y:(x^2+y^2-1),(x,-0.5,0.5),(y,-1.1,0),contours=1,cmap=["lightgrey","white"]);
  CPA2=contour_plot(lambda x,y:((x-1)^2+y^2-1),(x,0.5,1),(y,-1.1,0),contours=1,cmap=["lightgrey","white"]);
  CPA3=contour_plot(lambda x,y:((x+1)^2+y^2-1),(x,-1,-0.5),(y,-1.1,0),contours=1,cmap=["lightgrey","white"]);
  CPA4=contour_plot(lambda x,y:-y,(x,-1,1),(y,-1.1,2),contours=1,cmap=["lightgrey","white"]);
  CP1B = plot((lambda t: cos(t),lambda t: -sin(t)),(pi/3,2*pi/3),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP2B = plot((lambda t: cos(t)+1,lambda t: -sin(t)),(pi/2,2*pi/3),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP3B = plot((lambda t: cos(t)-1,lambda t: -sin(t)),(pi/3,pi/2),parametric=True,aspect_ratio=1,thickness=3,color="black");
  CP=plot(cplot(lambda t: t-I*sqrt(3)/2),(-1,1),parametric=True,color="red",thickness=1.5,zorder=10)
  CP8=plot(cplot(lambda t: t-(1-1/8)*I*sqrt(3)/2),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP4=plot(cplot(lambda t: t-(1-1/4)*I*sqrt(3)/2),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP2=plot(cplot(lambda t: t-(1-1/2)*I*sqrt(3)/2),(-1,1),parametric=True,color="red",thickness=0.5,zorder=10)
  CP1=plot((lambda t: t, lambda t: 0),            (-1,1),parametric=True,color="red",thickness=1,zorder=10)
  PF = point((-1/2,-sqrt(3)/2),color="red",size=40,zorder=10);
  PG = point((1/2,-sqrt(3)/2),color="red",size=40,zorder=10);

  PLOT=CPA4+CPA1+CPA2+CPA3+FD+CP+CP1+CP2+CP4+CP8+CP1B+CP2B+CP3B+CP+PF+PG;
  PLOT.set_aspect_ratio(1);
  
  return PLOT;

