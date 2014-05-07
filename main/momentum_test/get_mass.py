#!/usr/bin/python

import os, math, random, struct
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import numpy as np

bootstrapsize = 1000
T = 48 # temporal extend of the lattice
start_cfg = 750
delta_cfg = 8
nb_cfg = 101

################################################################################
# extract a correlation function and resorts it
def extract_corr_fct(name='', start_cfg=0, delta_cfg=0, nb_cfg=0, T=0, gamma=5): 
  re = []
  im = []
  for x in range(start_cfg, start_cfg+delta_cfg*nb_cfg, delta_cfg):
    filename = name + "%04d" % x + '.dat'
    f = open(filename, "rb") # Open a file
    print "reading from file: ", f.name
    f.seek(2*8*T*gamma)
    for t in range(0, T):
      re.insert(t, struct.unpack('d', f.read(8))) # returns a tuple -> convers.
      im.insert(t, struct.unpack('d', f.read(8))) # returns a tuple -> convers.
    f.close(); # close the file
  
  # conversion of the tuple to list and reorganise
  corr = [complex(0.0, 0.0)]*nb_cfg*T
  t = 0
  for x, y in zip(re, im):
    corr[(t%T)*nb_cfg + t/T] = complex(x[0], y[0])
    t += 1

  return corr
################################################################################
# writing the correlatin function in new order
def write_corr_fct(name, re_im, corr):
  f = open(name, "wb") # Open a file
  print "writing to file: ", f.name
  if re_im == 'real':
    for x in corr:
      f.write(struct.pack('d', x.real))
  elif re_im == 'imag':
    for x in corr:
      f.write(struct.pack('d', x.imag))
  else:
    print 'wring re_im -> must be real or imag'
  f.close()
################################################################################

# simple mean value
def mean(X):
  return sum(X)/ float(len(X))

# simple standard error
def std_error(X):
  m = mean(X)
  s = 0.0 
  for x in X:
    s += math.sqrt((m - x)*(m - x)) 
  return s/float(len(X)-1)

# bootstrapping
def bootstrap(X):
  boot = []
  random.seed(1227)
  for i in range(0, bootstrapsize):
    boot_dummy = 0.0 
    for j in range(0, len(X)):
      boot_dummy += X[random.randint(0, len(X)-1)]  
    boot.append(boot_dummy/len(X))
  return boot

# symmetrising and bootstrapping
def sym_and_boot(X):
  boot = []
  for t in range(0, T/2):
    data = []
    if t > 0:
      for a, b in zip(X[t*nb_cfg:(t+1)*nb_cfg], X[(T-t)*nb_cfg:(T-t+1)*nb_cfg]):
        data.append( (a+b)/2.0 )
    else:
      data = X[t*nb_cfg:(t+1)*nb_cfg]
    boot.extend(bootstrap(data))
  return boot

# mass computation
def compute_mass(boot):
  mass = []
  sigma = []
  val = []
  for t in range(1, T/2-1):
    for b in range(0, bootstrapsize):
      a = (boot[(t-1)*bootstrapsize + b] + boot[(t+1)*bootstrapsize + b])/(2.0*boot[t*bootstrapsize + b])
      if a >= 1:
        mass.append( math.acosh(a) )
      else:
        mass.append(100)
    sigma.append(std_error(mass[(t-1)*bootstrapsize:t*bootstrapsize])) 
    val.append(mean(mass[(t-1)*bootstrapsize:(t)*bootstrapsize]))
    print t, val[-1], sigma[-1]
  return mass, sigma, val

# compute the mean correlator with error
def return_mean_corr(boot):
  for t in range(1, T/2-1):
    print t, mean(boot[(t-1)*bootstrapsize:t*bootstrapsize]), std_error(boot[(t-1)*bootstrapsize:t*bootstrapsize])


################################################################################
# extracting the mass
name = '../4ptFunction_A80/data/C2_pi+-_conf'
corr = extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
write_corr_fct('./raw_data/pion_mom_0.dat', 'real', corr)
pi = []
for c in corr:
  pi.append(c.real)
print '\ncomputing the effective mass of the pion:\n'
boot = sym_and_boot(pi)
#pi2 = []
#for a in boot:
#  pi2.append(a)
mass, sigma, val = compute_mass(boot)


################################################################################
################ computing the fit #############################################
# fit function - just a constant for the mass
def fit_func(t, A):
  return A
# trying different fit intervals
for lo in range(8, 12):
  for up in range(20, 23):

    print "\nResult from bootstrapsample fit in interval:"
    print lo, up

    tt = []
    for t in range(0, T-1):
    	tt.append(t)
    x = np.asarray(tt[lo:up], float)
    
    error = np.asarray(sigma[lo:up], float)
    
    A = []
    # performing fit for each bootstrapsample
    for b in range(0, bootstrapsize):
      mass2 = []
      for t in range(lo, up):
        mass2.append(mass[t*bootstrapsize + b])
      y = np.asarray(mass2, float)
      popt, pcov = curve_fit(fit_func, x, y, p0 = [0.19], sigma=error)
      A.append(popt[0])
    print mean(A), std_error(A)

    # fitting the mean value again to extract the chi^2 value
    y = np.asarray(val[lo:up], float)
    fitfunc = lambda p,x: p[0] 
    errfunc = lambda p, x, y, error: (y-fitfunc(p,x))/error
    p,cov,infodict,mesg,ier = leastsq(errfunc, [0.19], args=(x, y, error), full_output=1)
    chisq=sum(infodict['fvec']*infodict['fvec'])

    print "Chi^2/d.o.f.:"
    print chisq/(up-lo-2)




