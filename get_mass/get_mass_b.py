#!/usr/bin/python
import os, math, random, struct, sys
import matplotlib
matplotlib.use('QT4Agg')
#import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
#from matplotlib.backends.backend_pdf import PdfPages
#pdfplot = PdfPages('./test.pdf')
#import matplotlib.mlab as mlab
from scipy.stats import chi2

bootstrapsize = 1000
T = 48 # temporal extend of the lattice
start_cfg = int(sys.argv[1])
delta_cfg = 8
nb_cfg = 20
dirac_one = 5
dirac_two = 5
disp_one = 0
disp_two = 0
sinh = True
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
      s += (m - x)*(m - x)
  return np.sqrt(s/float(len(X)-1))

# simple standard error
def std_error2(X):
  m = mean(X)
  s = 0.0
  for x in X:
      s += (m - x)*(m - x)
  return m, np.sqrt(s/float(len(X)-1))

# computation of all kind of errors
def error_computation(X):
  std_err = std_error2(X)
  X.sort()
  median = X[len(X)/2]
  length = int(len(X)*0.159)
  low_err = abs(median - X[+length])
  upp_err = abs(median - X[-length])
  print '%.7e\t%.5e\t%.7e\t%.5e\t%.5e' % \
        (std_err[0], std_err[1], median, low_err, upp_err)
  return std_err[0], std_err[1], median, low_err, upp_err

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

# bootstrap without symmetrisation
def boot_wo_sym(data):
  boot = []
  tmp_dat = []
  for t in range(1,T):
    tmp_dat = data[(t-1)*nb_cfg:t*nb_cfg]
    print t, len(tmp_dat)
    boot.extend(bootstrap(tmp_dat))
  return boot

# symmetrising and bootstrapping
def sym_and_boot(X, anti):
  print "input length: ", len(X)
  boot = []
  for t in range(0, T/2):
    data = []
    if t > 0:
      for a, b in zip(X[t*nb_cfg:(t+1)*nb_cfg], X[(T-t)*nb_cfg:(T-t+1)*nb_cfg]):
        #no symmetrisation for sinh
        if anti:
          data.append( (a-b)/2.0 )
        else:  
          data.append( (a+b)/2.0 )
    else:
      data = X[t*nb_cfg:(t+1)*nb_cfg]
    boot.extend(bootstrap(data))
  print "bootstrap length:", len(boot) 
  return boot

# mass computation
def compute_mass(boot,anti):
  mass = []
  sigma = []
  val = []
  for t in range(1, T/2-1):
    for b in range(0, bootstrapsize):
      a = (boot[(t-1)*bootstrapsize + b] + boot[(t+1)*bootstrapsize + b])/(2.0*boot[t*bootstrapsize + b])
      if anti:
          mass.append( math.asinh(a) )
      else:
        if a >= 1:
          mass.append( math.acosh(a) )
        else:
          mass.append(100)
    sigma.append(std_error(mass[(t-1)*bootstrapsize:t*bootstrapsize])) 
    val.append(mean(mass[(t-1)*bootstrapsize:(t)*bootstrapsize]))
#    print t, val[-1], sigma[-1]
  return mass, sigma, val

# compute the mean correlator with error
def return_mean_corr(boot):
  print "length is:", len(boot)
  t_tmp, val_tmp, err_tmp = [], [], []
  for t in range(1, T/2):
     t_tmp.append(t)
     val_tmp.append(mean(boot[(t-1)*bootstrapsize:t*bootstrapsize]))
     err_tmp.append(std_error(boot[(t-1)*bootstrapsize:t*bootstrapsize]))
  return t_tmp, val_tmp, err_tmp
# plot correlator
def plot_corr(a, b, c):
  plt.grid(True)
  plt.title("Correlator")
  plt.xlabel("t")
  plt.ylabel("value")
  plt.plot(a,b,'bo')
  plt.errorbar(a,b,yerr=c, linestyle="None")
  pdfplot.savefig()
  pdfplot.close()

################################################################################
# extracting the mass
name = "/hiskp2/werner/LapH/correlators/A100_test/dirac_%02d"%dirac_one+"_%02d"%dirac_two+"_p_0_0_displ_%d" %disp_one + "_%d" %disp_two + "/C2_pi+-_conf"
corr = extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
of_name = "dirac_%02d"%dirac_one+"_%02d"%dirac_two+"_displ_%d" %disp_one + "_%d" %disp_two + ".dat"
write_corr_fct(of_name, 'real', corr)
pi = []
for c in corr:
  pi.append(c.imag)
print '\ncomputing the effective mass of the pion:\n'
boot = boot_wo_sym(pi)
#boot = sym_and_boot(pi,sinh)
t_tmp, val_tmp, err_tmp = return_mean_corr(boot)
for a,b,c in zip(t_tmp, val_tmp, err_tmp):
  print a, b, c
#pi2 = []
#for a in boot:
#  pi2.append(a)
mass, sigma, val = compute_mass(boot,sinh)
################################################################################
################ computing the fit #############################################
# fit function - just a constant for the mass
def plot_results(X, mean, dof, fitparameter="none"):
  X_mean, X_stderr, X_median, X_up, X_lo = error_computation(X)
  X_sigma = 0.5*(X_up + X_lo)
  # plotting the results
  n, bins, patches = plt.hist(X, bins = 50, normed=1, histtype='stepfilled')
  plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
  # plot gaussian distributions      
  x = np.linspace(X_mean-4*X_stderr, X_mean+4*X_stderr, 1000)
  p1, = plt.plot(x, chi2.pdf(x, dof), 'r')
  p2, = plt.plot(x, chi2.pdf(x, X_median), 'b')
  p3, = plt.plot(x, chi2.pdf(x, mean), 'g')
  plt.title(fitparameter)
  plt.grid(True)
  plt.legend([p1, p2, p3], ["mean, std", r'median, $1\sigma$', 'mean fit'])
  pdfplot.savefig()
  pdfplot.close()


def fit_func(t, A):
  return A

# trying different fit intervals
fitfunc1 = lambda p,x: p[0] 
errfunc = lambda p, x, y, error: (y-fitfunc1(p,x))/error

for lo in range(12, 13):
 for up in range(22, 23):

   print "\nResult from bootstrapsample fit in interval:"
   print lo, up

   tt = []
   for t in range(0, T-1):
   	tt.append(t)
   x = np.asarray(tt[lo:up], float)

   error = np.asarray(sigma[lo:up], float)
   
   results, chisquare, chisquare_full = [], [], []
   # performing fit for each bootstrapsample
   for b in range(0, bootstrapsize):
     mass2 = []
     for t in range(lo, up):
       mass2.append(mass[t*bootstrapsize + b])
     y = np.asarray(mass2, float)
     p,cov,infodict,mesg,ier = leastsq(errfunc, [0.19], args=(x, y, error), full_output=1)
     chisquare.append(float(sum(infodict['fvec']**2.)))
     chisquare_full.extend(np.asarray((infodict['fvec']), float))
     popt, pcov = curve_fit(fit_func, x, y, p0 = [0.19], sigma=error)
     results.append(p)
   print "mass = ", mean(y)
   print "dmass = ", std_error(y)
   print "dof =", float(up-lo-1)
   print "chi2 =", mean(chisquare)

   # fitting the mean value to extract the chi^2 value
   y = np.asarray(val[lo:up], float)
   p,cov,infodict,mesg,ier = leastsq(errfunc, [0.12], args=(x, y, error), full_output=1)
   chisq = sum(infodict['fvec']**2.)

   plt.grid(True)
   x1 = np.linspace(lo, up, 1000)
   fitval = mean(results)
   y1 = []
   for i in x1:
     y1.append(fitval)
   y1 = np.asarray(y1, float)
   plt.plot(x1, y1, 'r')
   y1 = []
   for i in x1:
     y1.append(p[0])
   y1 = np.asarray(y1, float)
   plt.plot(x1, y1, 'g')

   x = np.asarray(tt[5:up], float)
   y = np.asarray(val[5:up], float)
   error = np.asarray(sigma[5:up], float)
   plt.errorbar(x, y, error, fmt='x' + 'b')
   #plt.plot(x, y, 'r')
   pdfplot.savefig()
   plt.clf()

   plt.grid(True)
   x = np.asarray(tt[5:up], float)
   y = np.asarray(val_tmp[5:up], float)
   error = np.asarray(err_tmp[5:up], float)
   plt.errorbar(x, y, error, fmt='x' + 'b')
   plt.plot(x, y, 'r')
   pdfplot.savefig()
   plt.clf()
   # plotting the results of the fit
   # plot_results(chisquare, chisq, float(up-lo-1), fitparameter="chisquare")
   # plot_results(chisquare_full, float(up-lo-1), fitparameter="asdf")
   # error_computation(results)

   print "Chi^2/d.o.f.:", chisq
   print "err =", cov

pdfplot.close()

