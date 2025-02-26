import numpy as np
import gvar as gv
import lsqfit

E      = [ 0.39260, 0.38940, 0.38649]
DeltaE = [ 0.00140, 0.00140, 0.00051]
L      = [      14,      16,      24]

E      = [ 0.416, 0.39260, 0.38940, 0.38649]
DeltaE = [ 0.019, 0.00140, 0.00140, 0.00051]
L      = [    12,      14,      16,      24]


y = {'data': gv.gvar(E,DeltaE)}
x = {'data': np.array(L)}

def fcn(x, p):                        # fit function of x and parameters p
  ans = {}
  ans['data'] = p['m']*(1 + p['A']*gv.exp( - p['m'] * x['data'] )/(p['m'] * x['data'])**(3/2))
  return ans

prior = {}
prior['A'] = gv.gvar(10, 100)
prior['m'] = gv.gvar(0.38, 0.05)
fit = lsqfit.nonlinear_fit(data=(x, y), prior=prior, fcn=fcn)
print(fit.format(maxline=True)) 