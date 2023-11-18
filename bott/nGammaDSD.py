import numpy as np
def dm_lwc(nw,lwc,rho):
    dm=(lwc*1e-3*4**4/(nw*np.pi*rho))**(0.25)
    return dm

from scipy.special import gamma as gam

def fmu(mu):
    return 6/4**4*(4+mu)**(mu+4)/gam(mu+4)

def getMassDist(nw,dm,f_mu,mu,rho):
    dD=0.02
    lambd=(4+mu)/dm
    d=np.arange(500)*dD+dD/2
    md=nw*f_mu*(d/dm)**mu*np.exp(-lambd*d)*(0.1*d)**3/6\
          *np.pi*dD/10*rho*1e3
    lwc=np.sum(md)
    return lwc,md,d