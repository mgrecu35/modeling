import coal_bott as cb
rq0_in=10.
xmw_in=1.0
nbins=400
dt=1.0
g_out,r_out,dlnr_out = cb.coad1d_init(rq0_in,xmw_in,nbins,dt)
t=0.0
cb.set_g_initial(2*g_out)
while t<3600:
    g,t = cb.integrate(nbins,dt,t)

import matplotlib.pyplot as plt

plt.plot(r_out[:200],g_out[:200])
plt.xscale('log')
plt.yscale('log')

plt.plot(r_out[:200],g[:200])
