
import cffi
ffibuilder = cffi.FFI()

header = """
extern void get_dims_(int *nt, int *np, int *nfract);
extern void read_tables_(float *T, float *P, float *fract, float *qvmap, float *fract_map, float *qvimap, float *fract_ice_map);
"""
module = """
from my_plugin2 import ffi
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import os
with nc.Dataset("/Users/mgrecu/WDomains/MCS_OK/qvsat_map.nc") as f:
    T=f['T'][:]
    P=f['P'][:]
    fract=f['fract'][:]
    qvmap=f['qvs_map'][:]
    qv_ice_map=f['qv_ice_map'][:]
    fract_map=f['fract_map'][:]
    fract_ice_map=f['fract_ice_map'][:]
    nt=T.shape[0]
    npres=P.shape[0]
    nfract=fract.shape[-1]

fract=fract.astype('f4')
print(fract)

@ffi.def_extern()
def get_dims_(nt1,np1,nf1):
    nt1[0]=T.shape[0]
    np1[0]=P.shape[0]
    nf1[0]=nfract
    return

@ffi.def_extern()
def read_tables_(T1,P1,fract1,qvmap1,fract_map1,qvimap1,fract_ice_map1):
    x=1
    Tl=np.frombuffer(ffi.buffer(T1,nt*4),np.dtype('f4')).reshape(nt)
    Pl=np.frombuffer(ffi.buffer(P1,npres*4),np.dtype('f4')).reshape(npres)
    fractl=np.frombuffer(ffi.buffer(fract1,nfract*4),np.dtype('f4')).reshape(nfract)
    qvmapl=np.frombuffer(ffi.buffer(qvmap1,nt*npres*4),np.dtype('f4')).reshape(npres,nt)
    qvimapl=np.frombuffer(ffi.buffer(qvimap1,nt*npres*4),np.dtype('f4')).reshape(npres,nt)
    fract_mapl=np.frombuffer(ffi.buffer(fract_map1,nt*npres*nfract*4),np.dtype('f4')).reshape(nfract,npres,nt)
    fract_ice_mapl=np.frombuffer(ffi.buffer(fract_ice_map1,nt*npres*nfract*4),np.dtype('f4')).reshape(nfract,npres,nt)
    Tl[:]=T[:].data
    Pl[:]=P[:].data
    fractl[:]=fract[:].data
    print(fractl)
    qvmapl[:]=qvmap[:].data.T
    qvimapl[:]=qv_ice_map[:].data.T
    for k in range(nfract):
        fract_mapl[k,:,:]=fract_map[:,:,k].data.T
        fract_ice_mapl[k,:,:]=fract_ice_map[:,:,k].data.T
    #plt.pcolormesh(Tl,Pl,fract_ice_mapl[10,:,:],cmap='jet')
    #plt.show()
    return

"""
with open("plugin2.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("my_plugin2", r'''
    #include "plugin2.h"
''')

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target="lib_ml_interface.a", verbose=True)
