

import cffi
ffibuilder = cffi.FFI()

header = """
extern void emulator_interface_(float *th,float *qv_curr,float *qc_curr,float *qr_curr,float *qi_curr,float *qs_curr,float *qg_curr,float *qni_curr,float *qns_curr,float *qnr_curr,float *qng_curr,float *rho,float *pi_phy,float *p,float *dt,float *dz8w,int *nx,int *ny,int *nz, int *ireturn, int *ij);
extern void hello_world(void);
"""
module = """
from my_plugin import ffi
import numpy as np
import matplotlib.pyplot as plt
from keras.models import load_model
import sys
sys.path.append("/Users/mgrecu/WRF-4.3/")
model=load_model("unet1D_MP_model_1M.h5")

from call_emulator_gc import call_emulator

@ffi.def_extern()
def emulator_interface_(cp_th,cp_qv_curr,cp_qc_curr,cp_qr_curr,cp_qi_curr,cp_qs_curr,cp_qg_curr,cp_qni_curr,cp_qns_curr,cp_qnr_curr,cp_qng_curr,cp_rho,cp_pi_phy,cp_p,cp_dt,cp_dz8w,nx,ny,nz,ireturn,ij):
    #return 0
    nz_py=nz[0]
    ny_py=ny[0]
    nx_py=nx[0]
    dt=cp_dt[0]
    ireturn_py=ireturn[0]
    ij_py=ij[0]
    #print(ij)
    #print('nx_py,ny_py,nz_py=',nx_py,ny_py,nz_py,'in emulator_interface_')
    th = np.frombuffer(ffi.buffer(cp_th,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qv_curr = np.frombuffer(ffi.buffer(cp_qv_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qc_curr = np.frombuffer(ffi.buffer(cp_qc_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qr_curr = np.frombuffer(ffi.buffer(cp_qr_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qi_curr = np.frombuffer(ffi.buffer(cp_qi_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qs_curr = np.frombuffer(ffi.buffer(cp_qs_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qg_curr = np.frombuffer(ffi.buffer(cp_qg_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qni_curr = np.frombuffer(ffi.buffer(cp_qni_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qns_curr = np.frombuffer(ffi.buffer(cp_qns_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qnr_curr = np.frombuffer(ffi.buffer(cp_qnr_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    qng_curr = np.frombuffer(ffi.buffer(cp_qng_curr,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    rho = np.frombuffer(ffi.buffer(cp_rho,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    pi_phy = np.frombuffer(ffi.buffer(cp_pi_phy,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    p = np.frombuffer(ffi.buffer(cp_p,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    dz8w = np.frombuffer(ffi.buffer(cp_dz8w,nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    ireturn[0]=1
    iemulator=call_emulator(th,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,qni_curr,qns_curr,qnr_curr,qng_curr,rho,pi_phy,p,dt,dz8w,model,ij_py)
    #print('th=',th.mean(axis=(0,2)))
    #print('qv=',qv_curr.mean(axis=(0,2)))

@ffi.def_extern()
def hello_world():
    print("hello world")

"""
with open("plugin.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("my_plugin", r'''
    #include "plugin.h"
''')

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target="lib_interface.a", verbose=True)
