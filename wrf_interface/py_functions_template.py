<%
argnames=''
c_argnames=''
py_argnames=''
for var in varList:
    argnames+='cp_'+var+','
    c_argnames+='float *'+var+','
    py_argnames+=var+','
argnames+='nx,ny,nz,ireturn'
c_argnames+='int *nx,int *ny,int *nz, int *ireturn'                
%>

import cffi
ffibuilder = cffi.FFI()

header = """
extern void emulator_interface_(${c_argnames});
extern void hello_world(void);
"""
module = """
from my_plugin import ffi
import numpy as np
import matplotlib.pyplot as plt
from keras.models import load_model
model=load_model("unet1D_MP_ML_2.h5")
import sys
sys.path.append("/Users/mgrecu/WRF-4.3/")
from call_emulator import call_emulator

@ffi.def_extern()
def emulator_interface_(${argnames}):
    %for dim in dimens:
    ${dim}_py=${dim}[0]
    %endfor
    dt=cp_dt[0]
    ireturn_py=ireturn[0]
    #print('nx_py,ny_py,nz_py=',nx_py,ny_py,nz_py,'in emulator_interface_')
    %for var in varList:
    %if var!='dt':
    ${var} = np.frombuffer(ffi.buffer(cp_${var},nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(ny_py,nz_py,nx_py)
    %endif
    %endfor
    ireturn[0]=0
    return
    ireturn[0]=call_emulator(${py_argnames}model)
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
