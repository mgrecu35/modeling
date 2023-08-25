<%
argnames=''
c_argnames=''
for var in varList:
    argnames+='cp_'+var+','
    c_argnames+='float *'+var+','
argnames+='nx,ny,nz'
c_argnames+='int *nx,int *ny,int *nz'                
%>

import cffi
ffibuilder = cffi.FFI()

header = """
extern void emulator_interface_(${c_argnames});
"""
module = """
from my_plugin import ffi
import numpy as np
import matplotlib.pyplot as plt


@ffi.def_extern()
def emulator_interface_(${argnames}):
    %for dim in dimens:
    ${dim}_py=${dim}[0]
    %endfor
    %for var in varList:
    %if var!='dt':
    ${var} = np.frombuffer(ffi.buffer(cp_${var},nz_py*ny_py*nx_py*4),np.dtype('f4')).reshape(nx_py,ny_py,nz_py)
    %endif
    %endfor
"""
with open("plugin.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("my_plugin", r'''
    #include "plugin.h"
''')

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target="lib_interface.a", verbose=True)
