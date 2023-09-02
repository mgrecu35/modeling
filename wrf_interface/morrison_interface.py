import cffi
ffibuilder = cffi.FFI()

header = """
extern void hello_world_(int *array, int n);
extern void morrison_interface_(void);
"""

module = """
from my_plugin import ffi
import numpy as np

@ffi.def_extern()
def hello_world_(array,n):
    print("Hello World!")
    input_array=np.frombuffer(ffi.buffer(array, n*4), np.dtype('i4'))
    result=np.sum(input_array)
    input_array[:]*=3;
    #print(input_array)
    #print(input_array.ctypes.data)
    #array=ffi.cast("int *",input_array.ctypes.data)
    print(result)

@ffi.def_extern()
def morrison_interface_():
    print("Hello Morisson!")
"""

with open("plugin.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("my_plugin", r'''
    #include "plugin.h"
''')

ffibuilder.embedding_init_code(module)
ffibuilder.compile(target="lib_interface.a", verbose=True)
