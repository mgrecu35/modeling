import pickle
import numpy as np
inputScalers=pickle.load(open("inputScalers_1M_cell.pkl","rb"))
outputScalers=pickle.load(open("outputScalers_1M_cell.pkl","rb"))
import sys
sys.path.append("/Users/mgrecu/WDomains/MCS_OK")
from fscale import fscale1,update_var1, fscale1d,update_var1d
inputVars=["qc","qr","qi","qs","qg","qns","qng","qnr","qni","qv","press","temp","dz"]
outputVars=["qc_tend","qr_tend","qi_tend","qs_tend","qg_tend","qns_tend","qng_tend","qnr_tend","qni_tend","qv_tend","temp_tend"]
qcm,qcs=inputScalers["qc"]
qrm,qrs=inputScalers["qr"]
qim,qis=inputScalers["qi"]
qsm,qss=inputScalers["qs"]
qgm,qgs=inputScalers["qg"]
qvm,qvs=inputScalers["qv"]
pressm,presss=inputScalers["press"]
tempm,temps=inputScalers["temp"]
dzm,dzs=inputScalers["dz"]



qc_tendm,qc_tends=outputScalers["qc_tend"]
qr_tendm,qr_tends=outputScalers["qr_tend"]
qi_tendm,qi_tends=outputScalers["qi_tend"]
qs_tendm,qs_tends=outputScalers["qs_tend"]
qg_tendm,qg_tends=outputScalers["qg_tend"]
qv_tendm,qv_tends=outputScalers["qv_tend"]
temp_tendm,temp_tends=outputScalers["temp_tend"]

