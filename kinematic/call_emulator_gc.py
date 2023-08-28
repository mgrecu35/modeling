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


import time
def call_emulator(th,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,\
                  qg_curr,qni_curr,qns_curr,qnr_curr,qng_curr,rho,pi_phy,p,dt,dz8w,model,ij):
    #print(model.summary())
    #print('th=',th.mean(axis=(0,2)))
    #stop
    #print(th.shape)
    nx,nz,ny=th.shape
    nt=nx*ny
    xinput=np.zeros((nt,72,9),dtype=np.float32)
    it=0
    temp=th[:,:,:]*pi_phy[:,:,:]
    dtemp=temp-th*(p/1e5)**0.2857
    #rdcp=287.04/1005.7
    #print(dtemp.mean(axis=(0,2)))
    t0=time.time()
    #print(xinput.shape)
    xinput[:,:,0]=fscale1d(qc_curr[:,:72,:],qcm,qcs)
    xinput[:,:,1]=fscale1d(qr_curr[:,:72,:],qrm,qrs)
    xinput[:,:,2]=fscale1d(qi_curr[:,:72,:],qim,qis)
    xinput[:,:,3]=fscale1d(qs_curr[:,:72,:],qsm,qss)
    xinput[:,:,4]=fscale1d(qg_curr[:,:72,:],qgm,qgs)
    xinput[:,:,5]=fscale1d(qv_curr[:,:72,:],qvm,qvs)
    xinput[:,:,6]=fscale1d(p[:,:72,:],pressm,presss)
    xinput[:,:,7]=fscale1d(temp[:,:72,:],tempm,temps)
    xinput[:,:,8]=fscale1d(dz8w[:,:72,:],dzm,dzs)

    
            
    t1=time.time()
    ypred=model.predict(xinput,verbose=0)
    dt0=dt
    qvmin=qv_curr.min()
    qc_curr[:,:72,:]=update_var1d(qc_curr[:,:72,:],ypred[:,:,0],dt0,qc_tendm,qc_tends)
    qr_curr[:,:72,:]=update_var1d(qr_curr[:,:72,:],ypred[:,:,1],dt0,qr_tendm,qr_tends)
    qi_curr[:,:72,:]=update_var1d(qi_curr[:,:72,:],ypred[:,:,2],dt0,qi_tendm,qi_tends)
    qs_curr[:,:72,:]=update_var1d(qs_curr[:,:72,:],ypred[:,:,3],dt0,qs_tendm,qs_tends)
    qg_curr[:,:72,:]=update_var1d(qg_curr[:,:72,:],ypred[:,:,4],dt0,qg_tendm,qg_tends)
    #qns_curr[:,:72,:]=update_var(qns_curr[:,:72,:],ypred[:,:,5],dt0,qns_tendm,qns_tends)
    #qng_curr[:,:72,:]=update_var(qng_curr[:,:72,:],ypred[:,:,6],dt0,qng_tendm,qng_tends)
    #qnr_curr[:,:72,:]=update_var(qnr_curr[:,:72,:],ypred[:,:,7],dt0,qnr_tendm,qnr_tends)
    #qni_curr[:,:72,:]=update_var(qni_curr[:,:72,:],ypred[:,:,8],dt0,qni_tendm,qni_tends)
    qv_curr=update_var1d(qv_curr[:,:72,:],ypred[:,:,5],dt0,qv_tendm,qv_tends)
    temp_updated=update_var1d(temp[:,:72,:],ypred[:,:,6],dt0,temp_tendm,temp_tends)
    #print('temp_tend_max_min',ypred[:,6].max(),ypred[:,6].min(),'ij=',ij)
    #ind=np.argmax(ypred[:,10])
    qv_curr[qv_curr<qvmin]=qvmin
    qr_curr[qr_curr<0]=0
    qi_curr[qi_curr<0]=0
    qs_curr[qs_curr<0]=0
    qg_curr[qg_curr<0]=0
    th[:,:72,:]=temp_updated[:,:72,:]/pi_phy[:,:72,:]
    
    a=np.nonzero(ypred!=ypred)
    #print()
    t2=time.time()
    print('nans=',len(a[0]),"time taken for prediction=",t1-t0,t2-t1,"dt=",dt)
    return 1
