import pickle
import numpy as np
inputScalers=pickle.load(open("inputScalers.pkl","rb"))
outputScalers=pickle.load(open("outputScalers.pkl","rb"))
import sys
sys.path.append("/Users/mgrecu/WDomains/MCS_OK")
from fscale import fscale
inputVars=["qc","qr","qi","qs","qg","qns","qng","qnr","qni","qv","press","temp","dz"]
qcm=inputScalers["qc"].mean_
qcs=inputScalers["qc"].scale_
qrm=inputScalers["qr"].mean_
qrs=inputScalers["qr"].scale_
qim=inputScalers["qi"].mean_
qis=inputScalers["qi"].scale_
qsm=inputScalers["qs"].mean_
qss=inputScalers["qs"].scale_
qgm=inputScalers["qg"].mean_
qgs=inputScalers["qg"].scale_
qnsm=inputScalers["qns"].mean_
qnss=inputScalers["qns"].scale_
qngm=inputScalers["qng"].mean_
qngs=inputScalers["qng"].scale_
qnrm=inputScalers["qnr"].mean_
qnrs=inputScalers["qnr"].scale_
qnim=inputScalers["qni"].mean_
qnis=inputScalers["qni"].scale_
qvm=inputScalers["qv"].mean_
qvs=inputScalers["qv"].scale_
pressm=inputScalers["press"].mean_
presss=inputScalers["press"].scale_
tempm=inputScalers["temp"].mean_
temps=inputScalers["temp"].scale_
dzm=inputScalers["dz"].mean_
dzs=inputScalers["dz"].scale_

import time
def call_emulator(th,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,\
                  qg_curr,qni_curr,qns_curr,qnr_curr,qng_curr,rho,pi_phy,p,dt,dz8w,model):
    #print(model.summary())
    #print('th=',th.mean(axis=(0,2)))
    #stop
    print(th.shape)
    nx,nz,ny=th.shape
    nt=nx*ny
    xinput=np.zeros((nt,72,13),dtype=np.float32)
    it=0
    temp=th[:,:,:]*pi_phy[:,:,:]
    dtemp=temp-th*(p/1e5)**0.2857
    #rdcp=287.04/1005.7
    #print(dtemp.mean(axis=(0,2)))
    t0=time.time()
    xinput[:,:,0]=fscale(qc_curr[:,:72,:],qcm,qcs)
    xinput[:,:,1]=fscale(qr_curr[:,:72,:],qrm,qrs)
    xinput[:,:,2]=fscale(qi_curr[:,:72,:],qim,qis)
    xinput[:,:,3]=fscale(qs_curr[:,:72,:],qsm,qss)
    xinput[:,:,4]=fscale(qg_curr[:,:72,:],qgm,qgs)
    xinput[:,:,5]=fscale(qns_curr[:,:72,:],qnsm,qnss)
    xinput[:,:,6]=fscale(qng_curr[:,:72,:],qngm,qngs)
    xinput[:,:,7]=fscale(qnr_curr[:,:72,:],qnrm,qnrs)
    xinput[:,:,8]=fscale(qni_curr[:,:72,:],qnim,qnis)
    xinput[:,:,9]=fscale(qv_curr[:,:72,:],qvm,qvs)
    xinput[:,:,10]=fscale(p[:,:72,:],pressm,presss)
    xinput[:,:,11]=fscale(temp[:,:72,:],tempm,temps)
    xinput[:,:,12]=fscale(dz8w[:,:72,:],dzm,dzs)

    for it in range(-nx*ny):
        i=it//ny
        j=it%ny
        xinput[it,:,0]=(qc_curr[i,:72,j]-inputScalers["qc"].mean_)/inputScalers["qc"].scale_
        xinput[it,:,1]=(qr_curr[i,:72,j]-inputScalers["qr"].mean_)/inputScalers["qr"].scale_
        xinput[it,:,2]=(qi_curr[i,:72,j]-inputScalers["qi"].mean_)/inputScalers["qi"].scale_
        xinput[it,:,3]=(qs_curr[i,:72,j]-inputScalers["qs"].mean_)/inputScalers["qs"].scale_
        xinput[it,:,4]=(qg_curr[i,:72,j]-inputScalers["qg"].mean_)/inputScalers["qg"].scale_
        xinput[it,:,5]=(qns_curr[i,:72,j]-inputScalers["qns"].mean_)/inputScalers["qns"].scale_
        xinput[it,:,6]=(qng_curr[i,:72,j]-inputScalers["qng"].mean_)/inputScalers["qng"].scale_
        xinput[it,:,7]=(qnr_curr[i,:72,j]-inputScalers["qnr"].mean_)/inputScalers["qnr"].scale_
        xinput[it,:,8]=(qni_curr[i,:72,j]-inputScalers["qni"].mean_)/inputScalers["qni"].scale_
        xinput[it,:,9]=(qv_curr[i,:72,j]-inputScalers["qv"].mean_)/inputScalers["qv"].scale_
        xinput[it,:,10]=(p[i,:72,j]-inputScalers["press"].mean_)/inputScalers["press"].scale_
        xinput[it,:,11]=(temp[i,:72,j]-inputScalers["temp"].mean_)/inputScalers["temp"].scale_
        xinput[it,:,12]=(dz8w[i,:72,j]-inputScalers["dz"].mean_)/inputScalers["dz"].scale_
            
    t1=time.time()
    ypred=model.predict(xinput)
    a=np.nonzero(ypred!=ypred)
    print('nans=',len(a[0]))
    t2=time.time()
    print("time taken for prediction=",t1-t0,t2-t1)
    return 0
