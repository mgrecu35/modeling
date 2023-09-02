import pickle
import numpy as np
inputScalers=pickle.load(open("inputScalers.pkl","rb"))
outputScalers=pickle.load(open("outputScalers.pkl","rb"))
import sys
sys.path.append("/Users/mgrecu/WDomains/MCS_OK")
from fscale import fscale,update_var
inputVars=["qc","qr","qi","qs","qg","qns","qng","qnr","qni","qv","press","temp","dz"]
outputVars=["qc_tend","qr_tend","qi_tend","qs_tend","qg_tend","qns_tend","qng_tend","qnr_tend","qni_tend","qv_tend","temp_tend"]
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

qss[qss<1e-4]=1e-4
qvs[qvs<1e-4]=1e-4
qcs[qcs<1e-4]=1e-4
qrs[qrs<1e-4]=1e-4
qis[qis<1e-4]=1e-4
qgs[qgs<1e-4]=1e-4
qnss[qnss<1]=1
qngs[qngs<1]=1
qnrs[qnrs<1]=1
qnis[qnis<1]=1


qc_tendm=outputScalers["qc_tend"].mean_
qc_tends=outputScalers["qc_tend"].scale_
qr_tendm=outputScalers["qr_tend"].mean_
qr_tends=outputScalers["qr_tend"].scale_
qi_tendm=outputScalers["qi_tend"].mean_
qi_tends=outputScalers["qi_tend"].scale_
qs_tendm=outputScalers["qs_tend"].mean_
qs_tends=outputScalers["qs_tend"].scale_
qg_tendm=outputScalers["qg_tend"].mean_
qg_tends=outputScalers["qg_tend"].scale_
qns_tendm=outputScalers["qns_tend"].mean_
qns_tends=outputScalers["qns_tend"].scale_
qng_tendm=outputScalers["qng_tend"].mean_
qng_tends=outputScalers["qng_tend"].scale_
qnr_tendm=outputScalers["qnr_tend"].mean_
qnr_tends=outputScalers["qnr_tend"].scale_
qni_tendm=outputScalers["qni_tend"].mean_
qni_tends=outputScalers["qni_tend"].scale_
qv_tendm=outputScalers["qv_tend"].mean_
qv_tends=outputScalers["qv_tend"].scale_
temp_tendm=outputScalers["temp_tend"].mean_
temp_tends=outputScalers["temp_tend"].scale_


import time
def call_emulator(th,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,\
                  qg_curr,qni_curr,qns_curr,qnr_curr,qng_curr,rho,pi_phy,p,dt,dz8w,model,ij):
    #print(model.summary())
    #print('th=',th.mean(axis=(0,2)))
    #stop
    #print(th.shape)
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
    dt0=dt
    qvmin=qv_curr.min()
    #qc_curr[:,:72,:]=update_var(qc_curr[:,:72,:],ypred[:,:,0],dt0,qc_tendm,qc_tends)
    #qr_curr[:,:72,:]=update_var(qr_curr[:,:72,:],ypred[:,:,1],dt0,qr_tendm,qr_tends)
    #qi_curr[:,:72,:]=update_var(qi_curr[:,:72,:],ypred[:,:,2],dt0,qi_tendm,qi_tends)
    #qs_curr[:,:72,:]=update_var(qs_curr[:,:72,:],ypred[:,:,3],dt0,qs_tendm,qs_tends)
    #qg_curr[:,:72,:]=update_var(qg_curr[:,:72,:],ypred[:,:,4],dt0,qg_tendm,qg_tends)
    #qns_curr[:,:72,:]=update_var(qns_curr[:,:72,:],ypred[:,:,5],dt0,qns_tendm,qns_tends)
    #qng_curr[:,:72,:]=update_var(qng_curr[:,:72,:],ypred[:,:,6],dt0,qng_tendm,qng_tends)
    #qnr_curr[:,:72,:]=update_var(qnr_curr[:,:72,:],ypred[:,:,7],dt0,qnr_tendm,qnr_tends)
    #qni_curr[:,:72,:]=update_var(qni_curr[:,:72,:],ypred[:,:,8],dt0,qni_tendm,qni_tends)
    qv_curr=update_var(qv_curr[:,:72,:],ypred[:,:,9],dt0,qv_tendm,qv_tends)
    #print('qv=',(qv_curr_updated-qv_curr[:,:72,:]).max(axis=(0,2)))
    temp_updated=update_var(temp[:,:72,:],ypred[:,:,10],dt0,temp_tendm,temp_tends)
    print('temp_tend_max_min',ypred[:,:,10].max(),ypred[:,:,10].min(),'ij=',ij)
    ind=np.argmax(ypred[:,:,10])
    i0=ind//(72)
    if ypred[i0,:,10].max()>1e5:
        print('ypred=',ypred[i0,:,10].max(),'ij=',ij)
        ix=i0//ny
        iy=i0%ny
        xL=[qc_curr[ix,:,iy],qr_curr[ix,:,iy],qi_curr[ix,:,iy],qs_curr[ix,:,iy],qg_curr[ix,:,iy],qns_curr[ix,:,iy],\
            qng_curr[ix,:,iy],qnr_curr[ix,:,iy],qni_curr[ix,:,iy],qv_curr[ix,:,iy],p[ix,:,iy],temp[ix,:,iy],dz8w[ix,:,iy]]
        pickle.dump({"xL":xL,"ix":ix,"jy":iy,"ypred":ypred[i0,:],"xinput":xinput[i0,:,:]},open("xL_%i.pkl"%ij,"wb"))
    #print('temp=',(temp_updated-temp[:,:72,:]).max(axis=(0,2)),'ij=',ij)
    th[:,:72,:]=temp_updated[:,:72,:]/pi_phy[:,:72,:]
    qc_curr[qc_curr<0]=0
    qr_curr[qr_curr<0]=0
    qi_curr[qi_curr<0]=0
    qs_curr[qs_curr<0]=0
    qg_curr[qg_curr<0]=0
    qns_curr[qns_curr<0]=0
    qng_curr[qng_curr<0]=0
    qnr_curr[qnr_curr<0]=0
    qni_curr[qni_curr<0]=0
    qv_curr[qv_curr<qvmin]=qvmin
    a=np.nonzero(ypred!=ypred)
    print('nans=',len(a[0]))
    t2=time.time()
    print("time taken for prediction=",t1-t0,t2-t1,"dt=",dt)
    return 1
