dt1=10.0
import numpy as np
import matplotlib.pyplot as plt
import kinDriver as kd
import netCDF4 as nc
import wrf
import pickle

import glob
fs=sorted(glob.glob("/Users/mgrecu/WDomains/MCS_OK/save/wrfout_d03_2018-06-2*"))
print(fs)
#stop
def readWRF_File_Morr(f):
    with nc.Dataset(f) as wfile:
        xlong=wfile.variables['XLONG'][0,:,:]
        xlat=wfile.variables['XLAT'][0,:,:]
        press=wfile.variables['P'][:,:,:,:]+wfile.variables['PB'][:,:,:,:]
        theta=wfile.variables['T'][:,:,:,:]+300
        qv=wfile.variables['QVAPOR'][:,:,:,:]
        qr=wfile.variables['QRAIN'][:,:,:,:]
        qs=wfile.variables['QSNOW'][:,:,:,:]
        qg=wfile.variables['QGRAUP'][:,:,:,:]
        qc=wfile.variables['QCLOUD'][:,:,:,:]
        qi=wfile.variables['QICE'][:,:,:,:]
        qnice=wfile.variables['QNICE'][:,:,:,:]
        qnr=wfile.variables['QNRAIN'][:,:,:,:]
        qns=wfile.variables['QNSNOW'][:,:,:,:]
        qng=wfile.variables['QNGRAUPEL'][:,:,:,:]
        u=0.5*(wfile.variables['U'][:,:,:,1:]+wfile.variables['U'][:,:,:,0:-1])
        v=0.5*(wfile.variables['V'][:,:,1:,:]+wfile.variables['V'][:,:,0:-1,:])
        dbz=wrf.getvar(wfile,'dbz',timeidx=wrf.ALL_TIMES)
        height=wfile.variables['PH'][:,:,:,:]+wfile.variables['PHB'][:,:,:,:]
        height=height/9.81
        height_surf=wfile.variables['HGT'][0,:,:]

    return xlong,xlat,press,theta,qv,qs,qg,qr,qc,qi,u,v,height_surf,dbz,height,qns,qng,qnr,qnice


n1,n2=0,201
xtL=[]
ytL=[]
qcLc=[]
qrLc=[]
qiLc=[]
qsLc=[]
qgLc=[]
qnsLc=[]
qngLc=[]
qnrLc=[]
qniLc=[]
qvLc=[]
qvsLc=[]
pressLc=[]
tempLc=[]
dzLc=[]
qc_tendLc=[]
qr_tendLc=[]
qi_tendLc=[]
qs_tendLc=[]
qg_tendLc=[]
qns_tendLc=[]
qng_tendLc=[]
qnr_tendLc=[]
qni_tendLc=[]
qv_tendLc=[]
xL=[]
#press_tendL=[]
temp_tendLc=[]


fract=0.3
it=0
ncond=0
ntot=0
noprecip=0
def processFile(wrf_data,n1,n2,it,ncond,noprecip,ntot):
    xlong,xlat,press,theta,qv,qs,qg,qr,qc,qi,u,v,height_surf,\
        dbz,height,qns,qng,qnr,qnice=wrf_data
    dbzm=np.ma.array(dbz,mask=(dbz<0.0))
    for j0 in range(0,189):
        qs2d=qs[it,:,j0,n1:n2]
        qr2d=qr[it,:,j0,n1:n2]
        qc2d=qc[it,:,j0,n1:n2]
        qi2d=qi[it,:,j0,n1:n2]
        qg2d=qg[it,:,j0,n1:n2]
        qns2d=qns[it,:,j0,n1:n2]
        qnr2d=qnr[it,:,j0,n1:n2]
        qni2d=qnice[it,:,j0,n1:n2]
        qng2d=qng[it,:,j0,n1:n2]
        qv2d=qv[it,:,j0,n1:n2]
        press2d=press[it,:,j0,n1:n2]
        temp2d=theta[it,:,j0,n1:n2]*(press2d/1e5)**(287/1004)
        dtemp2d_wrf=(theta[it,:,j0,n1:n2]-theta[it,:,j0,n1:n2])*(press2d/1e5)**(287/1004)
        dz2d_wrf=height[it,1:,j0,n1:n2]-height[it,:-1,j0,n1:n2]
        dbz2d=dbzm[it,:,j0,n1:n2]
        w2d=qv2d.copy()*0.0
        dt1=10
        
        qc_tend2d,qi_tend2d,qs_tend2d,qr_tend2d,\
            ni_tend2d,ns_tend2d,nr_tend2d,t_tend2d,qv_tend2d,\
            qg_tend2d,ng_tend2d,\
            qs_stend2d,qr_stend2d,qg_stend2d,\
            nr_stend2d,ns_stend2d,ng_stend2d,qvs,qvi=kd.mphys_morrison_interface_2d(temp2d[:72,:].data,
                                                            press2d[:72,:].data,dz2d_wrf[:72,:].data,\
                                                            qc2d[:72,:],qv2d[:72,:].data,qr2d[:72,:].data,\
                                                            qi2d[:72,:].data,qni2d[:72,:].data,qs2d[:72,:].data,\
                                                            qg2d[:72,:].data,qns2d[:72,:].data,\
                                                            qnr2d[:72,:].data,qng2d[:72,:].data,w2d[:72,:].data,dt1)
        nz,nxs=press2d[:72,:].shape
        #stop
        a1=np.nonzero(qv2d[:72,:].data/qvs<1.000)
        #b1=np.nonzero()
        b1=np.nonzero((qr2d[:72,:]+qs2d[:72,:]+qg2d[:72,:]+qi2d[:72,:]+qc2d[:72,:])[a1]<1e-5)
        noprecip+=len(b1[0])
        #ntot+=nz*nxs
        # xvars='qr,qc,qi,qs,qg,press,temp,qv,qvs,qvi'
        for i1 in range(nz):
            for j1 in range(nxs):
                ntot+=1
                if qv2d[i1,j1]/qvs[i1,j1]>1.0001 or \
                    (qr2d[i1,j1]+qs2d[i1,j1]+qg2d[i1,j1]+qi2d[i1,j1]+qc2d[i1,j1])>1e-6 or \
                        qv2d[i1,j1]>qvi[i1,j1]:
                    qr_tendLc.append(qr_tend2d[i1,j1])
                    qc_tendLc.append(qc_tend2d[i1,j1])
                    qi_tendLc.append(qi_tend2d[i1,j1])
                    qg_tendLc.append(qg_tend2d[i1,j1])
                    qs_tendLc.append(qs_tend2d[i1,j1])
                    qv_tendLc.append(qv_tend2d[i1,j1])
                    xL.append([qr2d[i1,j1],qc2d[i1,j1],qi2d[i1,j1],qs2d[i1,j1],qg2d[i1,j1],\
                               press2d[i1,j1],temp2d[i1,j1],qv2d[i1,j1],qvs[i1,j1],qvi[i1,j1]])
                
                    ncond+=1
        noprecip=ntot-ncond
        continue
        qv_tendLc.extend(qv_tend2d[a1][b1][c1])
        qvLc.extend(qv2d[a1][b1][c1])
        qvsLc.extend(qvs[a1][b1][c1])
        qr_tendLc.extend(qr_tend2d[a1][b1][c1])
        qc_tendLc.extend(qc_tend2d[a1][b1][c1])
        qi_tendLc.extend(qi_tend2d[a1][b1][c1])
        qg_tendLc.extend(qg_tend2d[a1][b1][c1])
        qs_tendLc.extend(qs_tend2d[a1][b1][c1])
        ncond+=len(a1[0])
        ntot+=len(qv2d[:72,:].data.flatten())
        a1=np.nonzero(qc2d[:72,:].data+qr2d[:72,:].data+qi2d[:72,:].data+qs2d[:72,:].data+\
                      qg2d[:72,:].data<1e-8)
        noprecip+=len(a1[0])
        
    print('fracts=',ncond/ntot,noprecip/ntot)
    return ncond,noprecip,ntot
    
for fname in fs[0:1]:
    wrfdata=readWRF_File_Morr(fname)
    print(fname)
    for it in range(12):
        print('itime=',it)
        ncond,noprecip,ntot=processFile(wrfdata,n1,n2,it,ncond,noprecip,ntot)
        #stop


qv_tendLc=np.array(qv_tendLc)
qr_tendLc=np.array(qr_tendLc)
qc_tendLc=np.array(qc_tendLc)
qi_tendLc=np.array(qi_tendLc)
qg_tendLc=np.array(qg_tendLc)
qs_tendLc=np.array(qs_tendLc)
qvsLc=np.array(qvsLc)
qvLc=np.array(qvLc)
#print(np.corrcoef(qv_tendLc,qr_tendLc+qc_tendLc+qi_tendLc+qs_tendLc+qg_tendLc))

#ax=plt.subplot(111)
#plt.scatter(qc_tendLc+qr_tendLc,qv_tendLc)
#plt.xlabel('qc+qr tendencies [g/g/s]]')
#plt.ylabel('qv tendency [g/g/s]')
#ax.set_aspect(1.0)
#plt.title('Tendencies at saturation')
#plt.savefig('tendencies_at_saturation_2.png')
xL=np.array(xL)
import xarray as xr
ds=xr.Dataset({'qr_tend':(['n'],qr_tendLc),'qc_tend':(['n'],qc_tendLc),'qi_tend':(['n'],qi_tendLc),\
                'qg_tend':(['n'],qg_tendLc),'qs_tend':(['n'],qs_tendLc),\
                'qv_tend':(['n'],qv_tendLc),'xL':(['n','m'],xL)})

ds.to_netcdf('tendencies_training_data.nc',\
             encoding={'xL':{'zlib':True,'complevel':5},\
                    'qr_tend':{'zlib':True,'complevel':5},\
                    'qc_tend':{'zlib':True,'complevel':5},\
                    'qi_tend':{'zlib':True,'complevel':5},\
                    'qg_tend':{'zlib':True,'complevel':5},\
                    'qs_tend':{'zlib':True,'complevel':5},\
                    'qv_tend':{'zlib':True,'complevel':5}})
