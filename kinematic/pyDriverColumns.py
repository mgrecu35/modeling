dt1=10.0
import numpy as np
import matplotlib.pyplot as plt
import kinDriver as kd
import wrf
import pickle
#inputData=pickle.load(open('inputData.pkl','rb'))
import netCDF4 as nc
with nc.Dataset('input_data.nc','r') as f:
     inputData=f['input_data'][:,:,:]
#input_data=[qs2d,qr2d,qc2d,qi2d,qg2d,qns2d,qnr2d,qni2d,qng2d,qnv2d,press2d,temp2d,dz2d_wrf,temp_tend2d,]
if 1==1:
     qs2d=inputData[0,:,:]
     qr2d=inputData[1,:,:]
     qc2d=inputData[2,:,:]
     qi2d=inputData[3,:,:]
     qg2d=inputData[4,:,:]
     ns2d=inputData[5,:,:]
     nr2d=inputData[6,:,:]
     ni2d=inputData[7,:,:]
     ng2d=inputData[8,:,:]
     qv2d=inputData[9,:,:]
     press2d=inputData[10,:,:]
     temp2d=inputData[11,:,:]
     dz2d_wrf=inputData[12,:,:]
     w2d=inputData[12,:,:]*0.0
     dbz2d=inputData[13,:,:]
     t_tend2d_wrf=inputData[14,:,:]

#temp2d,press2d,dz2d_wrf,qv2d,qr2d,\
#    qi2d,ni2d,qs2d,qg2d,ns2d,nr2d,ng2d,w2d,dt1=inputData

j0=150
import glob
fs=sorted(glob.glob("/Users/mgrecu/WDomains/MCS_OK/wrfout_d03_2018-06-2*"))
print(fs)
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


n1,n2=50,175
xtL=[]
ytL=[]
qcL=[]
qrL=[]
qiL=[]
qsL=[]
qgL=[]
qnsL=[]
qngL=[]
qnrL=[]
qniL=[]
qvL=[]
pressL=[]
tempL=[]
dzL=[]
qc_tendL=[]
qr_tendL=[]
qi_tendL=[]
qs_tendL=[]
qg_tendL=[]
qns_tendL=[]
qng_tendL=[]
qnr_tendL=[]
qni_tendL=[]
qv_tendL=[]
#press_tendL=[]
temp_tendL=[]


fract=0.3
it=0
def processFile(wrf_data,n1,n2,it):
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
        dt1=5.0
        dt1=10
        qc_tend2d,qi_tend2d,qs_tend2d,qr_tend2d,\
            ni_tend2d,ns_tend2d,nr_tend2d,t_tend2d,qv_tend2d,\
            qg_tend2d,ng_tend2d,\
            qs_stend2d,qr_stend2d,qg_stend2d,\
            nr_stend2d,ns_stend2d,ng_stend2d=kd.mphys_morrison_interface_2d(temp2d[:72,:].data,
                                                            press2d[:72,:].data,dz2d_wrf[:72,:].data,\
                                                            qc2d[:72,:],qv2d[:72,:].data,qr2d[:72,:].data,\
                                                            qi2d[:72,:].data,qni2d[:72,:].data,qs2d[:72,:].data,\
                                                            qg2d[:72,:].data,qns2d[:72,:].data,\
                                                            qnr2d[:72,:].data,qng2d[:72,:].data,w2d[:72,:].data,dt1)
    
        r=np.random.rand(qc_tend2d.shape[1])
        a=np.nonzero(r<fract)
        qcL.extend(qc2d[:72,a[0]].transpose())
        qrL.extend(qr2d[:72,a[0]].transpose())
        qiL.extend(qi2d[:72,a[0]].transpose())
        qsL.extend(qs2d[:72,a[0]].transpose())
        qgL.extend(qg2d[:72,a[0]].transpose())
        qnsL.extend(qns2d[:72,a[0]].transpose())
        qngL.extend(qng2d[:72,a[0]].transpose())
        qnrL.extend(qnr2d[:72,a[0]].transpose())
        qniL.extend(qni2d[:72,a[0]].transpose())
        qvL.extend(qv2d[:72,a[0]].transpose())
        pressL.extend(press2d[:72,a[0]].transpose())
        tempL.extend(temp2d[:72,a[0]].transpose())
        dzL.extend(dz2d_wrf[:72,a[0]].transpose())
        qc_tendL.extend(qc_tend2d[:72,a[0]].transpose())
        qr_tendL.extend(qr_tend2d[:72,a[0]].transpose())
        qi_tendL.extend(qi_tend2d[:72,a[0]].transpose())
        qs_tendL.extend(qs_tend2d[:72,a[0]].transpose())
        qg_tendL.extend(qg_tend2d[:72,a[0]].transpose())
        qns_tendL.extend(ns_tend2d[:72,a[0]].transpose())
        qng_tendL.extend(ng_tend2d[:72,a[0]].transpose())
        qnr_tendL.extend(nr_tend2d[:72,a[0]].transpose())
        qni_tendL.extend(ni_tend2d[:72,a[0]].transpose())
        qv_tendL.extend(qv_tend2d[:72,a[0]].transpose())
        temp_tendL.extend(t_tend2d[:72,a[0]].transpose())


     
import xarray as xr
fname=fs[5]

#print(wrfdata[0].shape)

n1,n2=0,201


for fname in fs[4:-1]:
    wrfdata=readWRF_File_Morr(fname)
    print(fname)
    for it in range(12):
        print('itime=',it)
        processFile(wrfdata,n1,n2,it)
qcXR=xr.DataArray(qcL,dims=['nt','nz'])
qrXR=xr.DataArray(qrL,dims=['nt','nz'])
qiXR=xr.DataArray(qiL,dims=['nt','nz'])
qsXR=xr.DataArray(qsL,dims=['nt','nz'])
qgXR=xr.DataArray(qgL,dims=['nt','nz'])
qnsXR=xr.DataArray(qnsL,dims=['nt','nz'])
qngXR=xr.DataArray(qngL,dims=['nt','nz'])
qnrXR=xr.DataArray(qnrL,dims=['nt','nz'])
qniXR=xr.DataArray(qniL,dims=['nt','nz'])
qvXR=xr.DataArray(qvL,dims=['nt','nz'])
pressXR=xr.DataArray(pressL,dims=['nt','nz'])
tempXR=xr.DataArray(tempL,dims=['nt','nz'])
dzXR=xr.DataArray(dzL,dims=['nt','nz'])
qc_tendXR=xr.DataArray(qc_tendL,dims=['nt','nz'])
qr_tendXR=xr.DataArray(qr_tendL,dims=['nt','nz'])
qi_tendXR=xr.DataArray(qi_tendL,dims=['nt','nz'])
qs_tendXR=xr.DataArray(qs_tendL,dims=['nt','nz'])
qg_tendXR=xr.DataArray(qg_tendL,dims=['nt','nz'])
qns_tendXR=xr.DataArray(qns_tendL,dims=['nt','nz'])
qng_tendXR=xr.DataArray(qng_tendL,dims=['nt','nz'])
qnr_tendXR=xr.DataArray(qnr_tendL,dims=['nt','nz'])
qni_tendXR=xr.DataArray(qni_tendL,dims=['nt','nz'])
qv_tendXR=xr.DataArray(qv_tendL,dims=['nt','nz'])
temp_tendXR=xr.DataArray(temp_tendL,dims=['nt','nz'])


ds=xr.Dataset({'qc':qcXR,'qr':qrXR,'qi':qiXR,'qs':qsXR,'qg':qgXR,'qns':qnsXR,'qng':qngXR,'qnr':qnrXR,'qni':qniXR,'qv':qvXR,     'press':pressXR,'temp':tempXR,'dz':dzXR,'qc_tend':qc_tendXR,'qr_tend':qr_tendXR,'qi_tend':qi_tendXR,'qs_tend':qs_tendXR,'qg_tend':qg_tendXR,'qns_tend':qns_tendXR,'qng_tend':qng_tendXR,'qnr_tend':qnr_tendXR,'qni_tend':qni_tendXR,'qv_tend':qv_tendXR,'temp_tend':temp_tendXR})
encoding={'qc':{'zlib':True,'complevel':9},'qr':{'zlib':True,'complevel':9},'qi':{'zlib':True,'complevel':9},'qs':{'zlib':True,'complevel':9},'qg':{'zlib':True,'complevel':9},'qns':{'zlib':True,'complevel':9},'qng':{'zlib':True,'complevel':9},'qnr':{'zlib':True,'complevel':9},'qni':{'zlib':True,'complevel':9},'qv':{'zlib':True,'complevel':9},'press':{'zlib':True,'complevel':9},'temp':{'zlib':True,'complevel':9},'dz':{'zlib':True,'complevel':9},'qc_tend':{'zlib':True,'complevel':9},'qr_tend':{'zlib':True,'complevel':9},'qi_tend':{'zlib':True,'complevel':9},'qs_tend':{'zlib':True,'complevel':9},'qg_tend':{'zlib':True,'complevel':9},'qns_tend':{'zlib':True,'complevel':9},'qng_tend':{'zlib':True,'complevel':9},'qnr_tend':{'zlib':True,'complevel':9},'qni_tend':{'zlib':True,'complevel':9},'qv_tend':{'zlib':True,'complevel':9},'temp_tend':{'zlib':True,'complevel':9}}
ds.to_netcdf('trainingData.nc',encoding=encoding)

     
