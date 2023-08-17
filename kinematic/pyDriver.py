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

xlong,xlat,press,theta,qv,qs,qg,qr,qc,qi,u,v,height_surf,\
    dbz,height,qns,qng,qnr,qnice=readWRF_File_Morr(fs[7])
dbzm=np.ma.array(dbz,mask=(dbz<0.0))
n1,n2=50,175
for j0 in range(120,180):
     qs2d=qs[-1,:,j0,n1:n2]
     qr2d=qr[-1,:,j0,n1:n2]
     qc2d=qc[-1,:,j0,n1:n2]
     qi2d=qi[-1,:,j0,n1:n2]
     qg2d=qg[-1,:,j0,n1:n2]
     qns2d=qns[-1,:,j0,n1:n2]
     qnr2d=qnr[-1,:,j0,n1:n2]
     qni2d=qnice[-1,:,j0,n1:n2]
     qng2d=qng[-1,:,j0,n1:n2]
     qv2d=qv[-1,:,j0,n1:n2]
     press2d=press[-1,:,j0,n1:n2]
     temp2d=theta[-1,:,j0,n1:n2]*(press2d/1e5)**(287/1004)
     dtemp2d_wrf=(theta[-1,:,j0,n1:n2]-theta[-2,:,j0,n1:n2])*(press2d/1e5)**(287/1004)
     dz2d_wrf=height[-1,1:,j0,n1:n2]-height[-1,:-1,j0,n1:n2]
     dbz2d=dbzm[-1,:,j0,n1:n2]
     w2d=qv2d.copy()*0.0
     dt1=5.0
     
     #qc_tend2d,qi_tend2d,qni_tend2d,qr_tend2d,ni_tend2d,ns_tend2d,nr_tend2d,t_tend2d,qv_tend2d,qg_tend2d,ng_tend2d = mphys_morrison_interface_2d(t2d,p2d,dz2d_wrf,qv2d,qr2d,qi2d,ni2d,qs2d,qg2d,ns2d,nr2d,ng2d,w2d,dt1,[nx1,nz1])
     if 1==2:
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
     dt1=10
     qc_tend2d,qi_tend2d,qs_tend2d,qr_tend2d,\
          ni_tend2d,ns_tend2d,nr_tend2d,t_tend2d,qv_tend2d,\
          qg_tend2d,ng_tend2d=kd.mphys_morrison_interface_2d(temp2d[:70,:].data,press2d[:70,:].data,dz2d_wrf[:70,:].data,\
                                                             qv2d[:70,:].data,qr2d[:70,:].data,\
                                                             qi2d[:70,:].data,qni2d[:70,:].data,qs2d[:70,:].data,qg2d[:70,:].data,qns2d[:70,:].data,\
                                                             qnr2d[:70,:].data,qng2d[:70,:].data,w2d[:70,:].data,dt1)
     
     h1=dz2d_wrf[:,0].cumsum()
     plt.figure(figsize=(8,10))
     plt.subplot(211)
     dbz2dm=np.ma.array(dbz2d[:70,:],mask=(dbz2d[:70,:]<0.0))
     plt.pcolormesh(np.arange(w2d.shape[1]),h1[:70]/1e3,dbz2dm[:70,:],cmap='jet',vmin=0,vmax=60)
     plt.title('Reflectivity')
     plt.ylim(0,15.0)
     plt.colorbar()
     plt.subplot(212)
     plt.pcolormesh(np.arange(w2d.shape[1]),h1[:70]/1e3,t_tend2d,cmap='RdBu_r',vmin=-0.005,vmax=0.005)
     plt.title('Temperature tendency')
     plt.ylim(0,15.0)
     plt.colorbar()
     
     plt.savefig('tempTend_%3.3i.png'%j0)
     plt.close('all')
     
