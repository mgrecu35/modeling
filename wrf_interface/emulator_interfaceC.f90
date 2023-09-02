module scalers
  real :: qcm, qcs
  real :: qrm, qrs
  real :: qim, qis
  real :: qsm, qss
  real :: qgm, qgs
  real :: qvm, qvs
  real :: pressm, presss
  real :: tempm, temps
  real :: dzm, dzs
  real :: qc_tendm, qc_tends
  real :: qr_tendm, qr_tends
  real :: qi_tendm, qi_tends
  real :: qs_tendm, qs_tends
  real :: qg_tendm, qg_tends
  real :: qv_tendm, qv_tends
  real :: temp_tendm, temp_tends
  integer :: itime_ij(30)
  
end module scalers

!nx_py=ite-its+1
!ny_py=jte-jts+1
!nz_py=kte-kts+1
subroutine emulator_interface_f90(th,&
    qv_curr,qc_curr,qr_curr,&
    qi_curr,qs_curr,qg_curr,&
    qni_curr,qns_curr,qnr_curr,&
    qng_curr,rho,pi_phy,&
    p,dt,dz8w,&
    nx,ny,nz,ireturn,ij)
  use scalers
  use keras_def
  use MODULE_MP_MORR_TWO_MOMENT2, only : mp_MORR_TWO_MOMENT2, SIMPLE_MICRO3
implicit none
real,intent(inout)::th(nx,nz,ny)
real,intent(inout)::qv_curr(nx,nz,ny)
real,intent(inout)::qc_curr(nx,nz,ny)
real,intent(inout)::qr_curr(nx,nz,ny)
real,intent(inout)::qi_curr(nx,nz,ny)
real,intent(inout)::qs_curr(nx,nz,ny)
real,intent(inout)::qg_curr(nx,nz,ny)
real,intent(inout)::qni_curr(nx,nz,ny)
real,intent(inout)::qns_curr(nx,nz,ny)
real,intent(inout)::qnr_curr(nx,nz,ny)
real,intent(inout)::qng_curr(nx,nz,ny)
real,intent(inout)::rho(nx,nz,ny)
real,intent(inout)::pi_phy(nx,nz,ny)
real,intent(inout)::p(nx,nz,ny)
real,intent(inout)::dt
real,intent(inout)::dz8w(nx,nz,ny)
integer, intent(in)::nx,ny,nz,ij
integer, intent(inout) :: ireturn
integer :: i,k,j
real :: temp_curr

!real :: dt1
integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

integer itimestep
REAL::    RAINNC(nx,ny), RAINNCV(nx,ny), SR(nx,ny), refl_10cm(nx,nz,ny), ht(nx,ny),&
     qrcuten(nx,nz,ny), qscuten(nx,nz,ny), qicuten(nx,nz,ny), mu(nx,ny), w(nx,nz,ny),&
     temp(nx,nz,ny)
real, dimension(nx,nz,ny):: TH_old, QV_old, QC_old, QR_old, QI_old, QS_old, &
     QG_old, qNI_old, qNS_old, qNR_old, qNG_old, p_old

real, dimension(nx,nz,ny):: TH_em, QV_em, QC_em, QR_em, QI_em, QS_em, &
     QG_em, qNI_em, qNS_em, qNR_em, qNG_em

real, dimension(nz) :: QC3DTEN,QI3DTEN,QS3DTEN,QR3DTEN,         &
                    NI3DTEN,NS3DTEN,NR3DTEN, &
                    QG3DTEN,NG3DTEN,NG3D,   &
                    QGSTEN,QRSTEN,QISTEN,QSSTEN,QCSTEN,W3D,WVAR, &
                    T3DTEN,QV3DTEN,t3d
real :: snowrt, precrt

character(len=7) :: s1
character (len=13) :: s2   
integer :: iemulator
ids=1
ims=1
its=1
jds=1
jms=1
jts=1
kds=1
kms=1
kts=1
ide=nx
ime=nx
ite=nx
jde=ny
jme=ny
jte=ny
kde=nz
kme=nz
kte=nz

qrcuten=0
qscuten=0
w=0
mu=1
ht=0
if(itime_ij(ij)==0) then
     qr_curr=0
     qs_curr=0
     qg_curr=0
     qni_curr=0
     qns_curr=0
     qnr_curr=0
     qng_curr=0
     qc_curr=0
     qi_curr=0
endif
itime_ij(ij)=itime_ij(ij)+1
ireturn=1
!return
th_old=th
qv_old=qv_curr
qc_old=qc_curr
qr_old=qr_curr
qs_old=qs_curr
qg_old=qg_curr
qi_old=qi_curr
!return
!call updatefields(qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,th,p,dz8w,&
!     qv_em,qc_em,qr_em,qi_em,qs_em,qg_em,th_em,nx,ny,nz,ij)

do i=1,nx
     do j=1,ny
          t3d=th(i,:,j)*pi_phy(i,:,j)
          call SIMPLE_MICRO3(QC3DTEN,QI3DTEN,QS3DTEN,QR3DTEN,         &
            NI3DTEN,NS3DTEN,NR3DTEN,qc_curr(i,:,j),qi_curr(i,:,j), qs_curr(i,:,j),&
            qr_curr(i,:,j),qni_curr(i,:,j),qns_curr(i,:,j),qnr_curr(i,:,j), &
            T3DTEN,QV3DTEN,t3d,qv_curr(i,:,j),p(i,:,j),dz8w(i,:,j),&
            W3D,WVAR,PRECRT,SNOWRT,    &
            DT,                                                   &
            1,nz,           & ! ADD GRAUPEL
            QG3DTEN,NG3DTEN,qg_curr(i,:,j),NG3D,   &
            QGSTEN,QRSTEN,QISTEN,QSSTEN,QCSTEN)
          th(i,:,j)=t3d/pi_phy(i,:,j)
     enddo
enddo
return
if(1==-1) then
   call mp_MORR_TWO_MOMENT2(ITIMESTEP,                       &
     TH, QV_curr, QC_curr, QR_curr, QI_curr, QS_curr, &
     QG_curr, qNI_curr, qNS_curr, qNR_curr, qNG_curr, &
     RHO, Pi_phy, P, DT, DZ8w, HT, W,          &
     RAINNC, RAINNCV, SR,                    &
     qrcuten, qscuten, qicuten, mu           & ! hm added
     ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
     ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
     ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
     )
   !return
endif
iemulator=1
if(iemulator==-1) then
     th=th_em
     qv_curr=qv_em
     qc_curr=qc_em
     qr_curr=qr_em
     qi_curr=qi_em
     qs_curr=qs_em
     qg_curr=qg_em
     qni_curr=0
     qns_curr=0
     qnr_curr=0
     qng_curr=0
endif
ireturn=1
!return

if(mod(itime_ij(ij),14)==-1) then
   print*, 'time for another fresca', itime_ij(1:10), ij
   write(s1, '(A, I2.2)') 'tile_', ij
   write(s2, '(A, I3.3)') '/output_em', itime_ij(ij)/14
   open(11+ij,file=s1//s2)
   write(11+ij,*) nx, ny, nz
   close(11+ij)
   open(11+ij,file=s1//s2//'.bin',form='unformatted')
   write(11+ij) qv_old(:,1:72,:)
   write(11+ij) qc_old(:,1:72,:)
   write(11+ij) qr_old(:,1:72,:)
   write(11+ij) qi_old(:,1:72,:)
   write(11+ij) qs_old(:,1:72,:)
   write(11+ij) qg_old(:,1:72,:)
   write(11+ij) th_old(:,1:72,:)
   write(11+ij) qv_curr(:,1:72,:)
   write(11+ij) qc_curr(:,1:72,:)
   write(11+ij) qr_curr(:,1:72,:)
   write(11+ij) qi_curr(:,1:72,:)
   write(11+ij) qs_curr(:,1:72,:)
   write(11+ij) qg_curr(:,1:72,:)
   write(11+ij) th(:,1:72,:)
   write(11+ij) p(:,1:72,:)
   write(11+ij) dz8w(:,1:72,:)
   write(11+ij) qv_em(:,1:72,:)
   write(11+ij) qc_em(:,1:72,:)
   write(11+ij) qr_em(:,1:72,:)
   write(11+ij) qi_em(:,1:72,:)
   write(11+ij) qs_em(:,1:72,:)
   write(11+ij) qg_em(:,1:72,:)
   write(11+ij) th_em(:,1:72,:)
   close(11+ij)
end if

ireturn=1

end subroutine emulator_interface_f90
