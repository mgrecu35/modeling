subroutine test(x,y)
  real, intent(out)::y
  real :: x
  y=2*x
end subroutine test


subroutine mphys_morrison_interface_2d(t2d,p2d,dz2d_wrf,qv2d,qr2d,qi2d,&
     ni2d,qs2d,qg2d,ns2d,nr2d,ng2d,w2d,nx1,nz1,dt1,&
     qc_tend2d, qi_tend2d, qni_tend2d, & 
     qr_tend2d, ni_tend2d, ns_tend2d, & 
     nr_tend2d, t_tend2d, qv_tend2d, &
     qg_tend2d, ng_tend2d)
  use mphys_morr_two_moment
  implicit none
  integer :: nx1,nz1
  integer :: i, j, k, kts, kte
  real :: dt1
  real :: t2d(nz1,nx1), p2d(nz1,nx1), dz2d_wrf(nz1,nx1),qv2d(nz1,nx1),qc2d(nz1,nx1) &
       ,qr2d(nz1,nx1), qi2d(nz1,nx1), ni2d(nz1,nx1), qs2d(nz1,nx1)         &
       , qg2d(nz1,nx1), ns2d(nz1,nx1), nr2d(nz1,nx1), ng2d(nz1,nx1)        &
       , w2d(nz1,nx1)
    
  real,intent(out) :: qc_tend2d(nz1,nx1), qi_tend2d(nz1,nx1), qni_tend2d(nz1,nx1), & 
       qr_tend2d(nz1,nx1), ni_tend2d(nz1,nx1), ns_tend2d(nz1,nx1), & 
       nr_tend2d(nz1,nx1), t_tend2d(nz1,nx1), qv_tend2d(nz1,nx1), &
       qg_tend2d(nz1,nx1), ng_tend2d(nz1,nx1)
  
  real :: precprt2d, snowrt2d
  
  real :: effc2d(nz1,nx1), effi2d(nz1,nx1), effs2d(nz1,nx1),        &
       effr2d(nz1,nx1), effg2d(nz1,nx1)
  
  real :: qrcu2d(nz,nx1), qscu2d(nz,nx1), qicu2d(nz,nx1)

  real :: qgsten2d(nz1,nx1), qrsten2d(nz1,nx1), qisten2d(nz1,nx1),  &
         qnisten2d(nz1,nx1), qcsten2d(nz1,nx1)
!----1D variables---------------------------------------!
  real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
       ,qr1d(nz), qi1d(nz), ni1d(nz), qs1d(nz)         &
       , qg1d(nz), ns1d(nz), nr1d(nz), ng1d(nz)        &
       , w1d(nz)
  
  real :: wvar1d(nz)

  real :: qc_tend1d(nz), qi_tend1d(nz), qni_tend1d(nz), & 
       qr_tend1d(nz), ni_tend1d(nz), ns_tend1d(nz), & 
       nr_tend1d(nz), t_tend1d(nz), qv_tend1d(nz), &
       qg_tend1d(nz), ng_tend1d(nz)
  
  real :: precprt1d, snowrt1d
  
  real :: effc1d(nz), effi1d(nz), effs1d(nz),        &
       effr1d(nz), effg1d(nz)
  
  real :: qrcu1d(nz), qscu1d(nz), qicu1d(nz)
  
  real :: qgsten(nz), qrsten(nz), qisten(nz),  &
       qnisten(nz), qcsten(nz)
  
print*, micro_unset
       if (micro_unset)then
          call morr_two_moment_init
          micro_unset=.False.
       end if
print*, nz, nz1, nx, nx1
  
  do i=1,nx1
     !print*,'temp',t2d(:,i)
     !print*,'press',p2d(:,i)
     !print*, 'qv',qv2d(:,i)
     do k=1,nz1
        ! zero some of these for safety
        qc_tend1d(k)  = 0.
        qi_tend1d(k)  = 0.
        qni_tend1d(k) = 0.
        qr_tend1d(k)  = 0.
        ni_tend1d(k)  = 0.
        ns_tend1d(k)  = 0.
        nr_tend1d(k)  = 0.
        t_tend1d(k)   = 0.
        qv_tend1d(k)  = 0.
        qg_tend1d(k)  = 0.
        ng_tend1d(k)  = 0.
        effc1d(k)     = 1.
        effi1d(k)     = 0.
        effs1d(k)     = 0.
        effr1d(k)     = 0.
        effg1d(k)     = 0.
        qrcu1d(k)     = 0.
        qscu1d(k)     = 0.
        qicu1d(k)     = 0.
        qgsten(k)     = 0.
        qrsten(k)     = 0.
        qisten(k)     = 0.
        qnisten(k)    = 0.
        qcsten(k)     = 0.
        
        t1d(k) = t2d(k,i)
        p1d(k) = p2d(k,i)
        dz1d(k) = dz2d_wrf(k,i)
        qv1d(k) = qv2d(k,i)
        qc1d(k) = qc2d(k,i)
        qr1d(k) = qr2d(k,i)
        nr1d(k) = nr2d(k,i)
        qi1d(k) = qi2d(k,i)
        ni1d(k) = ni2d(k,i)
        qs1d(k) = qs2d(k,i)
        ns1d(k) = ns2d(k,i)
        qg1d(k) = qg2d(k,i)
        ng1d(k) = ng2d(k,i)
        
        wvar1d(k) = 0.5 ! hard-wired not coupled to forcing!
        w1d(k) = w2d(k,i)*0
     end do
     
    
     j=1
     kts=1
     kte=nz1
     !print*, qc1d
     !print*, qr1d
     !print*, qc1d
     !print*, t1d
     !print*, 'dt1',i,j
     call morr_two_moment_micro(qc_tend1d(1:nz1), qi_tend1d(1:nz1), qni_tend1d(1:nz1), &
          qr_tend1d(1:nz1), ni_tend1d(1:nz1), ns_tend1d(1:nz1), nr_tend1d(1:nz1),             &
          qc1d(1:nz1), qi1d(1:nz1), qs1d(1:nz1), qr1d(1:nz1), ni1d(1:nz1), &
          ns1d(1:nz1), nr1d(1:nz1),               &
          t_tend1d(1:nz1), qv_tend1d(1:nz1), t1d(1:nz1), qv1d(1:nz1), &
          p1d(1:nz1), dz1d(1:nz1), w1d(1:nz1), wvar1d(1:nz1), &
          precprt1d, snowrt1d,                                    &
          effc1d(1:nz1), effi1d(1:nz1), effs1d(1:nz1), effr1d(1:nz1), dt1,                    &
          1,1,j,j,kts,nz1,                                        &
          1,1,j,j,kts,nz1,                                        &
          qg_tend1d(1:nz1), ng_tend1d(1:nz1), qg1d(1:nz1), &
          ng1d(1:nz1), effg1d(1:nz1),               & !effg1d graupel effective radius (microns)
          qrcu1d(1:nz1), qscu1d(1:nz1), qicu1d(1:nz1),                                 & !cumulus tendencies
          qgsten(1:nz1), qrsten(1:nz1), qisten(1:nz1), qnisten(1:nz1), qcsten(1:nz1))                  !sedimentation tendencies
     ! t_tend1d, dt1
     qc_tend2d(:,i)=qc_tend1d(1:nz1)
     qi_tend2d(:,i)=qi_tend1d(1:nz1)
     qni_tend2d(:,i)=qni_tend1d(1:nz1)
     qr_tend2d(:,i)=qr_tend1d(1:nz1)
     ni_tend2d(:,i)=ni_tend1d(1:nz1)
     ns_tend2d(:,i)=ns_tend1d(1:nz1)
     nr_tend2d(:,i)=nr_tend1d(1:nz1)
     t_tend2d(:,i)=t_tend1d(1:nz1)
     qg_tend2d(:,i)=qg_tend1d(1:nz1)
     ng_tend2d(:,i)=ng_tend1d(1:nz1)
     qv_tend2d(:,i)=qv_tend1d(1:nz1)
          ! save tendencies
       !do k=1,nz
       !   dtheta_mphys(k,i)=t_tend1d(k)/exner(k,i)
       !   dqv_mphys(k,i)=qv_tend1d(k)
       !   dhydrometeors_mphys(k,i,1)%moments(1,1)= qc_tend1d(k)
       !   dhydrometeors_mphys(k,i,2)%moments(1,1)= qr_tend1d(k)
       !   dhydrometeors_mphys(k,i,2)%moments(1,2)= nr_tend1d(k)
       !   dhydrometeors_mphys(k,i,3)%moments(1,1)= qi_tend1d(k)
       !   dhydrometeors_mphys(k,i,3)%moments(1,2)= ni_tend1d(k)
       !   dhydrometeors_mphys(k,i,4)%moments(1,1)= qni_tend1d(k)
       !   dhydrometeors_mphys(k,i,4)%moments(1,2)= ns_tend1d(k)
       !   dhydrometeors_mphys(k,i,5)%moments(1,1)= qg_tend1d(k)
       !   dhydrometeors_mphys(k,i,5)%moments(1,2)= ng_tend1d(k)
       !end do
       
  end do
  
end subroutine mphys_morrison_interface_2d
