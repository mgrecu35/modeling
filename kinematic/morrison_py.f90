subroutine test(x,y)
  real, intent(out)::y
  real :: x
  y=2*x
end subroutine test


subroutine mphys_morrison_interface_2d(t2d,p2d,dz2d,qv2d,qr2d,qi2d,&
     ni2d,qs2d,qg2d,ns2d,nr2d,ng2d,w2d,nx1,nz1,dt1)
  use mphys_morr_two_moment
  implicit none
  integer :: nx1,nz1
  integer :: i, j, k, kts, kte
  real :: dt1
  real :: t2d(nz1,nx1), p2d(nz1,nx1), dz2d(nz1,nx1),qv2d(nz1,nx1),qc2d(nz1,nx1) &
       ,qr2d(nz1,nx1), qi2d(nz1,nx1), ni2d(nz1,nx1), qs2d(nz1,nx1)         &
       , qg2d(nz1,nx1), ns2d(nz1,nx1), nr2d(nz1,nx1), ng2d(nz1,nx1)        &
       , w2d(nz1,nx1)
    
  real :: qc_tend2d(nz1,nx1), qi_tend2d(nz1,nx1), qni_tend2d(nz1,nx1), & 
       qr_tend2d(nz1,nx1), ni_tend2d(nz1,nx1), ns_tend2d(nz1,nx1), & 
       nr_tend2d(nz1,nx1), t_tend2d(nz1,nx1), qv_tend2d(nz1,nx1), &
       qg_tend2d(nz1,nx1), ng_tend2d(nz1,nx1)
  
  real :: precprt2d, snowrt2d
  
  real :: effc2d(nz1,nx1), effi2d(nz1,nx1), effs2d(nz1,nx1),        &
       effr2d(nz1,nx1), effg2d(nz1,nx1)
  
  real :: qrcu2d(nz1,nx1), qscu2d(nz1,nx1), qicu2d(nz1,nx1)
  
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
  
  do i=1,nx
     do k=1,nz
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
        effc1d(k)     = 0.
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
        !dz1d(k) = dz(k)
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
        w1d(k) = w2d(k,i)
     end do
     
     ! Initialise microphysics 
     if (micro_unset)then
        call morr_two_moment_init
        micro_unset=.False.
     end if
     
     
     call morr_two_moment_micro(qc_tend1d, qi_tend1d, qni_tend1d, &
          qr_tend1d, ni_tend1d, ns_tend1d, nr_tend1d,             &
          qc1d, qi1d, qs1d, qr1d, ni1d, ns1d, nr1d,               &
          t_tend1d, qv_tend1d, t1d, qv1d, p1d, dz1d, w1d, wvar1d, &
          precprt1d, snowrt1d,                                    &
          effc1d, effi1d, effs1d, effr1d, dt1,                    &
          i,i,j,j,kts,kte,                                        &
          i,i,j,j,kts,kte,                                        &
          qg_tend1d, ng_tend1d, qg1d, ng1d, effg1d,               & !effg1d graupel effective radius (microns)
          qrcu1d, qscu1d, qicu1d,                                 & !cumulus tendencies
          qgsten, qrsten, qisten, qnisten, qcsten)                  !sedimentation tendencies

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
