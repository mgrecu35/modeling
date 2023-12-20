cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine radrat (iflag,cosz,npp1,pl,pa)
c nb is a number of blocks for calculation
      integer nbb,nq,nx,nz,lay,nadd,nxi,nw,nx1,nx2,nz2,nz4,nz15,nx3,mb1
      integer mb2,iflag,npp1
      parameter(NX=514,NZ=43,lay=88,nadd=7,nbb=1)
      parameter(nq=nz+nadd-1,nxi=nx-2,nw=nq-1,nx1=nxI/nbb/2,nx2=nxi/nbb)
      parameter(nz2=nz*2,nz4=nz*4,nz15=nz*15,nx3=nx*3)
      parameter(mb1=nz15,mb2=nz*7+nx3)

      real    cosz,pl(lay),pa(lay)
      integer iradave
      common/radflux/rflux(nxi/2,7)
      common/iptionr/ iradave
      common/iceopt/ ice913,ilif
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      real    rsw(nx,nz),rlw(nx,nz)
      common/slwave/ rsw,rlw
      real    tb(nz),qb(nz),zz0(nz),rho(nz),zz1(nz2),tb1(nz),
     $        qb1(nz),zz2(mb1)
      common/b5/ tb,qb,zz0,rho,zz1,tb1,qb1,zz2
      real    zz3(nz4),p0(nz),pi(nz),zz4(mb2)
      common/b6/ zz3,p0,pi,zz4
      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),pairsfc(nx),
     1        thairsf(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc
      real    rsirbm(nx1,1),rsirdf(nx1,1),rsuvbm(nx1,1),
     $        rsuvdf(nx1,1),taual(nx1,1,nw)
      common/albedo1/ rsirbm,rsirdf,rsuvbm,rsuvdf,taual
      integer ict,icb
      common/cloudp/ ict,icb
      real    sun_4(nx,4)
      common/surface/ sun_4
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    ta(lay),wa(lay),oa(lay),plu(nadd)
      real    fcld(nx1,1,nw),tauir1(nx1,1,nw),flx(nx1,1,nq)
      real    flc(nx1,1,nq),dfdts(nx1,1,nq),st4(nx1,1),tsq(nx1,1)
      real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)
      real    coolr1(nx,1,nw),heatr1(nx,1,nw)
      real    tausw1(nx1,1,nw,2),reff(nx1,1,nw,2)
      real    flx1(nx1,1,nq)
      real    cosz1(nx1,1)
      real    fdirir(nx1,1),fdifir(nx1,1)
      real    fdirpar(nx1,1),fdifpar(nx1,1)
      logical high
      data high/.false./
      data plu/ 0.01, 0.56, 3.04, 4.04, 7.18,  12.38, 22.46/
      integer i,k,irst,np,km,ib,i1,ii,i3,k1,ihalfp
      real    aco,cmc,sc0
      DATA IRST/0/
      save

      do k=1,nq
      do i=1,nx1
        flx(i,1,k)=0.0
        flx1(i,1,k)=0.0
        flc(i,1,k)=0.0
        dfdts(i,1,k)=0.0
      enddo
      enddo
      do k=1,nw
      do i=1,nx
        coolr1(i,1,k)=0.0
        heatr1(i,1,k)=0.0
      enddo
      enddo
      do i=1,nx1
        st4(i,1)=0.0
        fdirir(i,1)=0.0
        fdifir(i,1)=0.0
        fdirpar(i,1)=0.0
        fdifpar(i,1)=0.0
      enddo

      do i=1,lay
        pa(i)=0.0
      enddo
c
      IF (IRST.NE.0) GO TO 500
c
c     np=number of atmospheric layers; npp1=surface level
      np  =kl2+nadd
      npp1=np+1
c      cpi =4.*atan(1.)
c  solar constant and cosine of solar zenith angle
c      sc0=1365.
      sc0=1411.
      if (ilif .eq. 1) sc0=1365.
c     cosz=cos(51.74*cpi/180.)
        do i=1,nx1
          rsirbm(i,1)=0.07
          rsirdf(i,1)=0.07
          rsuvbm(i,1)=0.07
          rsuvdf(i,1)=0.07
        enddo
c
        ict=24
        icb=32
c 
c  assign co2 (cmc). units are parts/part
      cmc=300.e-6
      write(6,*) 'npp1(3)=',npp1
c     write(6,*) 'npp1(3)=',npp1,'  nadd = ',nadd
c  assign aerosol optical thickness 
      do 5 k=1,nw
      do 5 i=1,nx1
        taual(i,1,k)=0.0
    5 continue

      do k=1,kl2+1
        km=kl2+3-k
        pl(k+nadd)=1.e-3*p00(km)
      end do

      do k=1,nadd
        pl(k)=plu(k)
      end do

      do k=1,kl2
        km=kl2+2-k
        pa(k+nadd)=1.e-3*p0(km)
      end do

      do 13 k=1,nadd
        pa(k)=0.5*(pl(k)+pl(k+1))
        wa(k)=1.e-6
c        ta(k)=190.
   13 continue
c     print *,'       radrat.f          5  fito3.f'
      call fito3 (npp1,pa,ta,oa)
       print*
      print*,' k       pa        pl       taual        ao        ta
     *   wa'
      write (6,975) (k,pa(k),pl(k),taual(1,1,k),oa(k)
     1 ,ta(k),wa(k),k=1,np)    ! shie 2/22/01
c    1 ,ta(k),wa(k),k=1,npp1)

      k=npp1                   ! shie 2/22/01
      write (6,975) k,pa(k),pl(k),taual(1,1,np),oa(k) ! shie 2/22/01
     1 ,ta(k),wa(k)    ! shie 2/22/01
c
      IRST=1
c
  500 continue
  975 format(i3,3f10.3,e12.4,f12.5,f12.7)
        do i=1,nx1
          cosz1(i,1)=cosz
        enddo

      DO 1000 IB=1,NBB
        I1=(IB-1)*NX2+1
        CALL OPT4 (I1,NX1,NX2,PL,TA,WA,OA,TAUSW1,TAUIR1
     +                 ,FCLD,TAQ,WAQ,OAQ,PLQ,TSQ,REFF)

C ----------------- GCSS -----------------------------------

        do i=1,nx1
          rflux(i,1)=0.
          rflux(i,2)=0.
          rflux(i,3)=0.
          rflux(i,4)=0.
          rflux(i,5)=0.
          rflux(i,6)=0.
          rflux(i,7)=0.
        enddo

C ------------------------------------------------------------

        if (cosz.ge.0.005) then

          call sorad (nx1,1,1,plq,taq,waq,oaq,cmc,
     $                 tausw1,reff,taual,rsirbm,rsuvbm,
     1                 cosz1,flx1,fdirir,fdifir,fdirpar,fdifpar)
        end if

C ----------------- GCSS -----------------------------------

        do i=1,nx1
          rflux(i,2)=0.
        enddo

C ------------------------------------------------------------

        call irrad (nx1,1,1,tauir1,fcld,plq,taq,waq,oaq,cmc,tsq,
     $                  high,flx,flc,dfdts,st4)
c
c since flx1 is normalized, they should be multiplied
c by sc0*cosz
c since upward flux should be positive, heatr is multiplied by
c a minus sign
c save solar and long-wave radiative fluxes
C
C
        DO I=I1,I1+NX2-1,2
          II=1
          IF (IRADAVE .EQ. 0) II=(I-I1)/2+1
           i3=i+1
          sun_4(i3,1)=flx1(ii,1,nw+1)*sc0*cosz1(ii,1)
          sun_4(i3,2)=flx(ii,1,nw+1)
        enddo

C
        do k=1,nw
        DO I=I1,I1+NX2-1,2
          II=1
          IF (IRADAVE .EQ. 0) II=(I-I1)/2+1
           I3=I+1
          heatr1(i3,1,k)=(flx1(ii,1,k+1)-flx1(ii,1,k))*8.441874/
     1       (plq(ii,1,k+1)-plq(ii,1,k))
          heatr1(i3,1,k)=-heatr1(i3,1,k)*sc0*cosz1(ii,1)
          coolr1(i3,1,k)=(flx(ii,1,k+1)-flx(ii,1,k))*8.441874/
     1      (plq(ii,1,k+1)-plq(ii,1,k))
        enddo
        enddo
C
        do k1=nadd+1,np
          k=np+2-k1
        DO I=I1,I1+NX2-1,2
            i3=i+1
           aco=1./(86400.*pi(k))
          rsw(i3,k)=aco*heatr1(i3,1,k1)
          rlw(i3,k)=aco*coolr1(i3,1,k1)
        enddo
        enddo

1000   continue
CC
        IF (IRADAVE .EQ. 0) THEN
C
        DO I=2,NX-2,2
            i3=i+1
          sun_4(i3,1)=0.5*(sun_4(i3-1,1)+sun_4(i3+1,1))
          sun_4(i3,2)=0.5*(sun_4(i3-1,2)+sun_4(i3+1,2))
        enddo
        sun_4(nx-1,1)=.5*(sun_4(nx-2,1)+sun_4(2,1))
        sun_4(nx-1,2)=.5*(sun_4(nx-1,2)+sun_4(2,2))
        do k1=nadd+1,np
          k=np+2-k1
        DO I=2,NX-2,2
           i3=i+1
          rsw(i3,k)=0.5*(rsw(i3-1,k)+rsw(i3+1,k))
          rlw(i3,k)=0.5*(rlw(i3-1,k)+rlw(i3+1,k))
        enddo
          rsw(nx-1,k)=.5*(rsw(nx-2,k)+rsw(2,k))
          rlw(nx-1,k)=.5*(rlw(nx-2,k)+rlw(2,k))
        enddo
CC
       ELSE
C
        DO I=2,NX-2,2
            I3=I+1
          SUN_4(I3,1)=SUN_4(2,1)
          SUN_4(I3,2)=SUN_4(2,2)
        ENDDO
        SUN_4(NX-1,1)=SUN_4(2,1)
        SUN_4(NX-1,2)=SUN_4(2,2)
        DO K1=NADD+1,NP
          K=NP+2-K1
        DO I=2,NX-2,2
           I3=I+1
          RSW(I3,K)=RSW(2,K)
          RLW(I3,K)=RLW(2,K)
        ENDDO
          RSW(NX-1,K)=RSW(2,K)
          RLW(NX-1,K)=RLW(2,K)
        ENDDO
       ENDIF

c ------------------------- for GCSS workshop ----------------------------

       do i=1,nx1
         rflux(i,1)=rflux(i,1)*sc0*cosz1(i,1)
         rflux(i,3)=rflux(i,3)*sc0*cosz1(i,1)
c         rflux(i,3)=1411*sc0*cosz1(i,1)
         rflux(i,4)=flx1(i,1,1)*sc0*cosz1(i,1)-rflux(i,3)
         rflux(i,6)=flx1(i,1,npp1)*sc0*cosz1(i,1)-rflux(i,1)
       enddo


c -----------------------------------------------------------------------

      if(iflag.eq.2) then
        ihalfp=il2/2+1
        write(6,*) 'cosz=',cosz
        write(6,*) 'check rsw(i,k) and rlw(i,k) at central point'
        DO K1=NADD+1,NP
          K=NP+2-K1
          write(6,7821) k,pl(k),rsw(ihalfp,k),rlw(ihalfp,k),
     1        heatr1(ihalfp,1,k1),coolr1(ihalfp,1,k1)
        enddo
      endif
7821  format(2x,i6,2x,f10.3,2x,4e20.10)
      return
      end
