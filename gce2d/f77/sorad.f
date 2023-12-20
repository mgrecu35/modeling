cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***************  charney  /silo/z1mhy/peng/solar.f ***10/20/95******
      subroutine sorad (m,n,ndim,pl,ta,wa,oa,co2,taucld,reff,taual,
     1             rsirbm,rsuvbm,cosz,flx,fdirir,fdifir,fdirpar,fdifpar)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1
      PARAMETER (NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NXX-2,NBB=1,MM=NXI/NBB/2,N1=1)
      real rflux(mm,7)
      common /radflux/rflux
c********************************************************************
c  
c This routine computes solar fluxes due to the absoption by water
c  vapor, ozone, co2, o2, clouds, and aerosols and due to the
c  scattering by clouds, aerosols, and gases.
c
c This is a vectorized code. It computes the fluxes simultaneous for
c  (m x n) soundings, which is a subset of the (m x ndim) soundings.
c  In a global climate model, m and ndim correspond to the numbers of
c  grid boxes in the zonal and meridional directions, respectively.
c
c Ice and liquid cloud particles are allowed to co-exist in any of the
c  np layers. Two sets of cloud parameters are required as inputs, one
c  for ice paticles and the other for liquid particles.  These parameters
c  are optical thickness (taucld) and effective particle size (reff).
c
c If no information is available for reff, a default value of
c  10 micron for liquid water and 75 micron for ice can be used.
c
c Clouds are grouped into high, middle, and low clouds separated by the
c  level indices ict and icb.  For detail, see the subroutine cldscale.
c
c----- Input parameters:                           
c                                                   units      size
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c           meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   level pressure (pl)                              mb       m*ndim*(np+1)
c   layer temperature (ta)                           k        m*ndim*np
c   layer specific humidity (wa)                     gm/gm    m*ndim*np
c   layer ozone concentration (oa)                   gm/gm    m*ndim*np
c   co2 mixing ratio by volumn (co2)               parts/part   1
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c                index 1 for ice particles
c                index 2 for liquid drops
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c                index 1 for ice particles
c                index 2 for liquid drops
c   aerosol optical thickness (taual)                n/d      m*ndim*np 
c   solar ir surface albedo for beam                fraction   m*ndim
c                radiation (rsirbm)                
c   uv + par surface albedo for beam                     fraction   m*ndim
c                radiation (rsuvbm)                
c   cosine of solar zenith angle (cosz)            n/d        m*ndim
c
c----- Output parameters
c
c   all-sky flux (downward minus upward) (flx)     fraction   m*ndim*(np+1)
c   all-sky direct downward ir (0.7-10 micron)
c                flux at the surface (fdirir)      fraction   m*ndim
c   all-sky diffuse downward ir flux at
c                the surface (fdifir)              fraction   m*ndim
c   all-sky direct downward par (0.4-0.7 micron)
c                flux at the surface (fdirpar)     fraction   m*ndim
c   all-sky diffuse downward par flux at
c                the surface (fdifpar)             fraction   m*ndim
*
c----- Notes:
c
c    (1) The unit of flux is fraction of the incoming solar radiation
c        at the top of the atmosphere.  Therefore, fluxes should
c        be equal to flux multiplied by the extra-terrestrial solar
c        flux and the cosine of solar zenith angle.
c    (2) Clouds and aerosols can be included in any layers by specifying
c        taucld(i,j,k,*) and taual(i,j,k), k=1,np. 
c        For an atmosphere without clouds and aerosols,
c        set taucld(i,j,k,*)=taual(i,j,k)=0.0.
c    (3) Aerosol single scattering albedos and asymmetry
c        factors are specified in the subroutines solir and soluv.
c    (4) pl(i,j,1) is the pressure at the to of the model, and
c        pl(i,j,np+1) is the surface pressure.
c    (5) the pressure levels ict and icb correspond approximately
c        to 400 and 700 mb.
c        
c**************************************************************************

c-----input parameters

      integer m,n,ndim
      real    pl(m,ndim,np+1),ta(m,ndim,np),wa(m,ndim,np),oa(m,ndim,np)
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2)
      real    taual(m,ndim,np),rsirbm(m,ndim)
      real    rsuvbm(m,ndim),cosz(m,ndim),co2

c-----output parameters

      real    flx(m,ndim,np+1)
      real    fdirir(m,ndim),fdifir(m,ndim)
      real    fdirpar(m,ndim),fdifpar(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array
 
      integer i,j,k,mmm
      real    dp(MM,N1,np),wh(MM,N1,np),oh(MM,N1,np),scal(MM,N1,np)
      real    swh(MM,N1,np+1),so2(MM,N1,np+1),df(MM,N1,np+1)
      real    sdf(MM,N1),sclr(MM,N1),csm(MM,N1),taux,x,expmn

      save
 
c-----------------------------------------------------------------
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n 
       do i= 1, mmm 

         swh(i,j,1)=0. 
         so2(i,j,1)=0. 

c-----csm is the effective secant of the solar zenith angle
c     see equation (12) of Lacis and Hansen (1974, JAS)    
 
         csm(i,j)=35./sqrt(1224.*cosz(i,j)*cosz(i,j)+1.)

       enddo 
      enddo

      do k= 1, np
       do j= 1, n
         do i= 1, mmm

c-----compute layer thickness and pressure-scaling function. 
c     indices for the surface level and surface layer
c     are np+1 and np, respectively.
 
          dp(i,j,k)=pl(i,j,k+1)-pl(i,j,k)
          scal(i,j,k)=dp(i,j,k)*(.5*(pl(i,j,k)+pl(i,j,k+1))/300.)**.8
 
c-----compute scaled water vapor amount, unit is g/cm**2

          wh(i,j,k)=1.02*wa(i,j,k)*scal(i,j,k)*
     *              (1. + 0.00135*(ta(i,j,k)-240.))
          swh(i,j,k+1)=swh(i,j,k)+wh(i,j,k)

c-----compute ozone amount, unit is (cm-atm)stp.
 
          oh(i,j,k)=1.02*oa(i,j,k)*dp(i,j,k)*466.7

        enddo
       enddo
      enddo

c-----initialize fluxes for all-sky (flx) and flux reduction (df)

      do k=1, np+1
       do j=1, n
        do i=1, mmm
          flx(i,j,k)=0.
          df(i,j,k)=0.
        enddo
       enddo
      enddo

c-----compute solar ir fluxes

      call solir (m,n,ndim,wh,taucld,reff,taual,
     *            csm,rsirbm,flx,fdirir,fdifir)

c-----compute and update uv and par fluxes

      call soluv (m,n,ndim,oh,dp,taucld,reff,taual,
     *            csm,rsuvbm,flx,fdirpar,fdifpar)

c-----compute scaled amount of o2 (so2), unit is (cm-atm)stp.

      do k= 1, np
       do j= 1, n
        do i= 1, mmm
          so2(i,j,k+1)=so2(i,j,k)+165.22*scal(i,j,k)
        enddo
       enddo
      enddo

c-----compute flux reduction due to oxygen following
c      chou (J. climate, 1990). The fraction 0.0287 is the
c      extraterrestrial solar flux in the o2 bands.

       do k= 2, np+1
        do j= 1, n
         do i= 1, mmm
           x=so2(i,j,k)*csm(i,j)
           df(i,j,k)=df(i,j,k)+0.0287*(1.-expmn(-0.00027*sqrt(x)))
         enddo
        enddo
       enddo          

c-----compute scaled amounts for co2 (so2). unit is (cm-atm)stp.

      do k= 1, np
       do j= 1, n
        do i= 1, mmm
         so2(i,j,k+1)=so2(i,j,k)+co2*789.*scal(i,j,k)
        enddo
       enddo
      enddo

c-----compute and update flux reduction due to co2 following
c     chou (J. Climate, 1990)

      call flxco2(m,n,np,so2,swh,csm,df)

c-----adjust for the all-sky fluxes due to o2 and co2.  It is
c     assumed that o2 and co2 have no effects on solar radiation
c     below clouds.

      do j=1,n
       do i=1,mmm
          sdf(i,j)=0.0 
          sclr(i,j)=1.0 
       enddo
      enddo

      do k=1,np
       do j=1,n
        do i=1,mmm

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
         if(taux.gt.0.01 .and. sclr(i,j).eq.1.) then
          sdf(i,j)=df(i,j,k)
          sclr(i,j)=0.0
         endif

          flx(i,j,k+1)=flx(i,j,k+1)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
 
C --------------- GCSS -------------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
         endif
C ---------------------------------------------------------------
         


        enddo
       enddo
      enddo

c-----adjust for the direct downward ir flux.
      do j= 1, n
       do i= 1, mmm
           fdirir(i,j)=fdirir(i,j)-sdf(i,j)-df(i,j,np+1)*sclr(i,j)
       enddo
      enddo

      return
      end  
