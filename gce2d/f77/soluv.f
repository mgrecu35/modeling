
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine soluv (m,n,ndim,oh,dp,taucld,reff,taual,
     *                  csm,rsuvbm,flx,fdirpar,fdifpar)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1,nband
      PARAMETER(NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NXX-2,NBB=1,MM=NXI/NBB/2,N1=1)
      parameter (nband=8)
      real rflux(mm,7)
      common/radflux/rflux
c************************************************************************
c  compute solar fluxes in the uv+visible region. the spectrum is
c  grouped into 8 bands:
c  
c              Band     Micrometer
c
c       UV-C    1.     .175 - .225
c               2.     .225 - .245
c                      .260 - .280
c               3.     .245 - .260
c
c       UV-B    4.     .280 - .295
c               5.     .295 - .310
c               6.     .310 - .320
c      
c       UV-A    7.     .320 - .400
c      
c       PAR     8.     .400 - .700
c
c----- Input parameters:                            units      size
c
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c           meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   layer ozone content (oh)                      (cm-atm)stp m*n*np
c   layer pressure thickness (dp)                    mb       m*n*np
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   aerosol optical thickness (taual)                n/d      m*ndim*np
c   cosecant of the solar zenith angle (csm)         n/d      m*n
c   uv+par surface albedo for beam                 fraction   m*ndim
c           radiation (rsuvbm)
c
c----- output (updated) parameters:
c
c   all-sky net downward flux (flx)                fraction   m*ndim*(np+1)
c   all-sky direct downward par flux at
c          the surface (fdirpar)                   fraction   m*ndim
c   all-sky diffuse downward par flux at
c          the surface (fdifpar)                   fraction   m*ndim
c
c----- note: the following parameters must be specified by users:
c
c   aerosol single scattering albedo (ssaal)         n/d        1
c   aerosol asymmetry factor (asyal)                 n/d        1
c
*
c***********************************************************************

c-----input parameters

      integer m,n,ndim
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2)
      real    oh(m,n,np),dp(m,n,np),taual(m,ndim,np)
      real    rsuvbm(m,ndim),csm(m,n)

c-----output (updated) parameter

      real    flx(m,ndim,np+1)
      real    fdirpar(m,ndim),fdifpar(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----static parameters

c      integer nband

      real    hk(nband),xk(nband),ry(nband)
      real    asyal(nband),ssaal(nband),aig(3),awg(3)

c-----temporary array

      integer i,j,k,ib,mmm
      real    taurs,tauoz,tausto,ssatau,asysto,tauto,ssato,asyto
      real    taux,reff1,reff2,g1,g2,asycl(MM,N1,np)
      real    td(MM,N1,np+1),rr(MM,N1,np+1),tt(MM,N1,np+1),
     $        rs(MM,N1,np+1),ts(MM,N1,np+1)
      real    flxdn(MM,N1,np+1),fdndir(MM,N1),fdndif(MM,N1)
      real    asyclt(MM,N1)
      real   rr1t(MM,N1),tt1t(MM,N1),td1t(MM,N1),rs1t(MM,N1),ts1t(MM,N1)

c-----hk is the fractional extra-terrestrial solar flux.
c     the sum of hk is 0.47074.

      data hk/.00057, .00367, .00083, .00417,
     *        .00600, .00556, .05913, .39081/

c-----xk is the ozone absorption coefficient. unit: /(cm-atm)stp

      data xk /30.47, 187.2,  301.9,   42.83,
     *         7.09,  1.25,   0.0345,  0.0539/

c-----ry is the extinction coefficient for Rayleigh scattering.
c     unit: /mb.

      data ry /.00604, .00170, .00222, .00132,
     *         .00107, .00091, .00055, .00012/

c-----aerosol single-scattering albedo and asymmetry factor

      data ssaal/0.999,0.999,0.999,0.999,0.999,0.999,0.999,0.999/
      data asyal/0.850,0.850,0.850,0.850,0.850,0.850,0.850,0.850/

c-----coefficients for computing the asymmetry factor of ice clouds
c     from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

      data aig/.74625000,.00105410,-.00000264/

c-----coefficients for computing the asymmetry factor of liquid
c     clouds from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

      data awg/.82562000,.00529000,-.00014866/

      save

c-----initialize surface reflectances and transmittances
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
       do i= 1, mmm                    
         rr(i,j,np+1)=rsuvbm(i,j)
         rs(i,j,np+1)=rsuvbm(i,j)
         td(i,j,np+1)=0.0
         tt(i,j,np+1)=0.0
         ts(i,j,np+1)=0.0
       enddo
      enddo

c-----compute cloud asymmetry factor for a mixture of
c     liquid and ice particles.  unit of reff is micrometers.

      do k= 1, np

       do j= 1, n
        do i= 1, mmm

           asyclt(i,j)=1.0

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
          if (taux.gt.0.05) then

           reff1=min(reff(i,j,k,1),130.)
           reff2=min(reff(i,j,k,2),20.0)

           g1=(aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucld(i,j,k,1)
           g2=(awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucld(i,j,k,2)
           asyclt(i,j)=(g1+g2)/taux

          endif

        enddo
       enddo

       do j=1,n
        do i=1,mmm
           asycl(i,j,k)=asyclt(i,j)
        enddo
       enddo

      enddo
            
c-----integration over spectral bands

      do 100 ib=1,nband

       do 300 k= 1, np

        do j= 1, n
         do i= 1, mmm

c-----compute ozone and rayleigh optical thicknesses

          taurs=ry(ib)*dp(i,j,k)
          tauoz=xk(ib)*oh(i,j,k)
 
c-----compute total optical thickness, single scattering albedo,
c     and asymmetry factor

          tausto=taurs+tauoz+taual(i,j,k)+1.0e-8
          ssatau=ssaal(ib)*taual(i,j,k)+taurs
          asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)

c-----compute reflectance and transmittance

            taux=taucld(i,j,k,1)+taucld(i,j,k,2)
           tauto=tausto+taux
           ssato=(ssatau+taux)/tauto+1.0e-8
           ssato=min(ssato,0.999999)
           asyto=(asysto+asycl(i,j,k)*taux)/(ssato*tauto)

           call deledd (tauto,ssato,asyto,csm(i,j), 
     *                  rr1t(i,j),tt1t(i,j),td1t(i,j))

           call sagpol (tauto,ssato,asyto,rs1t(i,j),ts1t(i,j))

         enddo
        enddo

        do j=1,n
         do i=1,mmm
            rr(i,j,k)=rr1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            tt(i,j,k)=tt1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            td(i,j,k)=td1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            rs(i,j,k)=rs1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            ts(i,j,k)=ts1t(i,j)
         enddo
        enddo

 300  continue

c-----flux calculations

       call adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

       do k= 1, np+1
        do j= 1, n
         do i= 1, mmm
          flx(i,j,k)=flx(i,j,k)+flxdn(i,j,k)*hk(ib)

C --------- GCSS ---------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)+hk(ib)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)+rflux(i,2)*hk(ib)
         endif
C -------------------------------------------------------

         enddo
        enddo
       enddo

       if(ib.eq.8) then
         do j=1,n
          do i=1,mmm
           fdirpar(i,j)=fdndir(i,j)*hk(ib)
           fdifpar(i,j)=fdndif(i,j)*hk(ib)
         enddo
        enddo
       endif

 100  continue

      return
      end
