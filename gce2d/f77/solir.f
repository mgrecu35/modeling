

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solir (m,n,ndim,wh,taucld,reff,taual,
     $                  csm,rsirbm,flx,fdirir,fdifir)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1,nk,nband
      PARAMETER(NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NXX-2)
      PARAMETER(NBB=1,MM=NXI/NBB/2)
      PARAMETER(N1=1)
      parameter (nk=10,nband=3)
      real rflux(mm,7)
      common/radflux/rflux
c************************************************************************
c  compute solar flux in the infrared region. The spectrum is divided
c   into three bands:
c
c          band   wavenumber(/cm)  wavelength (micron)
c           1       1000-4400         2.27-10.0
c           2       4400-8200         1.22-2.27
c           3       8200-14300        0.70-1.22
c
c----- Input parameters:                            units      size
c
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c          meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   layer water vapor content (wh)                 gm/cm^2    m*n*np
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   aerosol optical thickness (taual)                n/d      m*ndim*np
c   cosecant of the solar zenith angle (csm)         n/d      m*n
c   near ir surface albedo for beam                fraction   m*ndim
c                radiation (rsirbm)
c
c----- output (updated) parameters:
c
c   all-sky flux (downward-upward) (flx)           fraction   m*ndim*(np+1)
c   all-sky direct downward ir flux at
c          the surface (fdirir)                    fraction   m*ndim
c   all-sky diffuse downward ir flux at
c          the surface (fdifir)                    fraction   m*ndim
c
c----- note: the following parameters must be specified by users:
c   aerosol single scattering albedo (ssaal)         n/d      nband
c   aerosol asymmetry factor (asyal)                 n/d      nband
c
c*************************************************************************

c-----input parameters

      integer m,n,ndim
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2),rsirbm(m,ndim)
      real    wh(m,n,np),taual(m,ndim,np),csm(m,n)

c-----output (updated) parameters

      real    flx(m,ndim,np+1)
      real    fdirir(m,ndim),fdifir(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----static parameters

      real    xk(nk),hk(nband,nk),ssaal(nband),asyal(nband)
      real    aia(nband,3),awa(nband,3),aig(nband,3),awg(nband,3)

c-----temporary array

      integer ib,ik,i,j,k,mmm
      real    ssacl(MM,N1,np),asycl(MM,N1,np)
      real    rr(MM,N1,np+1),tt(MM,N1,np+1),td(MM,N1,np+1),
     $        rs(MM,N1,np+1),ts(MM,N1,np+1)
      real    flxdn(MM,N1,np+1),fdndir(MM,N1),fdndif(MM,N1)
      real    tauwv,tausto,ssatau,asysto,tauto,ssato,asyto
      real    taux,reff1,reff2,w1,w2,g1,g2
      real    ssaclt(MM,N1),asyclt(MM,N1)
      real   rr1t(MM,N1),tt1t(MM,N1),td1t(MM,N1),rs1t(MM,N1),ts1t(MM,N1)

c-----water vapor absorption coefficient for 10 k-intervals.
c     unit: cm^2/gm

      data xk/            
     1  0.0010, 0.0133, 0.0422, 0.1334, 0.4217,            
     2  1.334,  5.623,  31.62,  177.8,  1000.0/  

c-----water vapor k-distribution function,
c     the sum of hk is 0.52926. unit: fraction

      data hk/
     1 .01074,.08236,.20673,  .00360,.01157,.03497,
     2 .00411,.01133,.03011,  .00421,.01143,.02260,
     3 .00389,.01240,.01336,  .00326,.01258,.00696,
     4 .00499,.01381,.00441,  .00465,.00650,.00115,
     5 .00245,.00244,.00026,  .00145,.00094,.00000/

c-----aerosol single-scattering albedo and asymmetry factor

      data ssaal/0.999, 0.999, 0.999/
      data asyal/0.850, 0.850, 0.850/
 
c-----coefficients for computing the single scattering albedo of
c     ice clouds from ssa=1-(aia(*,1)+aia(*,2)*reff+aia(*,3)*reff**2)

      data aia/
     1  .08938331, .00215346,-.00000260,
     2  .00299387, .00073709, .00000746,
     3 -.00001038,-.00000134, .00000000/

c-----coefficients for computing the single scattering albedo of
c     liquid clouds from ssa=1-(awa(*,1)+awa(*,2)*reff+awa(*,3)*reff**2)

      data awa/
     1  .01209318,-.00019934, .00000007,
     2  .01784739, .00088757, .00000845,
     3 -.00036910,-.00000650,-.00000004/

c-----coefficients for computing the asymmetry factor of ice clouds
c     from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

      data aig/
     1  .84090400, .76098937, .74935228,
     2  .00126222, .00141864, .00119715,
     3 -.00000385,-.00000396,-.00000367/

c-----coefficients for computing the asymmetry factor of liquid clouds
c     from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

      data awg/
     1  .83530748, .74513197, .79375035,
     2  .00257181, .01370071, .00832441,
     3  .00005519,-.00038203,-.00023263/

      save

c-----initialize surface fluxes, reflectances, and transmittances
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
       do i= 1, mmm
         fdirir(i,j)=0.0
         fdifir(i,j)=0.0
         rr(i,j,np+1)=rsirbm(i,j)
         rs(i,j,np+1)=rsirbm(i,j)
         td(i,j,np+1)=0.0
         tt(i,j,np+1)=0.0
         ts(i,j,np+1)=0.0
       enddo
      enddo

c-----integration over spectral bands

      do 100 ib=1,nband

c-----compute cloud single scattering albedo and asymmetry factor
c     for a mixture of ice and liquid particles.

       do k= 1, np

        do j= 1, n
         do i= 1, mmm

           ssaclt(i,j)=1.0
           asyclt(i,j)=1.0

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
          if (taux.gt.0.05) then

           reff1=min(reff(i,j,k,1),130.)
           reff2=min(reff(i,j,k,2),20.0)

           w1=(1.-(aia(ib,1)+(aia(ib,2)+
     *         aia(ib,3)*reff1)*reff1))*taucld(i,j,k,1)
           w2=(1.-(awa(ib,1)+(awa(ib,2)+
     *         awa(ib,3)*reff2)*reff2))*taucld(i,j,k,2)
           ssaclt(i,j)=(w1+w2)/taux

           g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
           g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
           asyclt(i,j)=(g1+g2)/(w1+w2)

          endif

         enddo
        enddo

        do j=1,n
         do i=1,mmm
            ssacl(i,j,k)=ssaclt(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            asycl(i,j,k)=asyclt(i,j)
         enddo
        enddo

       enddo

c-----integration over the k-distribution function

         do 200 ik=1,nk

          do 300 k= 1, np

           do j= 1, n
            do i= 1, mmm

             tauwv=xk(ik)*wh(i,j,k)
 
c-----compute total optical thickness, single scattering albedo,
c     and asymmetry factor.
 
             tausto=tauwv+taual(i,j,k)+1.0e-8
             ssatau=ssaal(ib)*taual(i,j,k)
             asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)
 
c-----compute reflectance and transmittance

              taux=taucld(i,j,k,1)+taucld(i,j,k,2)
              tauto=tausto+taux
c             tauto=max(tausto+taux,1.0e-8)
              ssato=(ssatau+ssacl(i,j,k)*taux)/tauto+1.0e-8
c             ssato=max( (ssatau+ssacl(i,j,k)*taux)/tauto,1.0e-8)
              ssato=min(ssato,0.999999)
              asyto=(asysto+asycl(i,j,k)*ssacl(i,j,k)*taux)/
     *              (ssato*tauto)

              call deledd (tauto,ssato,asyto,csm(i,j), 
     *                     rr1t(i,j),tt1t(i,j),td1t(i,j))

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

c-----flux calculations at each level using the two-stream adding method
 
       call adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

       do k= 1, np+1
        do j= 1, n
         do i= 1, mmm
          flx(i,j,k) = flx(i,j,k)+flxdn(i,j,k)*hk(ib,ik)

C --------- GCSS ---------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)+hk(ib,ik)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)+rflux(i,2)*hk(ib,ik)
         endif
C -------------------------------------------------------

         enddo
        enddo
       enddo

       do j= 1, n
        do i= 1, mmm
          fdirir(i,j) = fdirir(i,j)+fdndir(i,j)*hk(ib,ik)
          fdifir(i,j) = fdifir(i,j)+fdndif(i,j)*hk(ib,ik)
        enddo
       enddo

  200 continue
  100 continue
 
      return
      end
