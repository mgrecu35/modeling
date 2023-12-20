

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

      implicit none

      integer nx,nz,nadd,np,nxi,nbb,mm,n1
      PARAMETER(NX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NX-2)
      PARAMETER(NBB=1,MM=NXI/NBB/2)
      PARAMETER(N1=1)
      real rflux(mm,7)
      common/radflux/rflux
c*********************************************************************
c  compute upward and downward fluxes using a two-stream adding method
c  computations follow equations (3)-(5) of Chou (1992, JAS)
c  input parameters:
c     np: total number of layers 
c     rr:  reflection of a layer illuminated by beam radiation
c     tt:  diffuse transmission of a layer illuminated by beam radiation
c     td: direct beam tranmssion
c     ts: transmission of a layer illuminated by diffuse radiation
c     rs: reflection of a layer illuminated by diffuse radiation
c
c  output parameters:
c     flxdn:  net downward fluxes
c     fdndir: surface direct downward flux
c     fdndif: surface diffuse downward flux
c*********************************************************************

c-----input parameters

      integer m,n
      real    rr(m,n,np+1),tt(m,n,np+1),td(m,n,np+1),rs(m,n,np+1),
     $        ts(m,n,np+1)

c-----output parameters

      real    flxdn(m,n,np+1),fdndir(m,n),fdndif(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array

      integer i,j,k,mmm
      real    denm,xx,fupdif
      real    rssab(MM,N1,np+1),rabx(MM,N1,np+1),rsabx(MM,N1,np+1)
      real    tbab(MM,N1,np+1),tab(MM,N1,np+1)

      save

c-----layers are added one at a time, going down from the top layer,
c     tbab is the composite transmittance illuminated by beam radiation
c     tab is the composite diffuse transmittance illuminated by
c         beam radiation
c     rssab is the composite reflectance illuminated from below
c         by diffuse radiation
c     tab and rssab are computed from eqs. (4b) and (3b) of Chou
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
        do i= 1, mmm
          tbab(i,j,1)  = td(i,j,1)
          tab(i,j,1)   = tt(i,j,1)
          rssab(i,j,1) = rs(i,j,1)
        enddo
      enddo

      do k= 2, np
       do j= 1, n
        do i= 1, mmm
                 denm = ts(i,j,k)/( 1.-rssab(i,j,k-1)*rs(i,j,k) )
          tbab(i,j,k) = tbab(i,j,k-1)*td(i,j,k)
          tab(i,j,k)  = tbab(i,j,k-1)*tt(i,j,k)
     1                + (tbab(i,j,k-1)*rssab(i,j,k-1)
     2                * rr(i,j,k)+tab(i,j,k-1))*denm
          rssab(i,j,k)= rs(i,j,k)+ts(i,j,k)*rssab(i,j,k-1)*denm
        enddo
       enddo
      enddo

c-----layers are added one at a time, going up
c     rabx is the composite reflectance illuminated by beam radiation
c     rsabx is the composite reflectance illuminated from above
c         by diffuse radiation
c     rabx and rsabx are computed from eqs. (4a) and (3a) of Chou
 
      do j= 1, n
       do i= 1, mmm
         rabx(i,j,np+1) = rr(i,j,np+1)
         rsabx(i,j,np+1)= rs(i,j,np+1)
       enddo
      enddo

      do k= np, 1, -1
       do j= 1, n
        do i= 1, mmm
                 denm  = ts(i,j,k)/( 1.-rs(i,j,k)*rsabx(i,j,k+1) )
          rabx(i,j,k)  = rr(i,j,k)+(td(i,j,k)*rabx(i,j,k+1)
     *                 + tt(i,j,k)*rsabx(i,j,k+1))*denm
          rsabx(i,j,k) = rs(i,j,k)+ts(i,j,k)*rsabx(i,j,k+1)*denm
        enddo
       enddo
      enddo
 
c-----compute fluxes following eq (5) of Chou (1992)
 
c     fdndir is the direct  downward flux
c     fdndif is the diffuse downward flux
c     fupdif is the diffuse upward flux

      do k=2,np+1
       do j=1, n
        do i=1, mmm
                denm  = 1./(1.- rssab(i,j,k-1)*rsabx(i,j,k))
         fdndir(i,j)  = tbab(i,j,k-1)
                  xx  = tbab(i,j,k-1)*rabx(i,j,k)
         fdndif(i,j)  = (xx*rssab(i,j,k-1)+tab(i,j,k-1))*denm
              fupdif  = (xx+tab(i,j,k-1)*rsabx(i,j,k))*denm
         flxdn(i,j,k) = fdndir(i,j)+fdndif(i,j)-fupdif

C ------------- GCSS ---------------------------------------------------
         if(k.eq.(np+1))then
           rflux(i,2)=fdndir(i,j)+fdndif(i,j)
         endif
C ----------------------------------------------------------------------
         
        enddo
       enddo
      enddo

       do j=1, n
        do i=1, mmm
         flxdn(i,j,1) = 1.0-rabx(i,j,1)
        enddo
       enddo

      return
      end
