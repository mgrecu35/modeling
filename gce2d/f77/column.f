

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine column (m,n,np,pa,dt,sabs0,sabs,spre,stem)
c
c**************************************************************************
c-----compute column-integrated (from top of the model atmosphere)
c     absorber amount (sabs), absorber-weighted pressure (spre) and
c     temperature (stem).
c     computations of spre and stem follows eqs. (37) and (38).
c
c--- input parameters
c   number of soundings in zonal direction (m)
c   number of soundings in meridional direction (n)
c   number of atmospheric layers (np)
c   layer pressure (pa)
c   layer temperature minus 250K (dt)
c   layer absorber amount (sabs0)
c
c--- output parameters
c   column-integrated absorber amount (sabs)
c   column absorber-weighted pressure (spre)
c   column absorber-weighted temperature (stem)
c
c--- units of pa and dt are mb and k, respectively.
c    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
c**************************************************************************
      implicit none
      integer m,n,np

c---- input parameters -----

      real    pa(m,n,np),dt(m,n,np),sabs0(m,n,np)

c---- output parameters -----

      real    sabs(m,n,np+1),spre(m,n,np+1),stem(m,n,np+1)

      integer iradave
      common/iptionr/ iradave

      integer i,j,k,mmm

      save

c*********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
        do j=1,n
         do i=1,mmm
          sabs(i,j,1)=0.0
          spre(i,j,1)=0.0
          stem(i,j,1)=0.0
         enddo
        enddo

        do k=1,np
         do j=1,n
          do i=1,mmm
           sabs(i,j,k+1)=sabs(i,j,k)+sabs0(i,j,k)
           spre(i,j,k+1)=spre(i,j,k)+pa(i,j,k)*sabs0(i,j,k)
           stem(i,j,k+1)=stem(i,j,k)+dt(i,j,k)*sabs0(i,j,k)
          enddo
         enddo
        enddo

       return
       end
