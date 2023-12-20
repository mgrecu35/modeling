

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)
c**********************************************************************
c   compute exponentials for water vapor line absorption
c   in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer water vapor amount for line absorption (dh2o) 
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer  (h2oexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np

c---- input parameters ------

      real    dh2o(m,n,np),pa(m,n,np),dt(m,n,np)

c---- output parameters -----

      real    h2oexp(m,n,np,6)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      integer mw(8)
      real    xkw(8),aw(8),bw(8)

c---- temporary arrays -----

      real    xh


c-----xkw  are the absorption coefficients for the first
c     k-distribution function due to water vapor line absorption
c     (tables 4 and 7).  units are cm**2/g    
 
      data xkw / 29.55  , 4.167e-1, 1.328e-2, 5.250e-4,
     *            5.25e-4, 2.340e-3, 1.320e-0, 5.250e-4/
 
c-----mw are the ratios between neighboring absorption coefficients
c     for water vapor line absorption (tables 4 and 7).
 
      data mw /6,6,8,6,6,8,6,16/

c-----aw and bw (table 3) are the coefficients for temperature scaling
c     in eq. (25).
 
      data aw/ 0.0021, 0.0140, 0.0167, 0.0302,
     *         0.0307, 0.0154, 0.0008, 0.0096/
      data bw/ -1.01e-5, 5.57e-5, 8.54e-5, 2.96e-4,
     *          2.86e-4, 7.53e-5,-3.52e-6, 1.64e-5/

      integer i,j,k,ik,mmm

      save

c**********************************************************************
c    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
c    and bw.  therefore, h2oexp for these sub-bands are identical.
c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
 
        do k=1,np
         do j=1,n
          do i=1,mmm

c-----xh is   the scaled water vapor amount for line absorption
c     computed from (27).
 
           xh = dh2o(i,j,k)*(pa(i,j,k)*0.002)
     1        * ( 1.+(aw(ib)+bw(ib)* dt(i,j,k))*dt(i,j,k) )

c-----h2oexp is the water vapor transmittance of the layer (k2-1)
c     due to line absorption

           h2oexp(i,j,k,1) = exp(-xh*xkw(ib))

          enddo
         enddo
        enddo

        do ik=2,6

         if(mw(ib).eq.6) then

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             h2oexp(i,j,k,ik) = xh*xh*xh
            enddo
           enddo
          enddo

        elseif(mw(ib).eq.8) then

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             xh = xh*xh
             h2oexp(i,j,k,ik) = xh*xh
            enddo
           enddo
          enddo

        else

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             xh = xh*xh
             xh = xh*xh
             h2oexp(i,j,k,ik) = xh*xh
            enddo
           enddo
          enddo

        endif
       enddo


      return
      end
