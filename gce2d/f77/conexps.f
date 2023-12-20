


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine conexps(ib,m,n,np,dcont,conexp)
c**********************************************************************
c   compute exponentials for continuum absorption in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer scaled water vapor amount for continuum absorption (dcont) 
c
c---- output parameters
c  1 or 3 exponentials for each layer (conexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np

c---- input parameters ------

      real    dcont(m,n,np)

c---- updated parameters -----

      real    conexp(m,n,np,3)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      real    xke(8)

c-----xke are the absorption coefficients for the first
c     k-distribution function due to water vapor continuum absorption
c     (table 6).  units are cm**2/g
 
      data xke /  0.00,   0.00,   27.40,   15.8,
     *            9.40,   7.75,     0.0,    0.0/
 
      integer i,j,k,iq,mmm

      save

c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1

        do k=1,np
         do j=1,n
          do i=1,mmm
           conexp(i,j,k,1) = exp(-dcont(i,j,k)*xke(ib))
          enddo
         enddo
        enddo

       if (ib .eq. 3) then

c-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
c     are, respectively, double and quadruple that for sub-band 3c (iq=1)
c     (table 6).

        do iq=2,3
         do k=1,np
          do j=1,n
           do i=1,mmm
            conexp(i,j,k,iq) = conexp(i,j,k,iq-1) *conexp(i,j,k,iq-1)
           enddo
          enddo
         enddo
        enddo

       endif

      return
      end
