


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine co2kdis(m,n,np,k,co2exp,tco2,tran)
c**********************************************************************
c   compute co2 transmittances between levels k1 and k2 for m x n soundings
c   using the k-distribution method with linear pressure scaling.
c
c   computations follow eq. (34).
c
c---- input parameters
c   number of grid intervals in zonal direction (m)
c   number of grid intervals in meridional direction (n)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,n,np,k

c---- input parameters -----

      real    co2exp(m,n,np,6,2)

c---- updated parameters -----

      real    tco2(m,n,6,2),tran(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE
c---- static data -----

      real    gkc(6,2)

c---- temporary arrays -----

      real    xc

c-----gkc is the planck-weighted co2 k-distribution function 
c     in the band-wing and band-center regions given in table 7.
c     for computing efficiency, sub-bands 3a and 3c are combined.

      data gkc/  0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,
     2           0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/

c-----tco2 is the 6 exp factors between levels k1 and k2. 
c     xc is the total co2 transmittance given by eq. (53).

      integer i,j,mmm

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
         do j=1,n
          do i=1,mmm

c-----band-wings

           tco2(i,j,1,1)=tco2(i,j,1,1)*co2exp(i,j,k,1,1)
           xc=             gkc(1,1)*tco2(i,j,1,1)

           tco2(i,j,2,1)=tco2(i,j,2,1)*co2exp(i,j,k,2,1)
           xc=xc+gkc(2,1)*tco2(i,j,2,1)

           tco2(i,j,3,1)=tco2(i,j,3,1)*co2exp(i,j,k,3,1)
           xc=xc+gkc(3,1)*tco2(i,j,3,1)

           tco2(i,j,4,1)=tco2(i,j,4,1)*co2exp(i,j,k,4,1)
           xc=xc+gkc(4,1)*tco2(i,j,4,1)

           tco2(i,j,5,1)=tco2(i,j,5,1)*co2exp(i,j,k,5,1)
           xc=xc+gkc(5,1)*tco2(i,j,5,1)

           tco2(i,j,6,1)=tco2(i,j,6,1)*co2exp(i,j,k,6,1)
           xc=xc+gkc(6,1)*tco2(i,j,6,1)

c-----band-center region

           tco2(i,j,1,2)=tco2(i,j,1,2)*co2exp(i,j,k,1,2)
           xc=xc+gkc(1,2)*tco2(i,j,1,2)

           tco2(i,j,2,2)=tco2(i,j,2,2)*co2exp(i,j,k,2,2)
           xc=xc+gkc(2,2)*tco2(i,j,2,2)

           tco2(i,j,3,2)=tco2(i,j,3,2)*co2exp(i,j,k,3,2)
           xc=xc+gkc(3,2)*tco2(i,j,3,2)

           tco2(i,j,4,2)=tco2(i,j,4,2)*co2exp(i,j,k,4,2)
           xc=xc+gkc(4,2)*tco2(i,j,4,2)

           tco2(i,j,5,2)=tco2(i,j,5,2)*co2exp(i,j,k,5,2)
           xc=xc+gkc(5,2)*tco2(i,j,5,2)

           tco2(i,j,6,2)=tco2(i,j,6,2)*co2exp(i,j,k,6,2)
           xc=xc+gkc(6,2)*tco2(i,j,6,2)

           tran(i,j)=tran(i,j)*xc

          enddo
         enddo

      return
      end
