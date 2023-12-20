

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wvkdis(ib,m,n,np,k,h2oexp,conexp,th2o,tcon,tran)
c**********************************************************************
c   compute water vapor transmittance between levels k1 and k2 for
c   m x n soundings using the k-distribution method.
c
c   computations follow eqs. (34), (46), (50) and (52).
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of levels (np)
c  current level (k)
c  exponentials for line absorption (h2oexp) 
c  exponentials for continuum absorption (conexp) 
c
c---- updated parameters
c  transmittance between levels k1 and k2 due to
c    water vapor line absorption (th2o)
c  transmittance between levels k1 and k2 due to
c    water vapor continuum absorption (tcon)
c  total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,k

c---- input parameters ------

      real    conexp(m,n,np,3),h2oexp(m,n,np,6)

c---- updated parameters -----

      real    th2o(m,n,6),tcon(m,n,3),tran(m,n)

      integer iradave
      common/iptionr/ iradAVE

c---- static data -----
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------

      integer ne(8)
      real    fkw(6,8),gkw(6,3)

c---- temporary arrays -----
      real    trnth2o

c-----fkw is the planck-weighted k-distribution function due to h2o
c     line absorption given in table 4 of Chou and Suarez (1995).
c     the k-distribution function for the third band, fkw(*,3), is not used
 
      data fkw / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,
     2           0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,
     3           6*1.00,
     4           0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,
     5           0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,
     6           0.1846,0.2732,0.2353,0.1613,0.1146,0.0310,
     7           0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,
     8           0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/

c-----gkw is the planck-weighted k-distribution function due to h2o
c     line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
c     of band 3 given in table 7.  Note that the order of the sub-bands
c     is reversed.

      data gkw/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,
     2           0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,
     3           0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/

 
 
c-----ne is the number of terms used in each band to compute water vapor
c     continuum transmittance (table 6).
 
      data ne /0,0,3,1,1,1,0,0/


c-----tco2 are the six exp factors between levels k1 and k2 
c     tran is the updated total transmittance between levels k1 and k2


c-----th2o is the 6 exp factors between levels k1 and k2 due to
c     h2o line absorption. 

c-----tcon is the 3 exp factors between levels k1 and k2 due to
c     h2o continuum absorption.

c-----trnth2o is the total transmittance between levels k1 and k2 due
c     to both line and continuum absorption computed from eq. (52).


      integer i,j,mmm

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
         do j=1,n
          do i=1,mmm
           th2o(i,j,1) = th2o(i,j,1)*h2oexp(i,j,k,1)
           th2o(i,j,2) = th2o(i,j,2)*h2oexp(i,j,k,2)
           th2o(i,j,3) = th2o(i,j,3)*h2oexp(i,j,k,3)
           th2o(i,j,4) = th2o(i,j,4)*h2oexp(i,j,k,4)
           th2o(i,j,5) = th2o(i,j,5)*h2oexp(i,j,k,5)
           th2o(i,j,6) = th2o(i,j,6)*h2oexp(i,j,k,6)
          enddo
         enddo


      if (ne(ib).eq.0) then


         do j=1,n
          do i=1,mmm

           trnth2o      =(fkw(1,ib)*th2o(i,j,1)
     *                  + fkw(2,ib)*th2o(i,j,2)
     *                  + fkw(3,ib)*th2o(i,j,3)
     *                  + fkw(4,ib)*th2o(i,j,4)
     *                  + fkw(5,ib)*th2o(i,j,5)
     *                  + fkw(6,ib)*th2o(i,j,6))

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      elseif (ne(ib).eq.1) then


         do j=1,n
          do i=1,mmm

           tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)

           trnth2o      =(fkw(1,ib)*th2o(i,j,1)
     *                  + fkw(2,ib)*th2o(i,j,2)
     *                  + fkw(3,ib)*th2o(i,j,3)
     *                  + fkw(4,ib)*th2o(i,j,4)
     *                  + fkw(5,ib)*th2o(i,j,5)
     *                  + fkw(6,ib)*th2o(i,j,6))*tcon(i,j,1)

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      else

         do j=1,n
          do i=1,mmm

           tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)
           tcon(i,j,2)= tcon(i,j,2)*conexp(i,j,k,2)
           tcon(i,j,3)= tcon(i,j,3)*conexp(i,j,k,3)

           trnth2o      = (  gkw(1,1)*th2o(i,j,1)
     *                     + gkw(2,1)*th2o(i,j,2)
     *                     + gkw(3,1)*th2o(i,j,3)
     *                     + gkw(4,1)*th2o(i,j,4)
     *                     + gkw(5,1)*th2o(i,j,5)
     *                     + gkw(6,1)*th2o(i,j,6) ) * tcon(i,j,1)
     *                  + (  gkw(1,2)*th2o(i,j,1)
     *                     + gkw(2,2)*th2o(i,j,2)
     *                     + gkw(3,2)*th2o(i,j,3)
     *                     + gkw(4,2)*th2o(i,j,4)
     *                     + gkw(5,2)*th2o(i,j,5)
     *                     + gkw(6,2)*th2o(i,j,6) ) * tcon(i,j,2)
     *                  + (  gkw(1,3)*th2o(i,j,1)
     *                     + gkw(2,3)*th2o(i,j,2)
     *                     + gkw(3,3)*th2o(i,j,3)
     *                     + gkw(4,3)*th2o(i,j,4)
     *                     + gkw(5,3)*th2o(i,j,5)
     *                     + gkw(6,3)*th2o(i,j,6) ) * tcon(i,j,3)

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      endif



      return
      end
