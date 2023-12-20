

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine co2exps(m,n,np,dco2,pa,dt,co2exp)
c       
c**********************************************************************
c   compute co2 exponentials for individual layers.
c
c---- input parameters
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer co2 amount (dco2)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer (co2exp)
c**********************************************************************
      implicit none
      integer m,n,np

c---- input parameters -----

      real    dco2(m,n,np),pa(m,n,np),dt(m,n,np)

c---- output parameters -----

      real    co2exp(m,n,np,6,2)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      real    xkc(2),ac(2),bc(2),pm(2),prc(2)

c---- temporary arrays -----

      real    xc

c-----xkc is the absorption coefficients for the
c     first k-distribution function due to co2 (table 7).
c     units are 1/(cm-atm)stp.
 
      data xkc/2.656e-5,2.656e-3/
 
c-----parameters (table 3) for computing the scaled co2 amount
c     using (27).

      data prc/  300.0,   30.0/
      data pm /    0.5,   0.85/
      data ac / 0.0182, 0.0042/
      data bc /1.07e-4,2.00e-5/
 
      integer i,j,k,mmm

      save
c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
        do k=1,np
         do j=1,n
          do i=1,mmm

c-----compute the scaled co2 amount from eq. (27) for band-wings
c     (sub-bands 3a and 3c).

           xc = dco2(i,j,k)*(pa(i,j,k)/prc(1))**pm(1)
     1             *(1.+(ac(1)+bc(1)*dt(i,j,k))*dt(i,j,k))

c-----six exponential by powers of 8 (table 7).

           co2exp(i,j,k,1,1)=exp(-xc*xkc(1))

           xc=co2exp(i,j,k,1,1)*co2exp(i,j,k,1,1)
           xc=xc*xc
           co2exp(i,j,k,2,1)=xc*xc

           xc=co2exp(i,j,k,2,1)*co2exp(i,j,k,2,1)
           xc=xc*xc
           co2exp(i,j,k,3,1)=xc*xc

           xc=co2exp(i,j,k,3,1)*co2exp(i,j,k,3,1)
           xc=xc*xc
           co2exp(i,j,k,4,1)=xc*xc

           xc=co2exp(i,j,k,4,1)*co2exp(i,j,k,4,1)
           xc=xc*xc
           co2exp(i,j,k,5,1)=xc*xc

           xc=co2exp(i,j,k,5,1)*co2exp(i,j,k,5,1)
           xc=xc*xc
           co2exp(i,j,k,6,1)=xc*xc

c-----compute the scaled co2 amount from eq. (27) for band-center
c     region (sub-band 3b).

           xc = dco2(i,j,k)*(pa(i,j,k)/prc(2))**pm(2)
     1             *(1.+(ac(2)+bc(2)*dt(i,j,k))*dt(i,j,k))

           co2exp(i,j,k,1,2)=exp(-xc*xkc(2))

           xc=co2exp(i,j,k,1,2)*co2exp(i,j,k,1,2)
           xc=xc*xc
           co2exp(i,j,k,2,2)=xc*xc

           xc=co2exp(i,j,k,2,2)*co2exp(i,j,k,2,2)
           xc=xc*xc
           co2exp(i,j,k,3,2)=xc*xc

           xc=co2exp(i,j,k,3,2)*co2exp(i,j,k,3,2)
           xc=xc*xc
           co2exp(i,j,k,4,2)=xc*xc

           xc=co2exp(i,j,k,4,2)*co2exp(i,j,k,4,2)
           xc=xc*xc
           co2exp(i,j,k,5,2)=xc*xc

           xc=co2exp(i,j,k,5,2)*co2exp(i,j,k,5,2)
           xc=xc*xc
           co2exp(i,j,k,6,2)=xc*xc

          enddo
         enddo
        enddo

      return
      end
