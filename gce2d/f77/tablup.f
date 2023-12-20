

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tablup(k1,k2,m,n,np,nx,nh,nt,sabs,spre,stem,w1,p1,
     *                  dwe,dpe,coef1,coef2,coef3,tran)
c**********************************************************************
c   compute water vapor, co2, and o3 transmittances between levels k1 and k2
c   using table look-up for m x n soundings.
c
c   Calculations follow Eq. (40) of Chou and Suarez (1995)
c
c---- input ---------------------
c  indices for pressure levels (k1 and k2)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of atmospheric layers (np)
c  number of pressure intervals in the table (nx)
c  number of absorber amount intervals in the table (nh)
c  number of tables copied (nt)
c  column-integrated absorber amount (sabs)
c  column absorber amount-weighted pressure (spre)
c  column absorber amount-weighted temperature (stem)
c  first value of absorber amount (log10) in the table (w1) 
c  first value of pressure (log10) in the table (p1) 
c  size of the interval of absorber amount (log10) in the table (dwe)
c  size of the interval of pressure (log10) in the table (dpe)
c  pre-computed coefficients (coef1, coef2, and coef3)
c
c---- updated ---------------------
c  transmittance (tran)
c
c  Note:
c   (1) units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.
c   (2) units of spre and stem are, respectively, mb and K.
c   (3) there are nt identical copies of the tables (coef1, coef2, and
c       coef3).  the prupose of using the multiple copies of tables is
c       to increase the speed in parallel (vectorized) computations.
C       if such advantage does not exist, nt can be set to 1.
c   
c**********************************************************************
      implicit none
      integer k1,k2,m,n,np,nx,nh,nt

c---- input parameters -----

      real    w1,p1,dwe,dpe
      real    sabs(m,n,np+1),spre(m,n,np+1),stem(m,n,np+1)
      real    coef1(nx,nh,nt),coef2(nx,nh,nt),coef3(nx,nh,nt)

c---- update parameter -----

      real tran(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- temporary variables -----

      integer i,j,mmm,iw,ip,nn
      real    x1,x2,x3,we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2

      save

c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
      do j=1,n
       do i=1,mmm

        nn=mod(i,nt)+1

        x1=sabs(i,j,k2)-sabs(i,j,k1)
        x2=(spre(i,j,k2)-spre(i,j,k1))/x1
        x3=(stem(i,j,k2)-stem(i,j,k1))/x1

        we=(log10(x1)-w1)/dwe
        pe=(log10(x2)-p1)/dpe

        we=max(we,w1-2.*dwe)
        pe=max(pe,p1)

        iw=int(we+1.5)
        ip=int(pe+1.5)

        iw=min(iw,nh-1)
        iw=max(iw, 2)

        ip=min(ip,nx-1)
        ip=max(ip, 1)

        fw=we-float(iw-1)
        fp=pe-float(ip-1)

c-----linear interpolation in pressure

        pa = coef1(ip,iw-1,nn)*(1.-fp)+coef1(ip+1,iw-1,nn)*fp
        pb = coef1(ip,iw,  nn)*(1.-fp)+coef1(ip+1,iw,  nn)*fp
        pc = coef1(ip,iw+1,nn)*(1.-fp)+coef1(ip+1,iw+1,nn)*fp

c-----quadratic interpolation in absorber amount for coef1

        ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)

c-----linear interpolation in absorber amount for coef2 and coef3

        ba = coef2(ip,iw,  nn)*(1.-fp)+coef2(ip+1,iw,  nn)*fp
        bb = coef2(ip,iw+1,nn)*(1.-fp)+coef2(ip+1,iw+1,nn)*fp
        t1 = ba*(1.-fw) + bb*fw

        ca = coef3(ip,iw,  nn)*(1.-fp)+coef3(ip+1,iw,  nn)*fp
        cb = coef3(ip,iw+1,nn)*(1.-fp)+coef3(ip+1,iw+1,nn)*fp
        t2 = ca*(1.-fw) + cb*fw

c-----update the total transmittance between levels k1 and k2

        tran(i,j)= (ax + (t1+t2*x3) * x3)*tran(i,j)

       enddo
      enddo

      return
      end
