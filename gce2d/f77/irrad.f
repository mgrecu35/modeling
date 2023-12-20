Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine irrad (m,n,ndim,taucl,ccld,pl,ta,wa,oa,co2,ts,
     *                  high,flx,flc,dfdts,st4)

      implicit none

      integer NXX,nz,nadd,nbb,nxi,mm,np,n1,nx,no,nc,nh,nt
      parameter (NXX=514,NZ=43,nadd=7,nbb=1)
      parameter (nxi=nxx-2)
      parameter (mm=nxi/nbb/2,np=nz-2+nadd)
      parameter (n1=1)
      parameter (nx=26,no=21,nc=24,nh=31,nt=7)
c
c---- input parameters ------
      integer m,n,ndim
      real    co2
      real    taucl(m,ndim,np),ccld(m,ndim,np),pl(m,ndim,np+1),
     $        ta(m,ndim,np),wa(m,ndim,np),oa(m,ndim,np),ts(m,ndim)

      logical high

c---- output parameters ------

      real    flx(m,ndim,np+1),flc(m,ndim,np+1),dfdts(m,ndim,np+1),
     $        st4(m,ndim)

      real rflux(mm,7)
      common/radflux/rflux
*
**********************************************************************
*
* This routine computes ir fluxes due to water vapor, co2, and o3.
*   Clouds in different layers are assumed randomly overlapped.
*  
* This is a vectorized code.  It computes fluxes simultaneously for
*   (m x n) soundings, which is a subset of (m x ndim) soundings.
*   In a global climate model, m and ndim correspond to the numbers of
*   grid boxes in the zonal and meridional directions, respectively.
*
* Detailed description of the radiation routine is given in
*   Chou and Suarez (1994).
*
* There are two options for computing cooling rate profiles.
*
*   if high = .true., transmission functions in the co2, o3, and the
*   three water vapor bands with strong absorption are computed using
*   table look-up.  cooling rates are computed accurately from the
*   surface up to 0.01 mb.
*   if high = .false., transmission functions are computed using the
*   k-distribution method with linear pressure scaling.  cooling rates
*   are not calculated accurately for pressures less than 20 mb.
*   the computation is faster with high=.false. than with high=.true.
*
* The IR spectrum is divided into eight bands:
*   
*   bnad     wavenumber (/cm)   absorber         method
*
*    1           0 - 340           h2o            K/T
*    2         340 - 540           h2o            K/T
*    3         540 - 800       h2o,cont,co2       K,S,K/T
*    4         800 - 980       h2o,cont           K,S
*    5         980 - 1100      h2o,cont,o3        K,S,T
*    6        1100 - 1380      h2o,cont           K,S
*    7        1380 - 1900          h2o            K/T
*    8        1900 - 3000          h2o            K 
*
* Note : "h2o" for h2o line absorption
*        "cont" for h2o continuum absorption
*        "K" for k-distribution method
*        "S" for one-parameter temperature scaling
*        "T" for table look-up
*
* The 15 micrometer region (540-800/cm) is further divided into
*   3 sub-bands :
*
*   subbnad   wavenumber (/cm)
*
*    1          540 - 620
*    2          620 - 720
*    3          720 - 800
*
*---- Input parameters                               units    size
*
*   number of soundings in zonal direction (m)        n/d      1
*   number of soundings in meridional direction (n)   n/d      1
*   maximum number of soundings in
*                 meridional direction (ndim)         n/d      1 
*   number of atmospheric layers (np)                 n/d      1
*   cloud optical thickness (taucl)                   n/d     m*ndim*np
*   cloud cover (ccld)                              fraction  m*ndim*np
*   level pressure (pl)                               mb      m*ndim*(np+1)
*   layer temperature (ta)                            k       m*ndim*np
*   layer specific humidity (wa)                      g/g     m*ndim*np
*   layer ozone mixing ratio by mass (oa)             g/g     m*ndim*np
*   surface temperature (ts)                          k       m*ndim  
*   co2 mixing ratio by volumn (co2)                  pppv     1
*   high                                                       1
*
* pre-computed tables used in table look-up for transmittance calculations:
*
*   c1 , c2, c3: for co2 (band 3)
*   o1 , o2, o3: for  o3 (band 5)
*   h11,h12,h13: for h2o (band 1)
*   h21,h22,h23: for h2o (band 2)
*   h71,h72,h73: for h2o (band 7)
*
*---- output parameters
*
*   net downward flux, all-sky   (flx)             w/m**2     m*ndim*(np+1)
*   net downward flux, clear-sky (flc)             w/m**2     m*ndim*(np+1)
*   sensitivity of net downward flux  
*       to surface temperature (dfdts)             w/m**2/k   m*ndim*(np+1)
*   emission by the surface (st4)                  w/m**2     m*ndim 
* 
* Notes: 
*
*   (1)  Water vapor continuum absorption is included in 540-1380 /cm.
*   (2)  Scattering by clouds is not included.
*   (3)  Clouds are assumed "gray" bodies.
*   (4)  The diffuse cloud transmission is computed to be exp(-1.66*taucl).
*   (5)  If there are no clouds, flx=flc.
*   (6)  plevel(1) is the pressure at the top of the model atmosphere, and
*        plevel(np+1) is the surface pressure.
*   (7)  Downward flux is positive, and upward flux is negative.
*   (8)  dfdts is always negative because upward flux is defined as negative.
*   (9)  For questions and coding errors, please contact with Ming-Dah Chou,
*        Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
*        Phone: 301-286-4012, Fax: 301-286-1759,
*        e-mail: chou@climate.gsfc.nasa.gov
*
c-----parameters defining the size of the pre-computed tables for transmittance
c     calculations using table look-up.
c
c     "nx" is the number of intervals in pressure
c     "no" is the number of intervals in o3 amount
c     "nc" is the number of intervals in co2 amount
c     "nh" is the number of intervals in h2o amount
c     "nt" is the number of copies to be made from the pre-computed
c          transmittance tables to reduce "memory-bank conflict"
c          in parallel machines and, hence, enhancing the speed of
c          computations using table look-up. 
c          If such advantage does not exist, "nt" can be set to 1.
c***************************************************************************

ctao
      integer iradave
      common/iptionr/ iradave
c
c---- static data -----

      real    cb(5,8)

c---- temporary arrays -----

      real pa(MM,N1,np),dt(MM,N1,np)
      real sh2o(MM,N1,np+1),swpre(MM,N1,np+1),swtem(MM,N1,np+1)
      real sco3(MM,N1,np+1),scopre(MM,N1,np+1),scotem(MM,N1,np+1)
      real dh2o(MM,N1,np),dcont(MM,N1,np),dco2(MM,N1,np),do3(MM,N1,np)
      real th2o(MM,N1,6),tcon(MM,N1,3),tco2(MM,N1,6,2)
      real h2oexp(MM,N1,np,6), conexp(MM,N1,np,3),co2exp(MM,N1,np,6,2)
      real clr(MM,N1,0:np+1),fclr(MM,N1)
      real blayer(MM,N1,0:np+1),dbs(MM,N1)
      real trant(MM,N1)
      real flxu(MM,N1,np+1),flxd(MM,N1,np+1)

      logical oznbnd
      logical co2bnd
      logical h2otbl
      logical conbnd

      real c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
      real o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
      real h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
      real h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
      real h71(nx,nh,nt),h72(nx,nh,nt),h73(nx,nh,nt)
      real dp,xx,w1,p1,dwe,dpe

c-----the following coefficients (table 2 of chou and suarez, 1995)
c     are for computing spectrally integtrated planck fluxes of
c     the 8 bands using eq. (22)
 
       data cb/
     1 -2.6844e-1,-8.8994e-2, 1.5676e-3,-2.9349e-6, 2.2233e-9,
     2  3.7315e+1,-7.4758e-1, 4.6151e-3,-6.3260e-6, 3.5647e-9,
     3  3.7187e+1,-3.9085e-1,-6.1072e-4, 1.4534e-5,-1.6863e-8,
     4 -4.1928e+1, 1.0027e+0,-8.5789e-3, 2.9199e-5,-2.5654e-8,
     5 -4.9163e+1, 9.8457e-1,-7.0968e-3, 2.0478e-5,-1.5514e-8,
     6 -1.0345e+2, 1.8636e+0,-1.1753e-2, 2.7864e-5,-1.1998e-8,
     7 -6.9233e+0,-1.5878e-1, 3.9160e-3,-2.4496e-5, 4.9301e-8,
     8  1.1483e+2,-2.2376e+0, 1.6394e-2,-5.3672e-5, 6.6456e-8/

c-----copy tables to enhance the speed of co2 (band 3), o3 (band5),
c     and h2o (bands 1, 2, and 7 only) transmission calculations
c     using table look-up.


      logical first
      data first /.true./
      integer i,j,k,ip,iw,it,ib,ik,iq,isb,k1,k2
      integer mmm


      include "h2o.tran3"
      include "co2.tran3"
      include "o3.tran3"

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1

      if (first) then

c-----tables co2 and h2o are only used with 'high' option

       if (high) then

        do iw=1,nh
         do ip=1,nx
          h11(ip,iw,1)=1.0-h11(ip,iw,1)
          h21(ip,iw,1)=1.0-h21(ip,iw,1)
          h71(ip,iw,1)=1.0-h71(ip,iw,1)
         enddo
        enddo

        do iw=1,nc
         do ip=1,nx
          c1(ip,iw,1)=1.0-c1(ip,iw,1)
         enddo
        enddo

c-----tables are replicated to avoid memory bank conflicts

        do it=2,nt
         do iw=1,nc
          do ip=1,nx
           c1 (ip,iw,it)= c1(ip,iw,1)
           c2 (ip,iw,it)= c2(ip,iw,1)
           c3 (ip,iw,it)= c3(ip,iw,1)
          enddo
         enddo
         do iw=1,nh
          do ip=1,nx
           h11(ip,iw,it)=h11(ip,iw,1)
           h12(ip,iw,it)=h12(ip,iw,1)
           h13(ip,iw,it)=h13(ip,iw,1)
           h21(ip,iw,it)=h21(ip,iw,1)
           h22(ip,iw,it)=h22(ip,iw,1)
           h23(ip,iw,it)=h23(ip,iw,1)
           h71(ip,iw,it)=h71(ip,iw,1)
           h72(ip,iw,it)=h72(ip,iw,1)
           h73(ip,iw,it)=h73(ip,iw,1)
          enddo
         enddo
        enddo

       endif

c-----always use table look-up for ozone transmittance

        do iw=1,no      
         do ip=1,nx
          o1(ip,iw,1)=1.0-o1(ip,iw,1)
         enddo
        enddo

        do it=2,nt
         do iw=1,no
          do ip=1,nx
           o1 (ip,iw,it)= o1(ip,iw,1)
           o2 (ip,iw,it)= o2(ip,iw,1)
           o3 (ip,iw,it)= o3(ip,iw,1)
          enddo
         enddo
        enddo

       first=.false.

      endif

c-----compute layer pressure (pa) and layer temperature minus 250K (dt)
 
      do k=1,np
       do j=1,n
        do i=1,mmm
         pa(i,j,k)=0.5*(pl(i,j,k)+pl(i,j,k+1))
         dt(i,j,k)=ta(i,j,k)-250.0
        enddo
       enddo
      enddo

c-----compute layer absorber amount

c     dh2o : water vapor amount (g/cm**2)
c     dcont: scaled water vapor amount for continuum absorption (g/cm**2)
c     dco2 : co2 amount (cm-atm)stp
c     do3  : o3 amount (cm-atm)stp
c     the factor 1.02 is equal to 1000/980
c     factors 789 and 476 are for unit conversion
c     the factor 0.001618 is equal to 1.02/(.622*1013.25) 
c     the factor 6.081 is equal to 1800/296

      do k=1,np
       do j=1,n
        do i=1,mmm

         dp           = pl(i,j,k+1)-pl(i,j,k)
         dh2o (i,j,k) = 1.02*wa(i,j,k)*dp+1.e-10
         dco2 (i,j,k) = 789.*co2*dp+1.e-10
         do3  (i,j,k) = 476.0*oa(i,j,k)*dp+1.e-10

c-----compute scaled water vapor amount for h2o continuum absorption
c     following eq. (43).

         xx=pa(i,j,k)*0.001618*wa(i,j,k)*wa(i,j,k)*dp
         dcont(i,j,k) = xx*exp(1800./ta(i,j,k)-6.081)+1.e-10

c-----compute effective cloud-free fraction, clr, for each layer.
c     the cloud diffuse transmittance is approximated by using a
c     diffusivity factor of 1.66.

         clr(i,j,k)=1.0-(ccld(i,j,k)*(1.-exp(-1.66*taucl(i,j,k))))

        enddo
       enddo
      enddo

c-----compute column-integrated h2o amoumt, h2o-weighted pressure
c     and temperature.  it follows eqs. (37) and (38).

       if (high) then

        call column(m,n,np,pa,dt,dh2o,sh2o,swpre,swtem)

       endif

c-----the surface (with an index np+1) is treated as a layer filled with
c     black clouds.

      do j=1,n
       do i=1,mmm
        clr(i,j,0)    = 1.0
        clr(i,j,np+1) = 0.0
        st4(i,j)      = 0.0
       enddo
      enddo

c-----initialize fluxes

      do k=1,np+1
       do j=1,n
        do i=1,mmm
         flx(i,j,k)  = 0.0
         flc(i,j,k)  = 0.0
         dfdts(i,j,k)= 0.0
         flxu(i,j,k) = 0.0
         flxd(i,j,k) = 0.0
        enddo
       enddo
      enddo

c-----integration over spectral bands

      do 1000 ib=1,8

c-----if h2otbl, compute h2o (line) transmittance using table look-up.
c     if conbnd, compute h2o (continuum) transmittance in bands 3, 4, 5 and 6.
c     if co2bnd, compute co2 transmittance in band 3.
c     if oznbnd, compute  o3 transmittance in band 5.
       h2otbl=high.and.(ib.eq.1.or.ib.eq.2.or.ib.eq.7)
       conbnd=ib.ge.3.and.ib.le.6
       co2bnd=ib.eq.3
       oznbnd=ib.eq.5

c-----blayer is the spectrally integrated planck flux of the mean layer
c     temperature derived from eq. (22)
c     the fitting for the planck flux is valid in the range 160-345 K.

       do k=1,np
        do j=1,n
         do i=1,mmm
          blayer(i,j,k)=ta(i,j,k)*(ta(i,j,k)*(ta(i,j,k)
     *                 *(ta(i,j,k)*cb(5,ib)+cb(4,ib))+cb(3,ib))
     *                 +cb(2,ib))+cb(1,ib)
         enddo
        enddo
       enddo

c-----the earth's surface, with an index "np+1", is treated as a layer

       do j=1,n
        do i=1,mmm
         blayer(i,j,0)   = 0.0
         blayer(i,j,np+1)=ts(i,j)*(ts(i,j)*(ts(i,j)
     *                   *(ts(i,j)*cb(5,ib)+cb(4,ib))+cb(3,ib))
     *                   +cb(2,ib))+cb(1,ib)

c-----dbs is the derivative of the surface planck flux with respect to
c     surface temperature (eq. 59).

         dbs(i,j)=ts(i,j)*(ts(i,j)*(ts(i,j)
     *           *4.*cb(5,ib)+3.*cb(4,ib))+2.*cb(3,ib))+cb(2,ib)

        enddo
       enddo

c-----compute column-integrated absorber amoumt, absorber-weighted
c     pressure and temperature for co2 (band 3) and o3 (band 5).
c     it follows eqs. (37) and (38).

c-----this is in the band loop to save storage

      if( high .and. co2bnd) then

        call column(m,n,np,pa,dt,dco2,sco3,scopre,scotem)

      endif

      if(oznbnd) then

        call column(m,n,np,pa,dt,do3,sco3,scopre,scotem)

      endif

c-----compute the exponential terms (eq. 32) at each layer for
c     water vapor line absorption when k-distribution is used

      if( .not. h2otbl) then

        call h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)

      endif

c-----compute the exponential terms (eq. 46) at each layer for
c     water vapor continuum absorption

      if( conbnd) then

        call conexps(ib,m,n,np,dcont,conexp)

      endif


c-----compute the  exponential terms (eq. 32) at each layer for
c     co2 absorption

      if( .not.high .and. co2bnd) then

        call co2exps(m,n,np,dco2,pa,dt,co2exp)

      endif

c-----compute transmittances for regions between levels k1 and k2
c     and update the fluxes at the two levels.

      do 2000 k1=1,np

c-----initialize fclr, th2o, tcon, and tco2

        do j=1,n
         do i=1,mmm
          fclr(i,j)=1.0
         enddo
        enddo

c-----for h2o line absorption

      if(.not. h2otbl) then
        do ik=1,6
         do j=1,n
          do i=1,mmm
           th2o(i,j,ik)=1.0
          enddo
         enddo
        enddo
      endif

c-----for h2o continuum absorption

       if (conbnd) then
         do iq=1,3
          do j=1,n
           do i=1,mmm
            tcon(i,j,iq)=1.0
           enddo
          enddo
         enddo
       endif

c-----for co2 absorption when using k-distribution method.
c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
c     are combined in computing the co2 transmittance.

       if (.not. high .and. co2bnd) then
         do isb=1,2
          do ik=1,6
           do j=1,n
            do i=1,mmm
             tco2(i,j,ik,isb)=1.0
            enddo
           enddo
          enddo
         enddo
       endif

c-----loop over the bottom level of the region (k2)

      do 3000 k2=k1+1,np+1

          do j=1,n
           do i=1,mmm
            trant(i,j)=1.0
           enddo
          enddo

       if(h2otbl) then

          w1=-8.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

c-----compute water vapor transmittance using table look-up

          if (ib.eq.1 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h11,h12,h13,trant)

          endif
          if (ib.eq.2 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h21,h22,h23,trant)

          endif
          if (ib.eq.7 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h71,h72,h73,trant)

          endif

       else

c-----compute water vapor transmittance using k-distribution.

        call wvkdis(ib,m,n,np,k2-1,h2oexp,conexp,th2o,tcon,trant)

       endif

       if(co2bnd) then

        if( high ) then

c-----compute co2 transmittance using table look-up method

          w1=-4.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

          call tablup(k1,k2,m,n,np,nx,nc,nt,sco3,scopre,scotem,
     *                w1,p1,dwe,dpe,c1,c2,c3,trant)

        else

c-----compute co2 transmittance using k-distribution method

          call co2kdis(m,n,np,k2-1,co2exp,tco2,trant)
        
        endif 

       endif 

c-----compute o3 transmittance using table look-up

       if (oznbnd) then

          w1=-6.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

          call tablup(k1,k2,m,n,np,nx,no,nt,sco3,scopre,scotem,
     *                w1,p1,dwe,dpe,o1,o2,o3,trant)

       endif

c-----fclr is the clear line-of-sight between levels k1 and k2.
c     in computing fclr, clouds are assumed randomly overlapped
c     using eq. (10).
 
      do j=1,n
       do i=1,mmm
        fclr(i,j) = fclr(i,j)*clr(i,j,k2-1)
       enddo
      enddo

c-----compute upward and downward fluxes


c-----add "boundary" terms to the net downward flux.
c     these are the first terms on the right-hand-side of
c     eqs. (56a) and (56b).
c     downward fluxes are positive.

      if (k2 .eq. k1+1) then
       do j=1,n
        do i=1,mmm
         flc(i,j,k1)=flc(i,j,k1)-blayer(i,j,k1)
         flc(i,j,k2)=flc(i,j,k2)+blayer(i,j,k1)
        enddo
       enddo
      endif

c-----add flux components involving the four layers above and below
c     the levels k1 and k2.  it follows eqs. (56a) and (56b).

      do j=1,n
       do i=1,mmm
        xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
        flc(i,j,k1) =flc(i,j,k1)+xx
        xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
        flc(i,j,k2) =flc(i,j,k2)+xx
       enddo
      enddo

c-----compute upward and downward fluxes for all-sky situation

      if (k2 .eq. k1+1) then
       do j=1,n
        do i=1,mmm
         flxu(i,j,k1)=flxu(i,j,k1)-blayer(i,j,k1)
         flxd(i,j,k2)=flxd(i,j,k2)+blayer(i,j,k1)
        enddo
       enddo
      endif

      do j=1,n
       do i=1,mmm
        xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
        flxu(i,j,k1) =flxu(i,j,k1)+xx*fclr(i,j)
        xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
        flxd(i,j,k2) =flxd(i,j,k2)+xx*fclr(i,j)
       enddo
      enddo


 3000 continue

c-----compute the partial derivative of fluxes with respect to
c     surface temperature (eq. 59).

      do j=1,n
       do i=1,mmm
        dfdts(i,j,k1) =dfdts(i,j,k1)-dbs(i,j)*trant(i,j)*fclr(i,j)
       enddo
      enddo
     
 2000 continue

c-----add contribution from the surface to the flux terms at the surface.

      do j=1,n
       do i=1,mmm
        dfdts(i,j,np+1) =dfdts(i,j,np+1)-dbs(i,j)
       enddo
      enddo

      do j=1,n
       do i=1,mmm
        flc(i,j,np+1)=flc(i,j,np+1)-blayer(i,j,np+1)
        flxu(i,j,np+1)=flxu(i,j,np+1)-blayer(i,j,np+1)
        st4(i,j)=st4(i,j)-blayer(i,j,np+1)
       enddo
      enddo


c     write(7,3211) ib, flxd(1,1,52),flxu(1,1,52)
c     write(7,3211) ib, flxd(1,1,np+1),flxu(1,1,np+1)
c     3211 format ('ib, fluxd, fluxu=', i3,2f12.3)

 1000 continue

      do k=1,np+1
       do j=1,n
        do i=1,mmm
         flx(i,j,k)=flxd(i,j,k)+flxu(i,j,k)
        enddo
       enddo
      enddo


C --------------------  GCSS Wrokshop -----------------------------------

      do i=1,mm
        rflux(i,2)=flxd(i,1,np+1)
        rflux(i,5)=flxu(i,1,1)
        rflux(i,7)=flxu(i,1,np+1)
      enddo

C -----------------------------------------------------------------------

      return
      end
