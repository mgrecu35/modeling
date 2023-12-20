cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine var_ave (t,tave)
      parameter (NZ=43)
c***********************************************************************
c*** this program computes the quantities displayed in fig.6 of      ***
c*** sui et al. (1994) jas v51. 711-728                              ***
c***-----------------------------------------------------------------***
c*** input variables:                                                ***
c***       kk........which level to display the thee field           ***
c***       t=qv.........density field                                ***
c***-----------------------------------------------------------------***
c*** output variables:                                               ***
c***       tave .....domain averaged and density weighted water vapor***
c***********************************************************************
      dimension t(nz)
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     1   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     2   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
c compute vertically averaged field
      tave  =0.0
      do 500 k=2,nz-1
        aa    =dz*rho(k)/am(k)
        tave  =tave  +aa*t(k)
 500  continue
      tave=10.*tave
      print*,'denq (mm)=',tave
      return
      end
