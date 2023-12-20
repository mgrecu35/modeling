Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine var_avet (t,tave)
c      parameter (NX=514,NZ=43,ITT=244)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc*** input variables:                                                ***
cc***       t.........field                                           ***
cc***-----------------------------------------------------------------***
cc*** output variables:                                               ***
cc***       tave  ....domain averaged and density weighted temperature***
cc***********************************************************************
c      real t(nz),t1(nz)
c      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
c      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
c     1   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
c     2   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
c      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
c     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
c      save
cc  compute vertically averaged field
c      denave = 0.0
c      tave=0.0
c      do k=2,nz-1
c        aa=dz*rho(k)/am(k)
c        denave=denave+aa
c        tave  =tave  +aa*t1(k)
c      enddo
c      tave  =tave/denave
c        print *,'tave (k)= ',tave
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dom_ave (time,tave,denq,thee)
      parameter (NX=514,NZ=43)
c***********************************************************************
c*** This program computes the quantities displayed in Fig.6 of      ***
c*** Sui et al. (1994) JAS V51. 711-728                              ***
c***-----------------------------------------------------------------***
c*** Input variables:                                                ***
c***       kk........which level to display the thee field           ***
c***       t.........density field                                   ***
c***       the.......potential temperature field purtabation         ***
c***       qv........water vapor field purtabation                   ***
c***       p.........purtabation pressure field                      ***
c***       pi........base state pi field                             ***
c***       ta........base potential temperature                      ***
c***       qa........base q (kg/kg)                                  ***
c***       id........=1 when 2D model is used (y=1)                 ***
c***-----------------------------------------------------------------***
c*** Output variables:                                               ***
c***       tave  ....domain averaged and density weighted temperature***
c***       denq.....domain averaged and density weighted water vapor***
c***       thee......equivalent potential temperature                ***
c***********************************************************************
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      COMMON/BB/ DT,D2T,RIL2,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,PSFC
      COMMON/B2T/ DPT1(NX,NZ)
      COMMON/B2Q/ DQV1(NX,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),BA(NZ),BB(NZ),TA(NZ),
     1   QA(NZ),TA1(NZ),QA1(NZ),COEF(NZ),C1(NZ),C2(NZ),C3(NZ),AM(NZ),
     2   AM1(NZ),UB(NZ),VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ)
      COMMON/B6/ TLS(NZ),QLS(NZ),FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),
     1   ST(NZ),SV(NZ),SQ(NZ),SC(NZ),SE(NZ),SQA(NZ),RI(NX),AR(NX),RX(NX)
      common/tvertical/ denave
c local variables
      real    t(1000)
      save
      sx  = real(nx-2)
      do 100 k=2,nz-1
        t(k)=ta1(k)*PI(k)
 100  continue
c compute vertically averaged field
      denave=0.0
      tave  =0.0
      denq  =0.0
      dzt   =0.0
      do 500 k=2,nz-1
        aa    =dz*rho(k)/am(k)
        dzt   =dzt+dz/am(k)
        denave=denave+aa
        tave  =tave  +aa*t(k)
        denq  =denq  +aa*qa1(k)
 500  continue
      print *,'DZT=',DZT
      print *,'denave=',denave
      tave  =tave/denave
c     convert DENQ to mm
      denq  =10.0*denq
c compute THETAE at 950 mb
      kk=3
      thd=0.0
      do 600 i=2,nx-1
        t1=ta1(kk)+dpt1(i,kk)
        q1=qa1(kk)+dqv1(i,kk)
        thd=thd+t1*exp(1.0e3*2.4925263*q1/(t1*pi(kk)))
 600  continue
      thee=thd/sx
        print *,'time            = ',time
        print *,'temperature ave = ',tave
        print *,'water vapor ave = ',denq
        print *,'thret at 950 mb = ',thee
      return
      end
