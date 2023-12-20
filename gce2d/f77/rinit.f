

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rinit (irs)
c     ******   initialize all variables to zeros   ******
      parameter (NX=514,NZ=43,NT=2880,lay=88,ITT=244,nnt=481)
      parameter (nxi=nx-2,nxm=nx-1,nzm=nz-1)
      parameter (nx3=3*nx,nx16=16*nx,nz2=2*nz,nz4=4*nz,nz9=9*nz)
      parameter (nz13=13*nz,nz15=15*nz,nz20=20*nz,nz23=23*nz,nz28=28*nz)
      parameter (nb=nt*nxi,nb1=5*3*nt,nb8=12*nz+5*nt)
      parameter (nb5=178,nb6=193,nb7=2*nz+6,nb9=48*nz*28+2*nz*4)
      parameter (nb12=2*nxm*nzm,nb13=13*nz*28,nb14=3*nz+1)
      parameter (nb15=2*nz*28)
c      parameter (nb16=5*nx*nz*4)
      parameter (nz40=15*nz+5*nz*itt,NB24=NZ*21*14)
      PARAMETER (NB25=21*NZ*14+2*NZ,NB26=16*NZ*14+5*NZ*14)
      parameter (nx2=2*nx,nz17=17*nz,nnt3=3*nnt,nnt7=7*nnt)
      parameter (nb27=2*nz*4*7,NZ41=2*ITT*NZ+2*NZ,nz6=6*nz)
      parameter (nz42=2*ITT*nz+3*nz,NZ43=4*ITT*NX,NZ44=4*NX+2*NZ)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz17)
      common/tqave/ tavet(nz2)
      common/tbudget/ avett(nnt3)
      common/bsfc/ tsfc2(nx2)
      common/tbudget1/ t_sfcq(nnt7)
      COMMON/BSTS20/ STS20(NB27)

      COMMON/Q1Q2Z/ Q1Q2(NZ41)
      COMMON/Q1Q2T/ AQ1ZT(NZ4)
      COMMON/UPPER1/ T_ADJUST(NZ42)
      COMMON/BLS99/ THS(NZ43)
      COMMON/DINRAD/ P00(NZ44)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    tcon(10)
      common/cont/ tcon
      real    tervr(9)
      common/rterv/ tervr
      real    trsnw(nb6)
      common/rsnw/ trsnw
      real    tb3cs(12)
      common/b3cs/ tb3cs
      real    rfa11(nb14)
      common/damp/ rfa11
      integer ivx(3),ivz(3)
      common/bxy/ ivx,ivz
      real    vx(8),vz(8)
      common/bx/ vx
      common/bz/ vz
      integer iib(3)
      common/bstart/ iib
      real    rrb(3),ao4(nb7),vb(14)
      common/rbstart/ rrb
      common/o4/ ao4
      common/bb/ vb
      real    v11(nx,nz),v12(nx,nz),v13(nx,nz),v14(nx,nz),
     $        v21(nx,nz),v22(nx,nz),v23(nx,nz),v24(nx,nz)
      common/b1t/ v11
      common/b1q/ v12
      common/b1c/ v13
      common/b1r/ v14
      common/b1u/ v21
      common/b1v/ v22
      common/b1w/ v23
      common/b1a/ v24
      real    vb1(nx,nz),vb2(nx,nz),vb3(nx,nz),vb4(nx,nz),
     $        vb5(nx,nz),vb6(nx,nz),vb7(nx,nz),vb8(nx,nz)
      common/b2t/ vb1
      common/b2q/ vb2
      common/b2c/ vb3
      common/b2r/ vb4
      common/b2u/ vb5
      common/b2v/ vb6
      common/b2w/ vb7
      common/b2a/ vb8
      real    qi (nx,nz),qs (nx,nz),qg (nx,nz),
     $        qi1(nx,nz),qs1(nx,nz),qg1(nx,nz)
      common/b1i/ qi
      common/b1s/ qs
      common/b1g/ qg
      common/b2i/ qi1
      common/b2s/ qs1
      common/b2g/ qg1
      real    v31(nx,nz),v32(nx,nz),v33(nx,nz),v34(nx,nz)
      common/bsat/ v31
      common/bsat1/ v32
      common/badv/ v33
      common/badv1/ v34
      real    sw1(nx,nz),sw2(nx,nz)
      common/slwave/ sw1,sw2
      real    vw(nx,nz),df0(nx,nz),trah(nx,nz),trav(nx,nz)
      common/bw/ vw,df0,trah,trav
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real    rst(nb),pcltop(nb),pclbot(nb),vw2(nb),vw3(nb1),VW22(NB)
c      real vw3(nb1)
      common/bw1/ rst,pcltop,pclbot 
      common/bw2/ vw2,VW22
      common/bw3/ vw3

      real    y4(nz9),y5(nz23),y6(nz13),x6(nx3),y8(nz15),y7(nz20),
     $        y9(nz20),ym(nz20),yb(nz40),ybu(nz28)
      common/b4/ y4
      common/b5/ y5
      common/b6/ y6,x6
      common/b8/ y8,y7
      common/b9/ y9,ym
      common/bb6/ yb
      common/bbtu/ ybu

      real    xa(nx16)
      common/ba/ xa

      real    s1(nb24),s2(nb25),s22(nb26)
      common/btt/ s1
      common/bcs/ s2
      common/bcss/ s22

      real    vrh1(nb8),vterv(6),vsnw(nb5)
      common/brh1/ vrh1
      common/bterv/ vterv
      common/bsnw/ vsnw

      real    sts3(nb9),sts4(nb13),stls1(nb15),rfft(nb12)
c      real sts5(nb16)
      common/bsts/ sts3
      common/bsts1/ sts4
c      common/bstsi/ sts5
      common/stls/ stls1
      common/rfft/ rfft

      real    sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)
      common/srflx/sfir,sfsw,shir,shsw,salpha,si
c
      COMMON/SUE/ PPRESS(NX,NZ)
      COMMON/SUE1/ TMODO(NX,NZ),QMODO(NX,NZ),RAINCO(NX,1)
      COMMON/SUE3/ T_BM(NB27)
      common/b66b/ s_dep(nz6)

c      COMMON/Q_BUGT/ Q1_G_H(NX,NZ),Q1_G_V(NX,NZ),
c     1               Q1A_G_H(NX,NZ),Q1A_G_V(NX,NZ),
c     2               Q1_D_H(NX,NZ),Q1_D_V(NX,NZ),
c     3               Q1A_D_H(NX,NZ),Q1A_D_V(NX,NZ),
c     4               Q2_G_H(NX,NZ),Q2_G_V(NX,NZ),
c     5               Q2A_G_H(NX,NZ),Q2A_G_V(NX,NZ),
c     6               Q2_D_H(NX,NZ),Q2_D_V(NX,NZ),
c     7               Q2A_D_H(NX,NZ),Q2A_D_V(NX,NZ),
c     8               Q1_HYD(NX,NZ),Q2_HYD(NX,NZ),
c     9               Q1A_HYD(NX,NZ),Q2A_HYD(NX,NZ),
c     9               Q1_RAD(NX,NZ),Q1A_RAD(NX,NZ),
c     9               IBUDSEC,RBUD

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,n

      save

ccccccccccccccccccc
      DO I=1,NX
        RAINCO(I,1)=0.
        DO K=1,NZ
          ppress(i,k)=0.
          tmodo(i,k)=0.
          qmodo(i,k)=0.
        ENDDO
      ENDDO
      DO I=1,NZ41
        Q1Q2(I)=0.
      ENDDO
      DO I=1,NZ4
        AQ1ZT(I)=0.
      ENDDO

      DO I=1,NZ42
        T_ADJUST(I)=0.
      ENDDO
      DO I=1,NZ43
        THS(I)=0.
      ENDDO
      DO I=1,NZ44
        P00(I)=0.
      ENDDO
ccccccccccccccccccc
       do 2 i=1,nb8
    2   vrh1(i)=0.0
       do 3 k=1,nz20
        y7(k)=0.
        y9(k)=0.
    3   ym(k)=0.
       do 4 k=1,nz15
    4   y8(k)=0.
      DO 74 I=1,NB24
  74    S1(I)=0.
      DO 76 I=1,NB25
  76    S2(I)=0.
      DO 78 I=1,NB26
  78    S22(I)=0.
       do 7 n=1,nb9
    7   sts3(n)=0.
       do 8 n=1,nb13
    8   sts4(n)=0.
       do n=1,nb27
        T_BM(n)=0.
        STS20(n)=0.
       enddo
       do 9 n=1,nb15
    9   stls1(n)=0.

      do k=1,nz17
        acoc(k)=0.
      enddo
      do k=1,nz2
        tavet(k)=0.
      enddo
      do ni=1,nnt3
        avett(ni)=0.
      enddo
      do ni=1,nnt7
        t_sfcq(ni)=0.
      enddo
      do i=1,nx2
        tsfc2(i)=0.
      enddo
c
      if(irs.eq.1) return
       do 11 i=1,nx3
   11   x6(i)=0.0
       do 12 i=1,3
        ivx(i)=0
        ivz(i)=0
        iib(i)=0
   12   rrb(i)=0.
       do 13 i=1,8
        vx(i)=0.
   13   vz(i)=0.
       do 14 i=1,14
   14   vb(i)=0.
       do 15 k=1,nz
       do 15 i=1,nx
        v11(i,k)=0.
        v12(i,k)=0.
        v13(i,k)=0.
        v14(i,k)=0.
        v21(i,k)=0.
        v22(i,k)=0.
        v23(i,k)=0.
        v24(i,k)=0.
        vb1(i,k)=0.
        vb2(i,k)=0.
        vb3(i,k)=0.
        vb4(i,k)=0.
        vb5(i,k)=0.
        vb6(i,k)=0.
        vb7(i,k)=0.
        vb8(i,k)=0.
        v31(i,k)=0.
        v32(i,k)=0.
        v33(i,k)=0.
        v34(i,k)=0.
        qi(i,k)=0.
        qg(i,k)=0.
        qs(i,k)=0.
        qi1(i,k)=0.
        qg1(i,k)=0.
        qs1(i,k)=0.
        sw1(i,k)=0.
        sw2(i,k)=0.

c        q1_g_h(i,k)=0.
c        q1_g_v(i,k)=0.
c        q1_d_h(i,k)=0.
c        q1_d_v(i,k)=0.
c        q2_g_h(i,k)=0.
c        q2_g_v(i,k)=0.
c        q2_d_h(i,k)=0.
c        q2_d_v(i,k)=0.
c        q1_hyd(i,k)=0.
c        q2_hyd(i,k)=0.
c        q1_rad(i,k)=0.
c        q1a_g_h(i,k)=0.
c        q1a_g_v(i,k)=0.
c        q1a_d_h(i,k)=0.
c        q1a_d_v(i,k)=0.
c        q2a_g_h(i,k)=0.
c        q2a_g_v(i,k)=0.
c        q2a_d_h(i,k)=0.
c        q2a_d_v(i,k)=0.
c        q1a_hyd(i,k)=0.
c        q2a_hyd(i,k)=0.
c        q1a_rad(i,k)=0.

   15   vw(i,k)=0.
       do 16 k=1,nz9
   16   y4(k)=0.
       do 17 k=1,nz23
   17   y5(k)=0.
       do 18 k=1,nz13
   18   y6(k)=0.
       do 19 i=1,nx16
   19   xa(i)=0.
       do 20 k=1,nz40
   20   yb(k)=0.
       do 22 n=1,nb
        rst(n)=0.
        pcltop(n)=0. 
        pclbot(n)=0.
        VW22(N)=0.
   22   vw2(n)=0.
       do 23 n=1,nb1
   23   vw3(n)=0.
       do 24 n=1,6
   24   vterv(n)=0.
       do 25 n=1,nb5
   25   vsnw(n)=0.
       do 26 n=1,nb7
   26   ao4(n)=0.
       do 27 n=1,nb12
   27   rfft(n)=0.
       do 28 n=1,10
   28   tcon(n)=0.
       do 29 n=1,9
   29   tervr(n)=0.
       do 30 n=1,nb6
   30   trsnw(n)=0.
       do 31 n=1,12
   31   tb3cs(n)=0.
       do 32 n=1,nb14
   32   rfa11(n)=0.
       do 33 i=1,4
       salpha(i)=0.
       si(i)=0.
       do 33 n=1,lay
        sfsw(n,i)=0.
        sfir(n,i)=0.
        shsw(n,i)=0.
   33   shir(n,i)=0.
      do k=1,nz6
         s_dep(k)=0.
      enddo
      return
      end
