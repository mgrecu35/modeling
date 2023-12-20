

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadvuw (x,x1,u,w,umod,wmod)

CC    ****   COMPUTE ADVECTION OF DIFFERENT TRACERS

      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),u(nx,nz),w(nx,nz),umod(nx,nz),wmod(nx,nz)
      real    x1(nx,nz)

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,epss,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   t(nx),tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,t,tp,tm,q,qp,qm,y8,y9,y10

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,im,imm,ip,k,km,kmm,kp
      real    a18,ar1,arm,arp,arr
c      real    a0
      real    dtti,dttk,scc
c      real    d1t25

      save
c
       EPS=1.E-10
       SCC=1.0
       A18=1./8.*DT*SCC
c       D1T25=.25*DT*SCC
       DO 5 K=1,KMAX
       DO 5 I=1,IMAX
        UMOD(I,K)=0.
    5   WMOD(I,K)=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 100 K=2,KLES
         KP=K+1
         KM=K-1
c        A0=DZR(K)
        ARR=RRHO(K)
        AR1=RHO1(K)
        ARP=RHO1(KP)
       DO 10 I=2,ILES
        Y4(I)=(X(I,KP)+X(I-1,KP))
        Y5(I)=(X(I,K)+X(I-1,K))
   10   Y6(I)=(X(I,KM)+X(I-1,KM))
       DO 15 I=2,ILES
         IP=I+1
         IM=I-1
         IMM=I-2
         IF(I.EQ.2) IMM=IL2
        IF (Y5(I) .GE. EPS) THEN
           UMOD(I,K)=SCC*(ABS(U(I,K))-U(I,K)*U(I,K)*DXXT(I))
     1                  *(X(I,K)-X(IM,K))/(X(IM,K)+X(I,K)+EPS)
     2      -A18*AM(K)*RDZ*U(I,K)*ARR*(Y4(I)-Y6(I))/(Y4(I)+Y6(I)+EPS)
     3                *(AR1*(W(I,K)+W(IM,K))+ARP*(W(I,KP)+W(IM,KP)))
     4         -U(I,K)*(X(IMM,K)-X(IM,K)-X(I,K)+X(IP,K))/(3.*(X(IMM,K)
     5                  +X(IM,K)+X(I,K)+X(IP,K)+EPS))
c    6         -D1T25*RHO(K)*U(I,K)*(RDX*(U(IP,K)-U(IM,K))
c    7                        +A0*(ARP*(W(I,KP)+W(IM,KP))
c    8                            -AR1*(W(I,K)+W(IM,K))))
        ENDIF
   15  CONTINUE
       UMOD(IMAX,K)=UMOD(2,K)
       UMOD(1,K)=UMOD(ILES,K)
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 300 K=3,KLES
         KP=K+1
         KM=K-1
         KMM=K-2
c       A0=AM1(K)*RDZ*RRHO1(K)
       ARR=RRHO1(K)
       AR1=RHO(K)
       ARM=RHO(KM)
      DO 30 I=2,ILES
        Y4(I)=X(I,K)
        Y5(I)=X(I,K)+X(I,KM)
   30   Y6(I)=X(I,KM)
      DO 35 I=2,ILES
        IP=I+1
        IM=I-1
       IF (Y5(I) .GE. EPS) THEN
         WMOD(I,K)=SCC*(ABS(W(I,K))-W(I,K)*W(I,K)*DZZT(K))
     1                *(X(I,K)-X(I,KM))/(X(I,K)+X(I,KM)+EPS)
     2             -A18*RDX*ARR*W(I,K)
     3                 *(AR1*(U(I,K)+U(IP,K))+ARM*(U(I,KM)+U(IP,KM)))
     4  *(Y4(IP)-Y4(IM)+Y6(IP)-Y6(IM))/(Y4(IP)+Y4(IM)+Y6(IP)+Y6(IM)+EPS)
     5          -W(I,K)*(X(I,KMM)-X(I,KM)-X(I,K)+X(I,KP))/(3.*(X(I,KMM)
     6                  +X(I,KM)+X(I,K)+X(I,KP)+EPS))
c    7       -D1T25*RHO1(K)*W(I,K)*(RDX*ARR*(AR1*(U(IP,K)-U(I,K))
c    8                                       +ARM*(U(IP,KM)-U(I,KM)))
c    9                 +A0*(RHO1(KP)*W(I,KP)-RHO1(KM)*W(I,KM)))
       ENDIF
   35 CONTINUE
       WMOD(IMAX,K)=WMOD(2,K)
       WMOD(1,K)=WMOD(ILES,K)
  300 CONTINUE
      DO 900 I=1,IMAX
       UMOD(I,1)=UMOD(I,2)
       UMOD(I,KMAX)=UMOD(I,KLES)
       WMOD(I,1)=0.0
  900  WMOD(I,KMAX)=0.0
CCCCC   NON-OSCILLATORY OPTION (SMOLARKIEWICZ AND GRABOWSKI, 1990)
      DO 301 K=1,KMAX
      DO 301 I=1,IMAX
       U(I,K)=0.
  301  W(I,K)=0.
      DO 200 K=2,KLES
         KP=K+1
         KM=K-1
          DTTK=DT*AM(K)*RDZ
       DO 20 I=2,ILES
   20   Y5(I)=(X(I,K)+X(I-1,K))
       DO 22 I=2,ILES
         IP=I+1
         IM=I-1
         Y1(I)=0.
         Y2(I)=0.
        IF (Y5(I) .GE. EPS) THEN
         Y1(I)=MAX(X(IM,K),X(I,K),X(IP,K),X1(IM,K),X1(I,K),X1(IP,K))
         Y2(I)=MIN(X(IM,K),X(I,K),X(IP,K),X1(IM,K),X1(I,K),X1(IP,K))
        ENDIF
   22  CONTINUE
       DO 24 I=2,ILES
         IP=I+1
         IM=I-1
         Y3(I)=0.
         Y4(I)=0.
        IF (Y5(I) .GE. EPS) THEN
          DTTI=DT*RDX
        Y3(I)=(Y1(I)-X(I,K))/(DTTI*(MAX(UMOD(I,K),0.)*X(IM,K)
     1                             -MIN(UMOD(IP,K),0.)*X(IP,K))
     2                       +DTTK*(MAX(WMOD(I,K),0.)*X(I,KM)
     3                             -MIN(WMOD(I,KP),0.)*X(I,KP))+EPS)
        Y4(I)=(X(I,K)-Y2(I))/(DTTI*(MAX(UMOD(IP,K),0.)*X(I,K)
     1                             -MIN(UMOD(I,K),0.)*X(I,K))
     2                       +DTTK*(MAX(WMOD(I,KP),0.)*X(I,K)
     3                             -MIN(WMOD(I,K),0.)*X(I,K))+EPS)
        ENDIF
   24  CONTINUE
        Y3(1)=Y3(ILES)
        Y4(1)=Y4(ILES)
       DO 26 I=2,ILES
         IM=I-1
        IF (Y5(I) .GE. EPS) THEN
           U(I,K)=MIN(1.,Y4(IM),Y3(I))*MAX(0.,UMOD(I,K))
     1           +MIN(1.,Y4(I),Y3(IM))*MIN(0.,UMOD(I,K))
        ENDIF
   26  CONTINUE
       U(IMAX,K)=U(2,K)
       U(1,K)=U(ILES,K)
  200 CONTINUE
      DO 28 I=1,NX
       U(I,1)=U(I,2)
   28  U(I,KMAX)=U(I,KLES)
       DO 400 I=2,ILES
         IM=I-1
         IP=I+1
          DTTI=DT*RDX
      DO 40 K=2,KLES
   40   Y5(K)=X(I,K)+X(I,K-1)
       DO 42 K=2,KLES
         KP=K+1
         KM=K-1
         Y1(K)=0.
         Y2(K)=0.
        IF (Y5(K) .GE. EPS) THEN
         Y1(K)=MAX(X(I,KM),X(I,K),X(I,KP),X1(I,KM),X1(I,K),X1(I,KP))
         Y2(K)=MIN(X(I,KM),X(I,K),X(I,KP),X1(I,KM),X1(I,K),X1(I,KP))
        ENDIF
   42  CONTINUE
       DO 44 K=2,KLES
         KP=K+1
         KM=K-1
         Y3(K)=0.
         Y4(K)=0.
        IF (Y5(K) .GE. EPS) THEN
          DTTK=DT*AM(K)*RDZ
         Y3(K)=(Y1(K)-X(I,K))/(DTTI*(MAX(UMOD(I,K),0.)*X(IM,K)
     1                              -MIN(UMOD(IP,K),0.)*X(IP,K))
     2                        +DTTK*(MAX(WMOD(I,K),0.)*X(I,KM)
     3                              -MIN(WMOD(I,KP),0.)*X(I,KP))+EPS)
         Y4(K)=(X(I,K)-Y2(K))/(DTTI*(MAX(UMOD(IP,K),0.)*X(I,K)
     1                              -MIN(UMOD(I,K),0.)*X(I,K))
     2                        +DTTK*(MAX(WMOD(I,KP),0.)*X(I,K)
     3                              -MIN(WMOD(I,K),0.)*X(I,K))+EPS)
        ENDIF
   44  CONTINUE
c       Y3(2)=Y3(3)
c       Y4(2)=Y4(3)
       DO 46 K=3,KLES
         KM=K-1
        IF (Y5(K) .GE. EPS) THEN
           W(I,K)=MIN(1.,Y4(KM),Y3(K))*MAX(0.,WMOD(I,K))
     1           +MIN(1.,Y4(K),Y3(KM))*MIN(0.,WMOD(I,K))
        ENDIF
   46  CONTINUE
  400 CONTINUE
      DO 48 K=3,KLES
       W(IMAX,K)=W(2,K)
       W(1,K)=W(ILES,K)
   48 CONTINUE

      DO 800 K=1,KMAX
      DO 800 I=1,IMAX
       UMOD(I,K)=U(I,K)
  800  WMOD(I,K)=W(I,K)

      RETURN
      END
