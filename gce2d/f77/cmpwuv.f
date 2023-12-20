


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmpwuv (id)

      implicit none
c     *******   compute new w, u and v   *****************************

      integer nx,nz,nx10,itt,nt,id,mt
c     parameter (NX=514,NZ=43,ITT=244)
      parameter (NX=514,NZ=43,NT=2880,ITT=244) !6/21/01
      parameter (nx10=10*nx)

      integer LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    uu1(nx,nz),vv1(nx,nz),ww1(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz)
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    qcl(nx,nz),qrn(nx,nz),qci(nx,nz),qcs(nx,nz),qcg(nx,nz)
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg

      real    pi(nx,nz)
      common/bsat/ pi

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real ppress(nx,nz)
      COMMON/SUE/ PPRESS

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
c
      integer     IT(NX),IV(NT),ICS(NX,4),IBZ(NX,4)
      COMMON/BCH/ IT,IV,ICS,IBZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ICS(I,3)=1:  STRATIFORM REGION                                  C
C      ICS(I,2)=1:  CONVECTIVE REGION                                  C
C      ICS(I,4)=1:  NO SFC RAIN BUT CLOUD ALOFT                        C
C      ICE(I,1)=1:  TOTAL MODEL DOMAIN                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

      real             SDIFF_VX(NZ,4),SDIFF_VZ(NZ,4),SGRID_VX(NZ,4),
     1                 SGRID_VZ(NZ,4),SDIFF_NV(NZ,4),SV_LARGE(NZ,4),
     2                 SPRE_V(NZ,4),SDT_V(NZ,4),TN_VWIND(NZ,4)
      COMMON/DEBUG_SV/ SDIFF_VX,SDIFF_VZ,SGRID_VX,
     1                 SGRID_VZ,SDIFF_NV,SV_LARGE,
     2                 SPRE_V,SDT_V,TN_VWIND

      real             SDIFF_UX(NZ,4),SDIFF_UZ(NZ,4),SGRID_UX(NZ,4),
     1                 SGRID_UZ(NZ,4),SDIFF_NU(NZ,4),SU_LARGE(NZ,4),
     2                 SPRE_U(NZ,4),SDT_U(NZ,4),TN_UWIND(NZ,4)
      COMMON/DEBUG_SU/ SDIFF_UX,SDIFF_UZ,SGRID_UX,
     1                 SGRID_UZ,SDIFF_NU,SU_LARGE,
     2                 SPRE_U,SDT_U,TN_UWIND

      real        UB_MEAN(NZ,4),VB_MEAN(NZ,4)

      COMMON/SUV/ UB_MEAN,VB_MEAN

      real               UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

      real             SSDIFF_VX(NZ,4),SSDIFF_VZ(NZ,4),SSGRID_VX(NZ,4),
     1                 SSGRID_VZ(NZ,4),SSDIFF_NV(NZ,4),SSV_LARGE(NZ,4),
     2                 SSPRE_V(NZ,4),SSDT_V(NZ,4),STN_VWIND(NZ,4)
      COMMON/DEBUG_SVS/ SSDIFF_VX,SSDIFF_VZ,SSGRID_VX,
     1                  SSGRID_VZ,SSDIFF_NV,SSV_LARGE,
     2                  SSPRE_V,SSDT_V,STN_VWIND

      real             SSDIFF_UX(NZ,4),SSDIFF_UZ(NZ,4),SSGRID_UX(NZ,4),
     1                 SSGRID_UZ(NZ,4),SSDIFF_NU(NZ,4),SSU_LARGE(NZ,4),
     2                 SSPRE_U(NZ,4),SSDT_U(NZ,4),STN_UWIND(NZ,4)
      COMMON/DEBUG_SUS/ SSDIFF_UX,SSDIFF_UZ,SSGRID_UX,
     1                  SSGRID_UZ,SSDIFF_NU,SSU_LARGE,
     2                  SSPRE_U,SSDT_U,STN_UWIND

      real         SUB_MEAN(NZ,4),SVB_MEAN(NZ,4)

      COMMON/SSUV/ SUB_MEAN,SVB_MEAN

      real qc_t,qc_tl,ww_t

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a1,a2,a11,a22,a33,a44,a55,a66,a77,a88,a99,ab11
      real    b11,b22,b33,b44,b55,b66,b77,b88,b99,ab22
      save
c     ******************************************************************
        a1=.5*cp*d2t*rdz
        a2=cp*d2t*rdx
        IF (LIPPS .EQ. 1) A1=CP*D2T*RDZ
       do 10 k=2,kles
        y1(k)=am1(k)*(tb(k-1)+tb(k))*a1
        y2(k)=tb(k)*a2
        IF (LIPPS .EQ. 1) then
          Y1(K)=AM1(K)*A1
          y2(k)=a2
        endif
   10   pi(1,k)=pi(iles,k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 1000 k=2,kles
       if(k.eq.2) go to 150
       do 100 i=2,iles
        w(i,k)=w(i,k)-y1(k)*(pi(i,k)-pi(i,k-1))
  100   ww1(i,k)=ww1(i,k)+eps*w(i,k)
  150  continue
       do 175 i=2,iles
          PPRESS(I,K)=PI(I,K)
c
c         PRE_U(I,K)=PRE_U(I,K)-Y2(K)*(PI(I,K)-PI(I-1,K))/d2t
         PRE_U(I,K)=-Y2(K)*(PI(I,K)-PI(I-1,K))/d2t
c
        u(i,k)=u(i,k)-y2(k)*(pi(i,k)-pi(i-1,k))
  175   uu1(i,k)=uu1(i,k)+eps*u(i,k)
 1000 continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      if (iuvbar .eq. 1) then
c        do k=2,kles
c        do i=2,iles
c          u(i,k)=u(i,k)+ubt(k)*dt
c          uu1(i,k)=uu1(i,k)+ubt(k)*dt
c          v(i,k)=v(i,k)+vbt(k)*dt
c          vv1(i,k)=vv1(i,k)+vbt(k)*dt
c        enddo
c        enddo
c      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 20 i=2,iles
       w(i,2)=0.
       ww1(i,2)=0.
   20 continue
      call boundy (u,uu1)
      call boundy (v,vv1)
      call boundy (w,ww1)
      do 30 i=1,imax
       w(i,2)=0.
       ww1(i,2)=0.
       w(i,kmax)=0.
       ww1(i,kmax)=0.
   30 continue
cc    ***********************
      do 40 k=1,kmax
      y1(k)=0.
      y2(k)=0.
      y3(k)=0.
   40 y4(k)=0.
      do 50 i=2,iles
      do 50 k=1,kmax
       y1(k)=y1(k)+u(i,k)
       y2(k)=y2(k)+uu1(i,k)
       y3(k)=y3(k)+v(i,k)
   50  y4(k)=y4(k)+vv1(i,k)
      do 60 k=1,kmax
       ub(k)=y1(k)*ril2
       ub1(k)=y2(k)*ril2
       vb(k)=y3(k)*ril2
       vb1(k)=y4(k)*ril2
   60 continue
c
      DO K=2,KLES
         DO I=2,ILES
c           DT_UWIND(I,K)=DT_UWIND(I,K)+(U(I,K)-UU1(I,K))/DT
c           DT_VWIND(I,K)=DT_VWIND(I,K)+(V(I,K)-VV1(I,K))/DT
           DT_UWIND(I,K)=(UU1(I,K)-UTEM_DT(I,K))/D2T
           DT_VWIND(I,K)=(VV1(I,K)-VTEM_DT(I,K))/D2T
         ENDDO
      ENDDO
C
CC
C
      IF (ID .EQ. 1) THEN

          do i=2,iles
             ICS(I,1)=1
          enddo

         DO MT=1,4

            DO K=2,KLES
            A11=0.
            A22=0.
            A33=0.
            A44=0.
            A55=0.
            A66=0.
            A77=0.
            A88=0.
            A99=0.
c
            B11=0.
            B22=0.
            B33=0.
            B44=0.
            B55=0.
            B66=0.
            B77=0.
            B88=0.
            B99=0.
c
            AB11=0.
            AB22=0.
               DO I=2,ILES
                  IF (ICS(I,MT) .EQ. 1) THEN
                     A11=A11+DIFF_VX(I,K)
                     A22=A22+DIFF_VZ(I,K)
                     A33=A33+GRID_VX(I,K)
                     A44=A44+GRID_VZ(I,K)
                     A55=A55+DIFF_NV(I,K)
                     A66=A66+V_LARGE(I,K)
                     A77=A77+PRE_V(I,K)
                     A88=A88+DT_VWIND(I,K)
                     A99=A99+1.
                     AB11=AB11+V(I,K)
CCCCC
                     B11=B11+DIFF_UX(I,K)
                     B22=B22+DIFF_UZ(I,K)
                     B33=B33+GRID_UX(I,K)
                     B44=B44+GRID_UZ(I,K)
                     B55=B55+DIFF_NU(I,K)
                     B66=B66+U_LARGE(I,K)
                     B77=B77+PRE_U(I,K)
                     B88=B88+DT_UWIND(I,K)
                     B99=B99+1.
                     AB22=AB22+U(I,K)
                  ENDIF
               ENDDO
C
               SDIFF_VX(K,MT)=SDIFF_VX(K,MT)+A11
               SDIFF_VZ(K,MT)=SDIFF_VZ(K,MT)+A22
               SGRID_VX(K,MT)=SGRID_VX(K,MT)+A33
               SGRID_VZ(K,MT)=SGRID_VZ(K,MT)+A44
               SDIFF_NV(K,MT)=SDIFF_NV(K,MT)+A55
               SV_LARGE(K,MT)=SV_LARGE(K,MT)+A66
               SPRE_V(K,MT)=SPRE_V(K,MT)+A77
               SDT_V(K,MT)=SDT_V(K,MT)+A88
               TN_VWIND(K,MT)=TN_VWIND(K,MT)+A99
c 7/27/01 tao, collect VB_MEAN, UB_MEAN every "id", 120 sec= 2 min
c                IF (A99 .EQ. 0.) A99=1.E20 ! tao 7/27/01
c              VB_MEAN(K,MT)=VB_MEAN(K,MT)+AB11/A99 ! tao 7/27/01, every "id",  120 sec= 2 min
c 8/4/01 tao
           IF (A99 .GE. 0.01) VB_MEAN(K,MT)=VB_MEAN(K,MT)+AB11/A99
CCCCC
               SDIFF_UX(K,MT)=SDIFF_UX(K,MT)+B11
               SDIFF_UZ(K,MT)=SDIFF_UZ(K,MT)+B22
               SGRID_UX(K,MT)=SGRID_UX(K,MT)+B33
               SGRID_UZ(K,MT)=SGRID_UZ(K,MT)+B44
               SDIFF_NU(K,MT)=SDIFF_NU(K,MT)+B55
               SU_LARGE(K,MT)=SU_LARGE(K,MT)+B66
               SPRE_U(K,MT)=SPRE_U(K,MT)+B77
               SDT_U(K,MT)=SDT_U(K,MT)+B88
               TN_UWIND(K,MT)=TN_UWIND(K,MT)+B99
           IF (B99 .GE. 0.01) UB_MEAN(K,MT)=UB_MEAN(K,MT)+AB22/B99
            ENDDO
         ENDDO
C
cc
c
         do mt=1,4
            do i=1,nx
c              ibx(i,mt)=0
               ibz(i,mt)=0 ! by shie 8/1/01
            enddo
         enddo
          DO K=2,KLES
            do i=2,iles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            qc_t=qcl(i,k)+qrn(i,k)+qci(i,k)+qcs(i,k)+qcg(i,k)
     1          +qcl(i-1,k)+qrn(i-1,k)+qci(i-1,k)+qcs(i-1,k)+qcg(i-1,k) ! shie
c    1            +qcl(i-1,k)+qrn(i-1,k)+qci(i-1,k)+qcs(i-1,k)+qcg(i-1,k)
c 8/6/01 tao, shie qc_tl (saturated) and qc_t (all) were found no significant
c         difference in impact, so comment out the sorting for qc_tl.
c           qc_tl=qcl(i,k)+qci(i,k)+qcl(i-1,k)+qci(i-1,k)
            ww_t=.0025*(w(i,k)+w(i,k-1)+w(i-1,k)+w(i-1,k-1))

              if (qc_t .ge. 2.e-5 .and. ww_t .ge. 1.) ibz(i,1)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .le. -0.5) ibz(i,2)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .ge. .001) ibz(i,3)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .le. -0.001) ibz(i,4)=1
c             if (qc_tl .ge. 2.e-5 .and. ww_t .ge. 1.) ibz(i,3)=1
c             if (qc_tl .ge. 2.e-5 .and. ww_t .le. -0.5) ibz(i,4)=1
            enddo
           ENDDO
c
              DO MT=1,4
c
         DO K=2,KLES
c
            A11=0.
            A22=0.
            A33=0.
            A44=0.
            A55=0.
            A66=0.
            A77=0.
            A88=0.
            A99=0.
c
            B11=0.
            B22=0.
            B33=0.
            B44=0.
            B55=0.
            B66=0.
            B77=0.
            B88=0.
            B99=0.
c
            AB11=0.
            AB22=0.
              do i=2,iles
                IF (IBZ(I,MT) .EQ. 1) THEN
                   A11=A11+DIFF_VX(I,K)
                   A22=A22+DIFF_VZ(I,K)
                   A33=A33+GRID_VX(I,K)
                   A44=A44+GRID_VZ(I,K)
                   A55=A55+DIFF_NV(I,K)
                   A66=A66+V_LARGE(I,K)
                   A77=A77+PRE_V(I,K)
                   A88=A88+DT_VWIND(I,K)
                   A99=A99+1.
                   AB11=AB11+V(I,K)
CCCCC
                   B11=B11+DIFF_UX(I,K)
                   B22=B22+DIFF_UZ(I,K)
                   B33=B33+GRID_UX(I,K)
                   B44=B44+GRID_UZ(I,K)
                   B55=B55+DIFF_NU(I,K)
                   B66=B66+U_LARGE(I,K)
                   B77=B77+PRE_U(I,K)
                   B88=B88+DT_UWIND(I,K)
                   B99=B99+1.
                   AB22=AB22+U(I,K)
                ENDIF
              ENDDO
C
               SSDIFF_VX(K,MT)=SSDIFF_VX(K,MT)+A11
               SSDIFF_VZ(K,MT)=SSDIFF_VZ(K,MT)+A22
               SSGRID_VX(K,MT)=SSGRID_VX(K,MT)+A33
               SSGRID_VZ(K,MT)=SSGRID_VZ(K,MT)+A44
               SSDIFF_NV(K,MT)=SSDIFF_NV(K,MT)+A55
               SSV_LARGE(K,MT)=SSV_LARGE(K,MT)+A66
               SSPRE_V(K,MT)=SSPRE_V(K,MT)+A77
               SSDT_V(K,MT)=SSDT_V(K,MT)+A88
               STN_VWIND(K,MT)=STN_VWIND(K,MT)+A99
c                IF (A99 .EQ. 0.) A99=1.E20
c              SVB_MEAN(K,MT)=SVB_MEAN(K,MT)+AB11/A99
         IF (A99 .GE. 0.01) SVB_MEAN(K,MT)=SVB_MEAN(K,MT)+AB11/A99 ! shie
CCCCC
               SSDIFF_UX(K,MT)=SSDIFF_UX(K,MT)+B11
               SSDIFF_UZ(K,MT)=SSDIFF_UZ(K,MT)+B22
               SSGRID_UX(K,MT)=SSGRID_UX(K,MT)+B33
               SSGRID_UZ(K,MT)=SSGRID_UZ(K,MT)+B44
               SSDIFF_NU(K,MT)=SSDIFF_NU(K,MT)+B55
               SSU_LARGE(K,MT)=SSU_LARGE(K,MT)+B66
               SSPRE_U(K,MT)=SSPRE_U(K,MT)+B77
               SSDT_U(K,MT)=SSDT_U(K,MT)+B88
               STN_UWIND(K,MT)=STN_UWIND(K,MT)+B99
c                IF (B99 .EQ. 0.) B99=1.E20
c              SUB_MEAN(K,MT)=SUB_MEAN(K,MT)+AB22/B99
         IF (B99 .GE. 0.01) SUB_MEAN(K,MT)=SUB_MEAN(K,MT)+AB22/B99 ! shie
           ENDDO
           ENDDO

      ENDIF

      return
      end
