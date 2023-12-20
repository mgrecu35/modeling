Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tmaxadv (x,small)
C   ****   FIND MAXIMUM VALUES' routine
      PARAMETER (NX=514,NZ=43)
      real    X(NX,NZ)
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      save
      SMALL=0.
C     ****   SEARCH FOR THE SMALLEST VALUE
      DO 10 K=2,KLES
      DO 10 I=2,ILES
       XA=X(I,K)
       SMALL=MIN (SMALL,XA)
   10 CONTINUE
      RETURN
      END
