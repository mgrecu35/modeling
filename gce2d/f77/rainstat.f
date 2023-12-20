cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rainstat (rst)

      implicit none
      integer nx,nt,nstr,nxi
      PARAMETER (NX=514,NT=2880,NSTR=100)
      PARAMETER (NXI=NX-2)

      real    rst(nt,nxi)

      real    area(nstr),accum(nstr),xinter(nstr)
      integer ik,it,ix
      real    ri,rmin,sum1,sum2

      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ----------------------------------------------------------------------
      RMIN = 0.01
      DO 1 IK = 1,NSTR
        XINTER(IK) = FLOAT(IK)
        AREA(IK)=0.
        ACCUM(IK)=0.
    1 CONTINUE
      DO 100 IT = 1,NT
      DO 100 IX = 1,NXI
        RI=RST(IT,IX)
        IF (RI .GE. RMIN) THEN
          IF (RI .LT. XINTER(1)) THEN
            ACCUM(1) = ACCUM(1) + RI
            AREA(1) = AREA(1) + 1.
          ENDIF 
          DO 10 IK = 2,NSTR
            IF (RI .GE. XINTER(IK-1) .AND. RI .LT. XINTER(IK)) THEN
              ACCUM(IK) = ACCUM(IK) + RI
              AREA(IK) = AREA(IK) + 1.
            ENDIF
   10     CONTINUE
        ENDIF
  100 CONTINUE
      SUM1 = 0.
      SUM2 = 0.
      DO IK = 1,NSTR
       SUM1 = SUM1 + ACCUM(IK)
       SUM2 = SUM2 + AREA(IK)
      ENDDO
      DO IK = 1,NSTR
        IF (SUM1 .EQ. 0.) THEN
          ACCUM(IK) = 0.
        ELSE
          ACCUM(IK) = 100.*ACCUM(IK)/SUM1
        ENDIF
        IF (SUM2 .EQ. 0.) THEN
          AREA(IK) = 0.
        ELSE
          AREA(IK) = 100.*AREA(IK)/SUM2
        ENDIF
      ENDDO
      DO 200 IX = 1,NSTR
        WRITE (6,101) IX,XINTER(IX),ACCUM(IX)
        WRITE (6,102) IX,AREA(IX)
  200 CONTINUE
  101 FORMAT(12HACCUM(IX) = ,I8,4X,2F12.3)
  102 FORMAT(12HAREA(IX)  = ,I8,4X,2F12.3)
      RETURN
      END
