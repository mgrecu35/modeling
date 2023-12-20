Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LKB(RR,RT,IFLAG)
C
C       TO DETERMINE THE LOWER BOUNDARY VALUE RT OF THE 
C       LOGARITHMIC PROFILES OF TEMPERATURE (IFLAG=1) 
C       OR HUMIDITY (IFLAG=2) IN THE ATMOSPHERE FROM ROUGHNESS 
C       REYNOLD NUMBER RR BETWEEN 0 AND 1000.  OUT OF RANGE
C       RR INDICATED BY RT=-999. BASED ON LIU ET AL.(1979)
C       JAS 36 1722-1723
c     New scalar RR relation from Moana Wave data.
C
      DIMENSION A(9,2),B(9,2),RAN(9)
      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
     10.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
     10.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./
      I=1
      IF(RR.LE.0..OR.RR.GE.1000.) GOTO 90
   10 CONTINUE
      IF(RR.LE.RAN(I)) GOTO 20
      I=I+1
      GOTO 10
   20 RT=A(I,IFLAG)*RR**B(I,IFLAG)
      GOTO 99
   90 RT=-999.
   99 RETURN
      END
