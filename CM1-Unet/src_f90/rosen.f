*=======================================================
*  ROSENKRANZ (1998) model  -- reference "Water vapor microwave
*  continuum absorption: a comparison of measurements and results"
*  To appear in Radio Science
*
* Top level interface (GasabsR98) by G. Petty
*
*   This subroutine calculates the mass extinction coefficient of the
*   dry and vapor components of air.
*
* Input:
*     F  = frequency (GHz),
*     Tk = absolute temperature (K)
*     Rhowv = water vapor density (kg/m**3).
*     Pa = Total air pressure (Pascals).
* Output:
*     Absair = extinction by dry air  (meters squared per kg *moist air*)
*     Abswv  = extinction by water vapor (meters squared per kg water vapor)
*
*  2/8/2001 - fixed division by zero when rhowv = 0  (G. Petty)
*


      SUBROUTINE GasabsR98(F,Tk,Rhowv,Pa,absair,abswv,ireturn)
cf2py intent(out) absair,abswv
* check for "reasonable" input values
      ireturn=1
      if (F .le. 0.0 .or. F .gt. 800.0) then
         return
         stop 'Frequency out of range in GasabsR98'
      endif
      if (Tk .lt. 100.0) then
         print*, Tk
         return
         stop 'Temperature out of range in GasabsR98'

      endif
      if (pa .lt. 10.0 .or. pa .gt. 1.2e5) then
         return
         stop 'Pressure out of range in GasabsR98'
      endif
      ireturn=0
* convert pressure from Pa to Mb
      pmb = pa/100.0

* convert vapor density from kg/m**3 to g/m**3

      vapden = rhowv*1000.0

* get volume extinction coefficients
      absair = absn2(tk,pmb,f) + o2abs(tk,pmb,vapden,f)
      abswv = abh2o(tk,pmb,vapden,f)
* 
* convert vapor density to vapor pressure
      e = rhowv*(tk*461.5) 
* calculate specific humidity
      q = 0.622*e/pa
* calculate virtual temperature
      tv = (1. + 0.61*q)*tk
* moist air density
      rhoair = pa/(tv*287.06)

* convert above from Np/km to m**2/kg
      
      !absair = 0.001*absair/rhoair
!      if (rhowv .eq. 0.0) then
!         abswv = 0.0
!      else
!         abswv = 0.001*abswv/rhowv
!      endif

      return
      end

*******************************************************************
*This file contains subroutines for computing atmospheric absorption at
*microwave wavelengths, as supplied by P. Rosenkranz (2/98) via his
*anonymous ftp server (mesa.mit.edu; login anonymous; go to phil/lpl_rt)AN
*Consolidated into one file by G. Petty
*
***********************************************************
*  Begin summaries
***********************************************************
*     FUNCTION ABSN2(T,P,F)
C     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
C             (NEPER/KM)
C     T = TEMPERATURE (K)
C     P = PRESSURE (MB)
C     F = FREQUENCY (GHZ)
C
*************************************************************************
*     FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
C     PRES   MILLIBARS PRESSURE           (3 TO 1000)
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          (0 TO 900)
**************************************************************
*     FUNCTION ABH2O(T,P,RHO,F)
C
C PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
C 
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      RHO     G/M**3    I   WATER VAPOR DENSITY
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
C
*****************************************************************
*     FUNCTION ABLIQ(WATER,FREQ,TEMP)
C     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
C     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
C     (INT. J. IR & MM WAVES V.12(17) JULY 1991
C     WATER IN G/M**3
C     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
C     TEMP IN KELVIN
C        PWR 8/3/92
C
**************************************************************
*
*  Begin actual function definitions
*
***************************************************************
      FUNCTION ABSN2(T,P,F)
C     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
C             (NEPER/KM)
C     T = TEMPERATURE (K)
C     P = PRESSURE (MB)
C     F = FREQUENCY (GHZ)
C
      TH = 300./T
      ABSN2 = 6.4E-14*P*P*F*F*TH**3.55
      RETURN
      END
*
*
      FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C      5/1/95  P. Rosenkranz
C
C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQ
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
C     PRES   MILLIBARS PRESSURE           (3 TO 1000)
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          (0 TO 900)
C
C     REFERENCE FOR EQUATIONS AND COEFFICIENTS:
C     P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
C      BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
C     AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
C     (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
C
      COMMON /O2COM1/ X,WB300,W300(40),F(40),Y300(40),S300(40),
     & V(40),BE(40)
C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
     2  59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
     3  56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
     4  55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
     5  53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
     6  52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7631,
     7  487.2494, 715.3932, 773.8397, 834.1453/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,
     &  .3351E-14,.3292E-14, .3721E-14,.3891E-14,
     &  .3640E-14,.4005E-14, .3227E-14,.3715E-14,
     &  .2627E-14,.3156E-14, .1982E-14,.2477E-14,
     &  .1391E-14,.1808E-14, .9124E-15,.1230E-14,
     &  .5603E-15,.7842E-15, .3228E-15,.4689E-15,
     &  .1748E-15,.2632E-15, .8898E-16,.1389E-15,
     &  .4264E-16,.6899E-16, .1924E-16,.3229E-16,
     &  .8191E-17,.1423E-16, .6460E-15, .7047E-14, .3011E-14,
     &  .1826E-14, .1152E-13, .3971E-14/
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,
     & 2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,
     & 2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,
     & 2*7.744, .048, .044, .049, .145, .141, .145/
C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.63, 1.646, 1.468, 1.449, 1.382, 1.360,
     & 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,
     & 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,
     & 2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.92, 3*1.81/
      DATA Y300/  -0.0233,  0.2408, -0.3486,  0.5227,
     & -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,
     &  0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,
     &  0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,
     &  0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,
     &  0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,
     &  0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,
     &  0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,
     &  0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,
     &  0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,
     &  0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP/217.
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO 32 K=1,40
      DF = W300(K)*DEN
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))**2 + DF*DF)
      SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))**2 + DF*DF)
32    SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
      O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      RETURN
      END

**************************************************************
      FUNCTION ABH2O(T,P,RHO,F)
C
C  NAME- ABH2O    LANGUAGE- FORTRAN 77
C
C PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
C 
      IMPLICIT NONE
C  CALLING SEQUENCE PARAMETERS-
C    SPECIFICATIONS
      REAL T,P,RHO,F,ABH2O
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      RHO     G/M**3    I   WATER VAPOR DENSITY
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
C
C   REFERENCES-
C   LINE INTENSITIES FROM HITRAN92 (SELECTION THRESHOLD=
C     HALF OF CONTINUUM ABSORPTION AT 1000 MB).
C   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED:
C     H.J.LIEBE AND T.A.DILLON, J.CHEM.PHYS. V.50, PP.727-732 (1969) &
C     H.J.LIEBE ET AL., JQSRT V.9, PP. 31-47 (1969)  (22GHz);
C     A.BAUER ET AL., JQSRT V.37, PP.531-539 (1987) & 
C     ASA WORKSHOP (SEPT. 1989) (380GHz);
C     AND A.BAUER ET AL., JQSRT V.41, PP.49-54 (1989) (OTHER LINES).
C   AIR-BROADENED CONTINUUM BASED ON LIEBE & LAYTON, NTIA 
C     REPORT 87-224 (1987); SELF-BROADENED CONTINUUM BASED ON 
C     LIEBE ET AL, AGARD CONF. PROC. 542 (MAY 1993), 
C     BUT READJUSTED FOR LINE SHAPE OF
C     CLOUGH et al, ATMOS. RESEARCH V.23, PP.229-241 (1989).
C
C   REVISION HISTORY-
C    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
C          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
C                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
C          OCT. 24, 95  PWR -ADD 1 LINE.
C          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
C                       REVISED CONTINUUM.
C
C   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),
     & WS(NLINES),XS(NLINES)
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON
C     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
C     LINE INTENSITIES AT 300K:
      DATA S1/ .1310E-13, .2273E-11, .8036E-13, .2694E-11, .2438E-10,
     & .2179E-11, .4624E-12, .2562E-10, .8369E-12, .3263E-11, .6659E-12,
     & .1531E-08, .1707E-10, .1011E-08, .4227E-10/
C     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,
     & 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
C     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00281, .0023, .00278, .00287, .0021, .00186,
     & .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
C     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,
     & .71, .68, .70/
C     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,
     & .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
C     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,
     & 1.0, .68, .84, .78/
C
      IF(RHO.LE.0.) THEN
        ABH2O = 0.
        RETURN
      ENDIF
      PVAP = RHO*T/217.
      PDA = P -PVAP
      DEN = 3.335E16*RHO
      TI = 300./T
      TI2 = TI**2.5
C
C      CONTINUUM TERMS
      CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5)*PVAP*F*F 
C
C      ADD RESONANCES
      SUM = 0.
      DO 30 I=1,NLINES
      WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
      WSQ = WIDTH*WIDTH
      S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
      DF(1) = F - FL(I)
      DF(2) = F + FL(I)
C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      BASE = WIDTH/(562500. + WSQ)
C  DO FOR POSITIVE AND NEGATIVE RESONANCES
      RES = 0.
      DO 20 J=1,2
      IF(ABS(DF(J)).LT.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
20    CONTINUE
30    SUM = SUM + S*RES*(F/FL(I))**2
      ABH2O = .3183E-4*DEN*SUM + CON
      RETURN
      END
**************************************************************
      FUNCTION ABLIQ(WATER,FREQ,TEMP)
C     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
C     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
C     (INT. J. IR & MM WAVES V.12(17) JULY 1991
C     WATER IN G/M**3
C     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
C     TEMP IN KELVIN
C        PWR 8/3/92
C
      COMPLEX EPS,RE
      IF(WATER.LE.0.) THEN
       ABLIQ = 0.
       RETURN
      ENDIF
      THETA1 = 1.-300./TEMP
      EPS0 = 77.66 - 103.3*THETA1
      EPS1 = .0671*EPS0
      EPS2 = 3.52 + 7.52*THETA1
      FP = (316.*THETA1 + 146.4)*THETA1 +20.20
      FS = 39.8*FP
      EPS = (EPS0-EPS1)/CMPLX(1.,FREQ/FP) +
     & (EPS1-EPS2)/CMPLX(1.,FREQ/FS) +EPS2
      RE = (EPS-1.)/(EPS+2.)
      ABLIQ = -.06286*AIMAG(RE)*FREQ*WATER
      RETURN
      END



      FUNCTION VAPOR(T)
C
C  COMPUTES SATURATION H2O VAPOR PRESSURE (OVER LIQUID)
C  USING LIEBE'S APPROXIMATION (CORRECTED)
C
C  T IN KELVIN
C  VAPOR IN MBAR
C                   PWR 4/8/92
C
      TH = 300./T
      VAPOR = 35.3*EXP(22.64*(1.-TH))*TH**5
      RETURN
      END
