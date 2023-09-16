!  SFM  04/06/2013  Code module added for WOlson
!
      SUBROUTINE RADTRAN(UMU, NLYR, tbout, BTEMP, LYRTEMP, LYRHGT, KEXT,
     $                  SALB, ASYM, FISOT, EMIS, EBAR)
C
C     CHRIS KUMMEROW
C     INCLUDES ASYMPTOTIC EXPRESSIONS FOR TERM3, TERM4, AND TERM5 IN
C     THE LIMIT OF SMALL EFFECTIVE OPTICAL DEPTHS; BILL OLSON FEB, 1995.
C
!     WSO 04/07/2013  large layer optical depth numerics (double precision) on outgoing radiance
!     WSO 04/07/2013  delta-eddington option
!     WSO 04/07/2013  there are many code changes relative to radtran.f, here

!f2py real, intent(out) :: tbout
      implicit none

      INTEGER *4 MAXLYR
      PARAMETER    ( MAXLYR = 250 )

      CHARACTER*1  POLN
      LOGICAL LAMBERT, PRNT(3)
      integer ier, imode, icode
      INTEGER NLYR, I, J, NANG, NN
      REAL  UMU, BTEMP, LYRTEMP(0:NLYR), LYRHGT(0:NLYR), 
     $      KEXT(NLYR), SALB(NLYR), ASYM(NLYR)
      REAL  FISOT, EMIS, EBAR

      logical delta
      real tbout
      real *8 fdelta
      real *8 kext_delta(maxlyr)
      real *8 salb_delta(maxlyr)
      real *8 asym_delta(maxlyr)

      real *8 l(maxlyr)
      real *8 h(maxlyr)
      real *8 b0(maxlyr)
      real *8 b1(maxlyr)
      real *8 w(2*maxlyr, 2*maxlyr)
      real *8 btempd
      real *8 bb(2*maxlyr)
      real *8 dp(maxlyr)
      real *8 dm(maxlyr)
      real *8 z(0:maxlyr)
      real *8 mu
      real *8 nu
      real *8 x(2*maxlyr)
      real *8 fisotd
      real *8 emisd
      real *8 ebard

      real *8 xa
      real *8 xb
      real *8 xc
      real *8 xd
      real *8 ya
      real *8 yb 
      real *8 yc 
      real *8 dz 
      real *8 xnu

      real *4 XIUP 
      real *4 I_IN(maxlyr+1, 100)
      real *4 IOUT(0:maxlyr)
      real *4 TERM1
      real *4 TERM2
      real *4 TERM3
      real *4 TERM4
      real *4 TERM5
      real *4 TERM4prime
      real *4 TERM5prime
      real *4 tb

!delta-Eddington option
      delta = .false.
      delta = .true.
      lambert = .false.
      if(.not. delta) then
        do i = 1, nlyr
          kext_delta(i) = dble(kext(i))
          salb_delta(i) = dble(salb(i))
          asym_delta(i) = dble(asym(i))
        end do
      else
        do i = 1, nlyr
          fdelta = dble(asym(i) * asym(i))
          kext_delta(i) = dble((1. - salb(i) * fdelta) * kext(i))
          salb_delta(i) = dble((1. - fdelta) * salb(i) /
     $     (1. - salb(i) * fdelta))
          asym_delta(i) = dble(asym(i) / (1. + asym(i)))
        end do
      end if

      fisotd = dble(fisot)
      emisd = dble(emis)
      ebard = dble(ebar)
      btempd = dble(btemp)

C
c      write(6,'(1x,"umu:  ",f7.2,"  nlyr: ",i5,"  tb: ",f7.2,
c     + "  btemp: ",f7.2,"  fisot: ",f7.2,
c     + /1x,"emis: ",f10.3,"  ebar: ",f10.3,"  lamb: ",l5, 
c     + /1x,"  lyrtemp(0): ",f7.2,"  lyrhgt(0): ",f7.2)')
c     + umu,nlyr,tb,btemp,fisot,emis,ebar,lambert,lyrtemp(0),
c     + lyrhgt(0)
c      do j=nlyr,1,-1
c      write(6,'(1x,"lyrtemp: ",f7.2,"  lyrhgt: ",f7.2,
c     + "  kext: ",f10.4,"  salb: ",f10.4,"  asym: ",f10.4)')
c     + lyrtemp(j),lyrhgt(j),kext(j),salb(j),asym(j)
c      end do
C                                CALCULATE SOME CONSTANTS
      Z(0) = dble(LYRHGT(0))
      DO 20  J = 1,NLYR
        Z(J)  = dble(LYRHGT(J))
        B0(J) = dble(LYRTEMP(J-1))
        B1(J) = dble((LYRTEMP(J) - LYRTEMP(J-1))/
     $   (LYRHGT(J) - LYRHGT(J-1)))
        L(J) = dsqrt(3.*KEXT_delta(J)*KEXT_delta(J)*
     $         (1. - SALB_delta(J))*
     $         (1. - SALB_delta(J)*ASYM_delta(J)))
        H(J) = dble(1.5*KEXT_delta(J)*
     $   (1. - SALB_delta(J)*ASYM_delta(J)))
  20  CONTINUE


C     FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
C     AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
C     TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE 
C     MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
C     ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

C     SET ALL MATRIX ELEMENTS TO ZERO	
      DO 45  I = 1,2*NLYR
       DO 45  J = 1,2*NLYR
        W(I,J) = 0.d0
  45  CONTINUE	

C     FILL IN THE NON-ZERO MATRIX ELEMENTS
      W(1,1)   = ((ebard - 2.d0)*L(1)/H(1)) + ebard
      W(1,2)   = ((2.d0 - ebard)*L(1)/H(1)) + ebard
      DO 50  I = 2,2*(NLYR-1),2
       W(I,I-1)   =  
     $  (1.d0 - L(I/2)/H(I/2))*dexp(+L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I,I  )   =  (1.d0 + L(I/2)/H(I/2))*
     $  dexp(-L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I,I+1)   = -(1.d0 - L(I/2+1)/H(I/2+1))
       W(I,I+2)   = -(1.d0 + L(I/2+1)/H(I/2+1))

       W(I+1,I-1) =  (1.d0 + L(I/2)/H(I/2))*
     $  dexp(+L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I+1,I)   =  (1.d0 - L(I/2)/H(I/2))*
     $  dexp(-L(I/2)*(Z(I/2)-Z(I/2-1)))
       W(I+1,I+1) = -(1.d0 + L(I/2+1)/H(I/2+1))
       W(I+1,I+2) = -(1.d0 - L(I/2+1)/H(I/2+1))
c        write(*,*) (w(i,i+j),j=-2,2)
  50  CONTINUE
      W(2*NLYR,2*NLYR-1) =  (1.d0 + L(NLYR)/H(NLYR))*dexp(+L(NLYR)*
     $                       (Z(NLYR)-Z(NLYR-1)))
      W(2*NLYR,2*NLYR)   =  (1.d0 - L(NLYR)/H(NLYR))*dexp(-L(NLYR)*
     $                       (Z(NLYR)-Z(NLYR-1)))

C     FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
      BB(1)    = ebard*btempd - ebard*B0(1) - (ebard - 2.d0)*B1(1)/H(1)
      DO 55  I = 2,2*(NLYR-1),2
       BB(I)   =  + B1(I/2)/H(I/2) - B1(I/2+1)/H(I/2+1)
       BB(I+1) =  - B1(I/2)/H(I/2) + B1(I/2+1)/H(I/2+1)
  55  CONTINUE
      BB(2*NLYR)  =  fisotd - B0(NLYR) - B1(NLYR)*(Z(NLYR) 
     $               - Z(NLYR-1) + 1.d0/H(NLYR))
C
C
C     MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
c      write(*,*) bb
c      write(*,*) b1
c      stop
      imode=0
      call band(imode,2*nlyr,w,bb,icode)
!      call linsysa(w,bb,x,2*nlyr,2*maxlyr,ier)
      
!     CALL LINPAK(NLYR, W, BB, RCOND)
!      do i=1,2*nlyr
!         write(*,*) x(i), bb(i)
!      enddo
!      stop
!

      do i=1,2*nlyr
!         bb(i)=x(i)
      enddo
!      do i = 1, nlyr
!        write(*, '("i: ", i5, "  bbim1: ", e15.6, 
!     $   "  bbi: ", e15.6)') i, bb(2*i-1), bb(2*i)
!      enddo
c     write(*,*) bb
c     stop
C     
c     WRITE(*,*) ' '
c     IF ( PRNT(2) ) WRITE(4,*) ' '
      DO 60  I = 1,NLYR
       DP(I) = BB(2*I-1)
       DM(I) = BB(2*I)
c       WRITE(6,653) I, DP(I), I, DM(I)      
c       IF ( PRNT(2) )  WRITE(4,653) I, DP(I), I, DM(I)
  60  CONTINUE
 653  FORMAT(10X,'D+ (',I2,') = ',E11.4,6X,'D-(',I2,') = ',E11.4)

C     AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE

      MU = dble(UMU)
      NU = dble(-UMU)
C     FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
C     *********************************************************************
      
      IF ( LAMBERT ) THEN 
C      CALCULATE THE DOWNWELLING FLUX THROUGH THE ATMOSPHERE AT 81 ANGLES
       NANG = 81
       DO 997  NN = 1,NANG
        XNU = -(2.d0*NN - 1.d0)/(NANG*2.d0)
        I_IN(NLYR+1,NN) = fisotd
C       LOOP THROUGH THE REMAINING LAYERS
        DO 100  J = NLYR,1,-1
C        CALCULATE RADIANCE FROM TOP OF LAYER "J"
         XA = B0(J) - 1.5d0*SALB_delta(J)*ASYM_delta(J)*XNU*B1(J)/H(J)
         XB = B1(J)
         XC = SALB_delta(J)*DP(J)*
     $    (1.d0 - 1.5d0*ASYM_delta(J)*XNU*L(J)/H(J))
         XD = SALB_delta(J)*DM(J)*
     $    (1.d0 + 1.5d0*ASYM_delta(J)*XNU*L(J)/H(J))
         YA = KEXT_delta(J)/XNU
         YB = YA + L(J)
         YC = YA - L(J)
         DZ = Z(J) - Z(J-1)

         TERM1 = I_IN(J+1,NN)*dexp(YA*DZ)
         TERM2 = XA*(1. - dexp(YA*DZ))
c     new tests
         if(dabs(ya*dz) .lt. 1.d-5) then
           term3=-xb*ya*dz*dz
         else
           TERM3 = XB/YA*(dexp(YA*DZ)*(1.d0 - YA*DZ) - 1.d0)
         end if
         if(dabs(yb*dz) .lt. 1.d-5) then
           term4=-xc*ya*dz
         else
           TERM4 = XC*YA/YB*(1.d0 - dexp(YB*DZ))
         end if
         if(dabs(yc*dz) .lt. 1.d-5) then
           term5=-xd*ya*dz
         else
           TERM5 = XD*YA/YC*(1.d0 - dexp(YC*DZ))
         end if
         I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
!        write(6,'(1x,"lev: ",i5,"  ang: ",i5,"  term1: ",f10.4,
!     +            "  term2: ",f10.4,"  term3: ",f10.4,"  term4: ",f10.4,
!     +            "  term5: ",f10.4,"  tbdn: ",f7.2)')
!     +   j,nn,term1,term2,term3,term4,term5,i_in(j,nn)
  100   CONTINUE
  997  CONTINUE
C
C      CALCULATE THE TOTAL DOWNWELLING FLUX REACHING THE SURFACE
       XIUP = 0.
       DO 47  NN = 1,NANG
        XIUP = XIUP + I_IN(1,NN)*(1./NANG)*(2.*NN-1.)/(2.*NANG)
  47   CONTINUE
       XIUP = 2.*XIUP
c       write(6,'(1x,"xiup: ",f7.2)') xiup
C
      ELSE

C      CALCULATE THE DOWNWELLING FLUX THROUGH THE ATMOSPHERE AT ANGLE MU
       NN = 22                      ! THIS IS A DUMMY INDEX
       I_IN(NLYR+1,NN) = fisotd
C      LOOP THROUGH THE REMAINING LAYERS
       DO 110  J = NLYR,1,-1
C       CALCULATE RADIANCE FROM TOP OF LAYER "J"
        XA = B0(J) - 1.5d0*SALB_delta(J)*ASYM_delta(J)*NU*B1(J)/H(J)
        XB = B1(J)
        XC = SALB_delta(J)*DP(J)*
     $   (1.d0 - 1.5d0*ASYM_delta(J)*NU*L(J)/H(J))
        XD = SALB_delta(J)*DM(J)*
     $   (1.d0 + 1.5d0*ASYM_delta(J)*NU*L(J)/H(J))
        YA = KEXT_delta(J)/NU
        YB = YA + L(J)
        YC = YA - L(J)
        DZ = Z(J) - Z(J-1)

        TERM1 = I_IN(J+1,NN)*dexp(YA*DZ)
        TERM2 = XA*(1. - dexp(YA*DZ))
c     new tests
        if(dabs(ya*dz) .lt. 1.d-5) then
          term3=-xb*ya*dz*dz
        else
          TERM3 = XB/YA*(dexp(YA*DZ)*(1.d0 - YA*DZ) - 1.d0)
        end if
        if(dabs(yb*dz) .lt. 1.d-5) then
          term4=-xc*ya*dz
        else
          TERM4 = XC*YA/YB*(1. - dexp(YB*DZ))
        end if
        if(abs(yc*dz) .lt. 1.d-5) then
          term5=-xd*ya*dz
        else
          TERM5 = XD*YA/YC*(1.d0 - dexp(YC*DZ))
        end if
        I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
!        write(*, '("down layer: ", i5, "  i_in: ", e15.6,
!     $   "  ya: ", e15.6, "  dz: ", e15.6, "  exp(ya*dz): ", e15.6,
!     $   "  yc: ", e15.6, "  exp(yc*dz): ", e15.6, "  term1: ", e15.6, 
!     $   "  term2: ", e15.6, "  term3: ", e15.6, 
!     $   "  term4: ", e15.6, "  term5: ", e15.6)') 
!     $   J, I_IN(J+1,NN), YA, DZ, dexp(YA*DZ), YC, dexp(YC*DZ),
!     $   TERM1, TERM2, TERM3, TERM4, TERM5
  110  CONTINUE
       XIUP = I_IN(1,22)
C
      ENDIF

C    
C
      IOUT(0) = emisd*btempd + (1.d0 - emisd)*XIUP
      DO 101  J = 1,NLYR
C      CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
       XA = B0(J) - 1.5d0*SALB_delta(J)*ASYM_delta(J)*MU*B1(J)/H(J)
       XB = B1(J)
       XC = SALB_delta(J)*DP(J)*
     $  (1.d0 - 1.5d0*ASYM_delta(J)*MU*L(J)/H(J))
       XD = SALB_delta(J)*DM(J)*
     $  (1.d0 + 1.5d0*ASYM_delta(J)*MU*L(J)/H(J))
       YA = KEXT_delta(J)/MU
       YB = YA + L(J)
       YC = YA - L(J)
       DZ = Z(J) - Z(J-1)
 
       TERM1 = IOUT(J-1)*dexp(-YA*DZ)
       TERM2 = XA*(1.d0 - dexp(-YA*DZ))
c     new tests
       if(dabs(ya*dz) .lt. 1.d-5) then
         term3=0.
       else
         TERM3 = XB/YA*(YA*DZ - 1.d0 + dexp(-YA*DZ))
       end if
       if(dabs(yb*dz) .lt. 1.d-5) then
         TERM4=xc*ya*dz*dexp(-ya*dz)
       else
         TERM4 = XC*YA/YB*(dexp((YB-YA)*DZ) - dexp(-YA*DZ))
       end if
       if(dabs(yc*dz) .lt. 1.d-5) then
         TERM5=xd*ya*dz*dexp(-ya*dz)
         TERM5prime=xd*ya*dz*dexp(-ya*dz)
       else
!modification for large optical depth in layer
         TERM5 = XD*YA/YC*dexp(-YA*DZ)*(dexp(YC*DZ) - 1.d0)
         TERM5prime = XD*YA/YC*(dexp(-L(J)*DZ) - dexp(-YA*DZ))
       end if
!       IOUT(J) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
!        write(*, '("up   layer: ", i5, "  i_ou: ", e13.6, 
!     $   "  ya: ", e13.6, "  exp: ", e13.6, "  L: ", e13.6,
!     $   "  yb: ", e13.6, "  exp: ", e13.6, "  term1: ", e13.6, 
!     $   "  term2: ", e13.6, "  term3: ", e13.6, 
!     $   "  term4: ", e13.6, "  term5: ", e13.6, "  term5prime: ", e13.6)')
!     $   J, IOUT(J-1), YA, dexp(-YA*DZ), L(J), YB, dexp(-L(J)*DZ),
!     $   TERM1, TERM2, TERM3, TERM4, TERM5, TERM5prime
       IOUT(J) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5prime
  101 CONTINUE
C
C
c      WRITE(*,*) ' '
c      IF ( PRNT(2) )  WRITE(4,*) ' '
      DO 44  J = 0,NLYR
c       WRITE(*,77) NLYR+1-J, I_IN(NLYR+1-J,22), NLYR-J, IOUT(NLYR-J)
c       IF ( PRNT(2) ) WRITE(4,77) NLYR+1-J, I_IN(NLYR+1-J,22), 
c     $                NLYR-J, IOUT(NLYR-J)
  44  CONTINUE
  77  FORMAT(10X,'Iin(',I2,') = ',F7.1,10X,'Iout(',I2,') = ',F7.1)

C
      TB = IOUT(NLYR)
      tbout=tb
!      write(*,*) tb
      RETURN
      END
