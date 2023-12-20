
C --------------------------------------------------------------------------

      FUNCTION RAN3 (idum)
      implicit none
      integer mbig,mseed,mz,idum
      real    ran3,fac
      PARAMETER  (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      integer i,k,iff,mj,mk,ii,inext,inextp,ma(55)
      SAVE  inext,inextp,ma,iff
      DATA  iff /0/
      IF( idum.LT.0 .OR. iff .EQ. 0 ) THEN
        iff= 1
        mj= mseed - ABS( idum )
c       mj=MOD( mj,mbig )
        mj= mj-mbig*(mj/mbig)
        ma(55)= mj
        mk= 1
        DO 10 i=1,54
c          ii= MOD( 21*i,55 )
          ii= (21*i)-55*((21*i)/55)
          ma(ii)= mk
          mk= mj-mk
          IF( mk.LT.mz ) mk= mk+mbig
          mj= ma(ii)
 10     CONTINUE
        DO 20 k=1,4
          DO 15 i=1,55
c            ma(i)= ma(i) - ma(1+MOD( i+30,55 ))
            ma(i)= ma(i) - ma(1+(i+30)-55*((i+30)/55))
            IF( ma(i).LT.mz ) ma(i)= ma(i)+mbig
 15       CONTINUE
 20     CONTINUE
        inext= 0
        inextp= 31
        idum= 1
      ENDIF
      inext= inext+1 
        IF( inext .EQ. 56 ) inext= 1
      inextp= inextp+1
        IF( inextp .EQ. 56 ) inextp= 1
      mj= ma(inext) - ma(inextp)
        IF( mj .LT. mz ) mj= mj+mbig
      ma(inext)= mj
      ran3= REAL( mj )*fac
      RETURN

      END
