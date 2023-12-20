c
      SUBROUTINE RINIT_UV2
      PARAMETER (NX=514,NZ=43)
      PARAMETER (NZ94=4*9*NZ,NZ24=2*4*NZ)
      real            V1(NX,NZ),V2(NX,NZ),V3(NX,NZ),
     1                V4(NX,NZ),V5(NX,NZ),V6(NX,NZ),
     2                V7(NX,NZ),V8(NX,NZ)
      COMMON/DEBUG_V/ V1,V2,V3,V4,V5,V6,V7,V8

      real            U1(NX,NZ),U2(NX,NZ),U3(NX,NZ),
     1                U4(NX,NZ),U5(NX,NZ),U6(NX,NZ),
     2                U7(NX,NZ),U8(NX,NZ)
      COMMON/DEBUG_U/ U1,U2,U3,U4,U5,U6,U7,U8

      real U9(NX,NZ),V9(NX,NZ)
      COMMON/DEBUG_UUVV/ U9,V9

      DO K=1,NZ
         DO I=1,NX
            V1(I,K)=0.
            V2(I,K)=0.
            V3(I,K)=0.
            V4(I,K)=0.
            V5(I,K)=0.
            V6(I,K)=0.
            V7(I,K)=0.
            V8(I,K)=0.
C
            U1(I,K)=0.
            U2(I,K)=0.
            U3(I,K)=0.
            U4(I,K)=0.
            U5(I,K)=0.
            U6(I,K)=0.
            U7(I,K)=0.
            U8(I,K)=0.

            U9(I,K)=0.
            V9(I,K)=0.
         ENDDO
      ENDDO
C
      RETURN
      END
