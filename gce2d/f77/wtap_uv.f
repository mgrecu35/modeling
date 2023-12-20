C
      SUBROUTINE WTAP_UV
      PARAMETER (NX=514,NZ=43)
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

C
      REAL BNXNZ(NX,NZ)
C
       CALL WTAP2 ( DIFF_VX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_VZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_VX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_VZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_NV,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( V_LARGE,BNXNZ,NX,NZ,4)
       CALL WTAP2 (   PRE_V,BNXNZ,NX,NZ,4)
       CALL WTAP2 (DT_VWIND,BNXNZ,NX,NZ,4)
CCCCCCC
       CALL WTAP2 ( DIFF_UX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_UZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_UX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_UZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_NU,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( U_LARGE,BNXNZ,NX,NZ,4)
       CALL WTAP2 (   PRE_U,BNXNZ,NX,NZ,4)
       CALL WTAP2 (DT_UWIND,BNXNZ,NX,NZ,4)
      RETURN
      END