c
      subroutine sat_zero
c     ***************************************
      parameter (NX=514,NZ=43)
      parameter (nx16=16*nx)
      common/bsat/ xx1(nx,nz)
      common/bsat1/ xx2(nx,nz)
      common/badv1/ xx3(nx,nz)
      common/b1c/ qc(nx,nz)
      common/b1r/ qr(nx,nz)
      common/b1i/ qi(nx,nz)
      common/b1s/ qs(nx,nz)
      common/b1g/ qg(nx,nz)
      common/ba/ y1(nx16)
       CMIN=1.e-20
       do 100 k=1,nz
       do 100 i=1,nx
          xx1(i,k)=0.
          xx2(i,k)=0.
          xx3(i,k)=0.
        if (qc(i,k) .le. CMIN) qc(i,k)=0.0
        if (qr(i,k) .le. CMIN) qr(i,k)=0.0
        if (qi(i,k) .le. CMIN) qi(i,k)=0.0
        if (qs(i,k) .le. CMIN) qs(i,k)=0.0
        if (qg(i,k) .le. CMIN) qg(i,k)=0.0
  100  continue
       do 10 i=1,nx16
          y1(i)=0.
   10  continue
      return
      end
