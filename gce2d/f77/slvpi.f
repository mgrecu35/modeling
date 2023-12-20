cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine slvpi (iflg)

      implicit none
c     ******   solve 2-d pressure equation  (fast fourier transform)

      integer nx,nz,nxm,nzm,nx12,iflg
      parameter (NX=514,NZ=43)
      parameter (nxm=nx-1,nzm=nz-1,nx12=12*nx)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5(11)
      common/bb/ dt,d2t,ril2,f5

      real    f(nx,nz),ee(nx,nz)
      common/bsat/ f
      common/bsat1/ ee

      real    aux(nx,nz),ff(nx,nz)
      common/badv/ aux
      common/badv1/ ff

      real    tb(nz),qb(nz),rho1(nz),rho(nz),cc(nz),aa(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,cc,aa,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx12)
      common/ba/ y1,y2,y3,y4,y5

      real    ar(nxm,nzm),ai(nxm,nzm)
      common/rfft/ ar,ai

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    beta(nx)
      integer i,ihf,ihf1,ihf2,imi,k,k1
      real    pix

      save

c     ******   initialize the para for solver
      if(iflg.eq.1) go to 100
      ihf=il2/2+1
      ihf1=ihf+1
      ihf2=ihf1+1
c     ******
      pix=8.*atan(1.)*ril2
       beta(1)=0.
      do 20 i=2,iles
   20  beta(i)=2.*dz2*(cos(pix*(i-2))-1.)*rdx2
       beta(imax)=0.
      call fft (0,il2)
      return
c     ****   forward transformation in x
  100 continue
      do 22 k=2,kles
      do 22 i=2,iles
       ar(i-1,k-1)=f(i,k)*ril2
       ai(i-1,k-1)=0.
   22 continue
      call fft (1,il2)
      do 200 k=2,kles
       do 24 i=2,ihf1
   24   f(i,k)=ar(i-1,k-1)
       do 26 i=ihf2,iles
   26   f(i,k)=ai(i-1,k-1)
  200 continue
c     ******      ******
      do 30 i=2,iles
       aux(i,2)=1./(aa(2)-beta(i))
       ee(i,2)=aa(2)*aux(i,2)
       ff(i,2)=-f(i,2)*aux(i,2)
   30 continue
      do 300 k=3,kl2
      do 300 i=2,iles
       aux(i,k)=1./(aa(k)+cc(k)-beta(i)-cc(k)*ee(i,k-1))
       ee(i,k)=aa(k)*aux(i,k)
  300  ff(i,k)=(-f(i,k)+cc(k)*ff(i,k-1))*aux(i,k)
      y1(2)=0.
      do 32 i=3,iles
       y1(i)=(-f(i,kles)+cc(kles)*ff(i,kl2))
     1       /((1.-ee(i,kl2))*cc(kles)-beta(i))
   32 continue
      do 34 i=2,iles
   34  f(i,kles)=y1(i)
      do 350 k=3,kles
       k1=kmax+1-k
      do 350 i=2,iles
       y1(i)=ee(i,k1)*y1(i)+ff(i,k1)
  350  f(i,k1)=y1(i)
c     ****   backward transformation in x    ******
      do 400 k=2,kles
       do 40 i=2,ihf1
   40   ar(i-1,k-1)=f(i,k)
       do 42 i=ihf2,iles
   42   ai(i-1,k-1)=-f(i,k)
  400 continue
      do 44 k=1,kl2
       ai(1,k)=0.
       ai(ihf,k)=0.
   44 continue
      do 440 i=ihf1,il2
       imi=imax-i
      do 440 k=1,kl2
       ar(i,k)=ar(imi,k)
  440  ai(imi,k)=-ai(i,k)
      call fft (1,il2)
      do 500 k=2,kles
      do 500 i=2,iles
  500  f(i,k)=ar(i-1,k-1)
      return
      end
