cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundy (x,x1)
      parameter (NX=514,NZ=43)
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    x(nx,nz),x1(nx,nz)
c     ******   periodic lateral boundary condition
      do 10 k=1,kmax
       x(1,k)=x(iles,k)
       x1(1,k)=x1(iles,k)
       x(imax,k)=x(2,k)
       x1(imax,k)=x1(2,k)
   10 continue
cc    ****   set boundary condition in z-direction
      do 20 i=1,imax
       x(i,1)=x(i,2)
       x1(i,1)=x1(i,2)
       x(i,kmax)=x(i,kles)
   20  x1(i,kmax)=x1(i,kles)
      return
      end
