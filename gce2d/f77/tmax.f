Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tmax (x,aa,aminut,big,small,id)
      implicit none
c   ****   find maximum values' routine
      integer nx,nz
      parameter (NX=514,NZ=43)
      character*4 aa
      integer id
      real    aminut,big,small
      real    x(nx,nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    cons(4)
      data cons/.01,1.,1000.,.0000001/

      integer i,k,ilo1,ilo2,klo1,klo2
      real    xa


      save

c
      big=0.
      small=0.
c     ****   search for the largest or/and smallest values
      do 10 k=2,kles
      do 10 i=2,iles
       xa=x(i,k)
       big=max (big,xa)
       small=min (small,xa)
   10 continue
c     ****   search for the location
      do 20 k=2,kles
      do 20 i=2,iles
       if(x(i,k).ge.big) then
        ilo1=i-1
        klo1=k-1
       endif
       if(x(i,k).le.small) then
        ilo2=i-1
        klo2=k-1
       endif
   20 continue
c     ************************************
      big=cons(id)*big
      small=cons(id)*small
      write(6,123) aminut,aa,big,ilo1,klo1,small,ilo2,klo2
  123 format(3x,5hamin=,f6.0,3x,a4,4x,4hbig=,f7.3,3x,2hi=,i3,3x,2hk=,
     1   i3,6x,6hsmall=,f7.3,3x,2hi=,i3,2x,2hk=,i3)
      return
      end
