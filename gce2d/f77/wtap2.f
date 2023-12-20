cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap2(x,x4,n1,n2,itape)
c     ******   write 2-d array data to data file   ******
      real x(n1,n2)
      real*4 x4(n1,n2)
      save
      do k=1,n2
        do i=1,n1
          x4(i,k)=x(i,k)
        enddo
      enddo
      write (itape) x4
      return
      end
