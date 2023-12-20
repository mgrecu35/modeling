cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap3(x,x4,n1,n2,n3,itape)
c     ******   write 3-d array data to data file   ******
      real x(n1,n2,n3)
      real*4 x4(n1,n2,n3)
      save
      do k=1,n3
        do j=1,n2
          do i=1,n1
            x4(i,j,k)=x(i,j,k)
          enddo
        enddo
      enddo
      write (itape) x4
      return
      end
