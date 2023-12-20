
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap1(x,x4,n,itape)
c     ******   write 1-d array data to data file   ******
      real x(n)
      real*4 x4(n)
      save
      do k=1,n
        x4(k)=x(k)
      enddo
      write (itape) x4
      return
      end

      subroutine isecwrite(x,x4,n,itape)


      end
