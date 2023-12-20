
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wtap0(x,itape)
c     ******   write 1-d array data to data file   ******
      real x
      real*4 x4
      save
      x4=x
      write (itape) x4
      return
      end
