C-----------------------------------------------------------------------------
      subroutine wtap (x,itape)
c     ******   write data to/from restart file   ******
      parameter (NX=514,NZ=43)
      real x(nx,nz)
      save
      if (itape.gt.6) go to 100
      write (itape) x
      return
  100  read (itape) x
      return
      end
