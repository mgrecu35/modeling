
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wstrad2 (itape)  ! newly added! 8/21/01 shie
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   8/21/01 by shie, write more radiation infromation to restart file

      parameter (NX=514,lay=88)

      common/srflx/ sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)

      common/surface/ sun_4(nx,4)

      save
      if(itape.gt.6) go to 100
       write(itape) sfir
       write(itape) sfsw
       write(itape) shir
       write(itape) shsw
       write(itape) salpha
       write(itape) si
       write(itape) sun_4
      return
  100  read(itape) sfir
       read(itape) sfsw
       read(itape) shir
       read(itape) shsw
       read(itape) salpha
       read(itape) si
       read(itape) sun_4
C
      RETURN
      END
