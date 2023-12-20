cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine wstrad (itape,cosz,month,iday,hrl)
      subroutine wstrad (itape,cosz,month,iday,hrl,rlat)  ! 8/21/01 shie
cccccccccccccccccc  write/read to restart file cccccccccccccccccccccc
      parameter (NX=514)
      real ri180(nx),riold180(nx)
      real cosz,hrl,rlat  ! 8/21/01 shie
      integer month,iday  ! 8/21/01 shie
      common/sfcri/ ri180,riold180,iricont
      save
      if(itape.gt.6) go to 100
       write(itape) month
       write(itape) iday
       write(itape) iricont
       write(itape) cosz
       write(itape) hrl
       write(itape) rlat  ! 8/21/01 shie
       write(itape) ri180
       write(itape) riold180
      return
  100  read(itape) month
       read(itape) iday
       read(itape) iricont
       read(itape) cosz
       read(itape) hrl
       read(itape) rlat  ! 8/21/01 shie
       read(itape) ri180
       read(itape) riold180
      return
      end
