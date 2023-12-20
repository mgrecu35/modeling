Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wstart (itape)
Ccccccccccccccc write to/from restart file cccccccccccccccccccccccccccccc
      parameter (NX=514,NZ=43,nnt=481)
      parameter (nz17=17*nz,nz4=4*nz,nz2=2*nz,nnt3=3*nnt,nx2=2*nx,
     1           nz9=9*nz,nnt7=7*nnt)
      common/bstart/ iw(3)
      common/rbstart/ rw(3)
c 8/21/01 shie found the dummy "sltm" was supposed to be as follows:
c  based on "gcem3d.F".  Put here for reference.
c COMMON/SLTM/ RLAT,RMONTH,RIDAY,HRL,SO0,COSZ,ICOSZ,TERMAN,RSFC
      common/sltm/ rlat(9) ! 8/21/01 shie, dummy!
      common/sun/ isun(2)  ! 8/21/01 shie, dummy!
      common/mbudget/ acoc(nz17)
      common/tqave/ tavet(nz2)
      common/tbudget/ avett(nnt3)
      common/tbudget1/ t_sfcq(nnt7)
      common/bsfc/ tsfc_1(nx2)
      common/gbs/ tlsw(nz9)
      common/gbs11/ tlsw1(nz4)
      common/timestat/ ndt_stat(3)
      COMMON/SFLUXS/ SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      save

      if(itape.gt.6) go to 100
        write(itape) iw
        write(itape) isun     ! 8/21/01 shie, dummy!
        write(itape) ndt_stat
        write(itape) rw
        write(itape) rlat     ! 8/21/01 shie, dummy!
c        write(itape) acoc
c        write(itape) tavet
        write(itape) avett
        write(itape) t_sfcq
c        write(itape) tsfc_1
c        write(itape) tlsw
c        write(itape) tlsw1
        write(itape) suw
        write(itape) svw
        write(itape) swt
        write(itape) swq

      return
  100   read(itape) iw
        read(itape) isun       ! 8/21/01 shie, dummy!
        read(itape) ndt_stat
        read(itape) rw
        read(itape) rlat       ! 8/21/01 shie, dummy!
c        read(itape) acoc
c        read(itape) tavet
        read(itape) avett
        read(itape) t_sfcq
c        read(itape) tsfc_1
c        read(itape) tlsw
c        read(itape) tlsw1
        read(itape) suw
        read(itape) svw
        read(itape) swt
        read(itape) swq
      return
      end
