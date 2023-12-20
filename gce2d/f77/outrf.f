cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outrf  (ircl,npp1,pl,pa)
      implicit none
      integer lay,iopt
      parameter(lay=88,iopt=8)

      integer ircl,npp1
      real    pl(lay),pa(lay)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)
      common/srflx/sfir,sfsw,shir,shsw,salpha,si

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    fir(4),fsw(4),hsw(4),hir(4)
      integer k,lc

      save


      write (iopt,605) ircl,(si(lc),lc=1,4)
      write (iopt,610)
      do 200  k=1,npp1
      do 100 lc=1,4
      if(si(lc).eq.0.0) go to 100
          fir(lc)=sfir(k,lc)/si(lc)
          fsw(lc)=sfsw(k,lc)/si(lc)
          hir(lc)=shir(k,lc)/si(lc)
          hsw(lc)=shsw(k,lc)/si(lc)
  100 continue
        write (iopt,600) k,pl(k),(fir(lc),lc=1,4),(fsw(lc),lc=1,4),
     +                     pa(k),(hir(lc),lc=1,4),(hsw(lc),lc=1,4)
  200 continue
      do 300 lc=1,4
      salpha(lc)=salpha(lc)/ircl
  300 continue
      write (iopt,615) salpha
  600 format(1x,i2,f8.2,8f8.2,f8.2,8f6.2)
  605 format(1x,'time count=',i4,10x,'cloud count=',4f10.0)
  610 format(2x,'k',6x,'pl',2x,'fi-tot',3x,'fi-cv',3x,'fi-cm',2x,
     +'fi-clr',2x,'fs-tot',3x,'fs-cv',3x,'fs-cm',2x,'fs-clr',
     + 4x,'pa',2x,'hi-tot',1x,'hi-cv',1x,'hi-cm',   'hi-clr',
     +            'hs-tot',1x,'hs-cv',1x,'hs-cm',   'hs-clr')
  615 format(1x,'albedo=',4f8.5)

      return
      end
