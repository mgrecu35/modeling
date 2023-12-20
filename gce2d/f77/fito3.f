  
********************* CLIRAD IR1  Date: October, 1994 ****************
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fito3 (npp1,p,at,ao)

      implicit none
c-----------------------------------------------------------------------
c--- Inputs: npp1 --- number of points in the vertical               ---
c---         p    --- pressure at the model levels                   ---
c-----------------------------------------------------------------------
c--- Outputs: at  --- tempreture at the model levels                 ---
c---          ao  --- ozone at the model levels                      ---
c-----------------------------------------------------------------------
      integer lay,n,npp1
      parameter (lay=88,n=61)

      real    p(lay),at(lay),ao(lay),pl(n),ta(n),wa(n),oa(n),pa(n)
      real    tl(n),wl(n),ol(n),w(n),tab(3),wk(n,4)
      integer itab(3),iop(2)
c----- ix=        1:trp; 2:mls; 3:mlw; 4:sas; 5:saw

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,ix,k,lun,int
      real    y

      save

      ix=1
      lun=90+ix

      OPEN( lun,FILE='trp.dat'
     1      ,FORM='FORMATTED',STATUS='OLD' )

      read (lun,51) (pl(i),tl(i),wl(i),ol(i),i=2,n)
   51 format (0p,2f10.3,1p,2e12.4)
      pl(1)=0.
      tl(1)=tl(2)
      wl(1)=wl(2)
      ol(1)=ol(2)
      do 25 i=1,n-1
        pa(i)=0.5*(pl(i)+pl(i+1))
        ta(i)=0.5*(tl(i)+tl(i+1))
        wa(i)=0.5*(wl(i)+wl(i+1))
        oa(i)=0.5*(ol(i)+ol(i+1))
   25 continue
      pa(n)=pl(n)
      ta(n)=tl(n)
      wa(n)=wl(n)
      oa(n)=ol(n)
      write (6,105)
      write (6,106) (i,pa(i),ta(i),oa(i),wa(i),i=1,n-1)

c-----------------------------------------------------------------------
c--- Starts interpolation                                            ---
c-----------------------------------------------------------------------
      iop(1)=4
      iop(2)=4
      int=1
      itab(1)=1
      itab(2)=0
      itab(3)=0
      call coeff (n,pa,ta,w,iop,int,wk)
      do 200 k=1,npp1-1
        y=p(k)
       call terp1 (n,pa,ta,w,y,int,tab,itab)
        at(k)=tab(1)
        print *,'k = ',k,'  y = ',y,'  at = ',at(k)
  200 continue
      call coeff (n,pa,oa,w,iop,int,wk)
      do 300 k=1,npp1-1
        y=p(k)
       call terp1 (n,pa,oa,w,y,int,tab,itab)
       ao(k)=tab(1)
       write(6,106) k,y,at(k),ao(k)
  300 continue
  105 format(5x,'k',7x,'pa',8x,'Ta',10x,'o3',8x,'qv')
  106 format(2x,i4,2f10.3,2e12.4)
      return
      end
