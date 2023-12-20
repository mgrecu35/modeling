
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zang (month,iday,cosz,rlat,hrl)

      implicit none
      integer month,iday
      real    cosz,rlat,hrl

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer imd(12)
      data imd/31,28,31,30,31,30,31,31,30,31,30,31/

      integer idayi,monthi
      real    ang,anx,cosphi,costph,costtp,delta,eq,fac
      real    phi,sinphi,sintph,sinttp,tphi,ttphi,x,xlat

      save

c
c      noczon=0
c      fixmon=12.
      if (month .gt. 12) month=month-12 
      if (hrl .ge. 24.) then
       hrl=hrl-24.0
       iday=iday+1
         if (iday .gt. imd(month)) then
          iday=iday-imd(month)
          month=month+1
         end if
      end if
       monthi=month
        idayi=iday
c      if (noczon .eq. 1) then
c      monthi=fixmon     
c      end if
       fac=.01745
       x=float((monthi-1)*30+idayi+monthi/2)
       xlat=rlat*fac
       print*,'month day   hour     cosz'
       phi=2.*3.14159*x/365.
       tphi=phi*2.
       ttphi=phi*3.
       cosphi=cos(phi)
       sinphi=sin(phi)
       costph=cos(tphi)
       sintph=sin(tphi)
       costtp=cos(ttphi)
       sinttp=sin(ttphi)
       delta=.006918-0.399912*cosphi+0.070257*sinphi-0.006758*
     *        costph+0.000907*sintph-0.002697*costtp+.00148*sinttp
       anx=(hrl-12.)*15.*fac
       eq=.000075+.001868*cosphi-.032077*sinphi-.014615*costph
     *     -.040849*sintph
       ang=anx+eq
       cosz=sin(xlat)*sin(delta)+cos(xlat)*cos(delta)*cos(ang)
      write (6,200) monthi,idayi,hrl,cosz
  200 format(1x,2i4,2f10.4)
      return
      end
