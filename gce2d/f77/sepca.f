cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sepca (isec,ics,iv,aco,aco1,aan,aan1,lconv,lanvl,lnspt)
cc   ***  gce model's convective and anvil separation   *****
      parameter (NX=514,NZ=43)
      parameter (nx13=16*nx-2*nx-nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    qc(nx,nz),qr(nx,nz),qi(nx,nz),qs(nx,nz),qg(nx,nz),w(nx,nz)

      common/b2c/ qc
      common/b2r/ qr
      common/b2i/ qi
      common/b2s/ qs
      common/b2g/ qg
      common/b2w/ w

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real y1(nx),y3(nx),tair(nz),y2(nx13)
      common/ba/ y1,y3,tair,y2

      integer ics(nx,4),iv(nx),jv(nx),it(nx),itt(nx)
cc
       do 4 i=1,imax
          ics(i,1)=1
    4  continue
c
       do 1 ist=2,4
       do 1 i=1,imax
    1    ics(i,ist)=0
       do 2 i=1,imax
         itt(i)=4
         if(ri(i).ge.0.001) itt(i)=1
         it(i)=4
         iv(i)=0
         jv(i)=0
    2    if (ri(i).lt.0.001) it(i)=1

      if(isec.eq.isec/14400*14400) then
        write(6,4567) 
 3456 format(/3x,128i1/)
 4567  format(/5x,'gce models convective-stratiform partition')
       write(6,3456) (it(i),i=2,iles,4)
       write(6,3456) (itt(i),i=2,iles,4)
      endif

      do 3 k=2,kles
    3   tair(k)=pi(k)*tb(k)-273.16
c   ***  churchill and houze's method (sfc rainfall only)  ********
       ri(1)=ri(iles)
       ri(imax)=ri(2)
        do 214 i=2,iles
         im=i-1
         imm=i-2
          if (i .eq. 2) imm=il2
         ip=i+1
         ipp=i+2
          if (i .eq. iles) ipp=3
  214   y1(i)=.4*(ri(imm)+ri(im)+ri(i)+ri(ip)+ri(ipp))
       nite=0
       do 245 i=2,iles
        if (ri(i) .ge. .001) nite=nite+1
  245  continue
        nstep=0
  248  nstep=nstep+1
        big=0.0
        ibig=2
       do 242 i=2,iles
        a1=ri(i)
        if (it(i) .le. 3) go to 242
        if (a1 .lt. big) go to 242
        ibig=i
        big=a1
  242  continue
        ip1=ibig+1
        im1=ibig-1
        rave=y1(ibig)
         if (ri(ibig) .lt. rave) go to 246
        k=3
        kp1=3
        km1=3
        go to 244
  246    k=2
         kp1=4
         km1=4
        if (it(ip1) .eq. 3) kp1=3
        if (it(im1) .eq. 3) km1=3
  244   it(ibig)=k
        if (it(ip1) .gt. 2) it(ip1)=kp1
        if (it(im1) .gt. 2) it(im1)=km1
       if (nstep .lt. nite) go to 248

      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c  ***  find the location and value of max sfc rainfall grid point
        big=0.0
        ibig=2
       do 251 i=2,iles
        a1=ri(i)
        if (a1 .lt. big) go to 251
        ibig=i
        big=a1
  251  continue

c        rianmax=max(20., min(25., big))
c        if(isec .le. 3*3600) rianmax=max(10., min(15.0, big))
c        if(isec .le. 6*3600) rianmax=max(15., min(20.0, big))
        rianmax=20.
       do 255 i=2,iles
         if (ri(i) .ge. rianmax) it(i)=3
  255    if (ri(i) .lt. 0.001) it(i)=1

      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c   ***  tao & simpson (1989) and tao et al (1993)   **************
c     ****   search for the largest value of w and cloud water
c      bigw=0.
c      do 10 k=2,kles
c      do 10 i=2,iles
c       xa=abs(w(i,k))
c       bigw=max (bigw,xa)
c   10 continue
c       wbig=max(300., min(300., 0.5*bigw), 0.25*bigw)
c      bigqc=0.
c      do 12 k=2,kles
c      do 12 i=2,iles
c       xa=qc(i,k)
c       bigqc=max (bigqc,xa)
c   12 continue
c       qcbig= max(0.50e-3, 0.5*bigqc)
c       qcbig1=max(1.00e-3, 0.5*bigqc)
       wbig=300.
       qcbig= 0.50e-3
       qcbig1=1.00e-3
cccccccccccccccccccccccccccccccccccccccccccccccc
         do 300 i=2,iles
          icloudw=0
          icloudl=0
          icloud=0
          icwwl=0
          do 30 k=2,kles
cc  ***   lower cloudy region   ****************
           if (tair(k) .ge. 0.0) then
            if (qc(i,k) .ge. qcbig) icloudw=1
            if (qc(i,k) .ge. qcbig1) icloudl=1
            if (w(i,k) .ge. wbig) icwwl=1
           else
cc  ***   middle and upper cloudy region
           endif
   30     continue
c   ***  middle-upper level w > 0.5*wbig m/s or low-level qc > 1 g/kg
          if (it(i).eq.1) then
            if (icwwl .eq. 1 .and. icloudw .eq. 1) it(i)=3
          endif
           if (it(i).eq.2) then
             if (icloudl .eq. 1) it(i)=3
           endif
  300   continue
c
         do 350 i=2,iles
          icloud=0
          do 35 k=2,kles
            rzm=rho(k)*1.e6
            y3(k)=rzm*(qc(i,k)+qi(i,k)+qr(i,k)+qs(i,k)+qg(i,k))
           if (y3(k) .ge. .01) icloud=1
   35     continue
          if (it(i).eq.1 .and. icloud.eq.1) jv(i)=1
  350   continue
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c   ***   ********************************************************
c   ***   it(i)=3  convective region
c   ***   it(i)=2  stratiform region
c   ***   it(i)=4  stratiform region but no sfc precipitation
c   ***   it(i)=1  cloud free region

       a1=0.
       a2=0.
       a3=0.
       a4=0.
       lconv=0
       lanvl=0
       lnspt=0
       do 900 i=2,iles
        if (it(i) .eq. 2) then
          ics(i,3)=1
          a2=a2+ri(i)
          a4=a4+qr(i,2)
          lanvl=lanvl+1
        endif
        if (it(i) .eq. 3) then
          ics(i,2)=1
          a1=a1+ri(i)
          a3=a3+qr(i,2)
          lconv=lconv+1
        endif
        if (it(i) .eq. 1 .and. jv(i) .eq. 1) then
          ics(i,4)=1
          it(i)=4
          lnspt=lnspt+1
        endif
        iv(i)=it(i)
  900  continue
       aco=aco+a1
       aco1=aco1+a3
       aan=aan+a2
       aan1=aan1+a4
      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif
      return
      end
