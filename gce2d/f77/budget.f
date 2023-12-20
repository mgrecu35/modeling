cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine budget (taa)

      implicit none
      integer nx,nz,itt
      parameter (NX=514,NZ=43,ITT=244)

      real    taa(nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1
      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    smf(nz),smu(nz),smd(nz),smn(nz),stf(nz),stu(nz),
     $   std(nz),stn(nz),sqf(nz),squ(nz),sqd(nz),sqn(nz),sdt(nz),
     $   sdq(nz),stv(nz),
     $   stf1(nz),stu1(nz),std1(nz),stn1(nz),sqf1(nz),squ1(nz),
     $   sqd1(nz),sqn1(nz),stf2(nz),stu2(nz),std2(nz),stn2(nz),sqf2(nz),
     $   squ2(nz),sqd2(nz),sqn2(nz),cld(nz),cldu(nz),cldd(nz),cldn(nz)
      common/b8/ smf,smu,smd,smn,stf,stu,std,stn,sqf,squ,sqd,sqn,sdt,
     $           sdq,stv,stf1,stu1,std1,stn1,sqf1,squ1,
     $           sqd1,sqn1,stf2,stu2,std2,stn2,sqf2,
     $           squ2,sqd2,sqn2,cld,cldu,cldd,cldn

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
      real    thom(nz,4,7),tdw(nz,4,7),tmlt(nz,4,7),saut(nz,4,7),
     $ saci(nz,4,7),sacw(nz,4,7),raci(nz,4,7),tacr(nz,4,7),raut(nz,4,7),
     $ racw(nz,4,7),sfw(nz,4,7),sfi(nz,4,7),gacs(nz,4,7),gacw(nz,4,7),
     $ gaci(nz,4,7),gacr(nz,4,7),gwet(nz,4,7),gaut(nz,4,7),racs(nz,4,7),
     $ sacr(nz,4,7),gfr(nz,4,7),smlt(nz,4,7),gmlt(nz,4,7),sdep(nz,4,7),
     $ ssub(nz,4,7),gsub(nz,4,7),pern(nz,4,7),d3ri(nz,4,7),d3ir(nz,4,7),
     $ d2sr(nz,4,7),d2rs(nz,4,7),gdry(nz,4,7),coc(nz,4,7),coe(nz,4,7),
     $ smf0(nz,4,7),qc0(nz,4,7),qr0(nz,4,7),qi0(nz,4,7),qs0(nz,4,7),
     $ qg0(nz,4,7),sqc0(nz,4,7),sqr0(nz,4,7),sqi0(nz,4,7),sqs0(nz,4,7),
     $ sqg0(nz,4,7),erns(nz,4,7),wgrs(nz,4,7),qsws(nz,4,7),
     $ tb0(nz,4),qb0(nz,4)
      common/bsts/ thom,tdw,tmlt,saut,saci,sacw,raci,tacr,raut,
     $             racw,sfw,sfi,gacs,gacw,gaci,gacr,gwet,gaut,racs,
     $             sacr,gfr,smlt,gmlt,sdep,ssub,gsub,pern,d3ri,d3ir,
     $             d2sr,d2rs,gdry,coc,coe,smf0,qc0,qr0,qi0,qs0,
     $             qg0,sqc0,sqr0,sqi0,sqs0,sqg0,erns,wgrs,qsws,
     $             tb0,qb0

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    q1(nz,4,7),q2(nz,4,7),sub(nz,4,7),otr(nz,4,7),pir(nz),
     $        fmlt(nz,4,7),evap(nz,4,7),cond(nz,4,7),dep(nz,4,7)
      real    d1(nz),d2(nz),d3(nz),d4(nz),d5(nz)

      integer iaminf,k,k1,m,l
      real    a00,a1,a11
      real    afc,alf,als,alv,asc,avc,bs1,bs2,bs3,bs4
      real    cp,dtday,dtmas,period,pird,pird1

      save

c
       m=4
       period=6.
        iaminf=360
       a00=60./20.
       a11=60.
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.004e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
       dtday=1./period
       dtmas=1./period
      do 10 k=1,nz
  10   pir(k)=pi(k)/512.
      do 15 k=2,nz-1
        pird=dtday*pir(k)*am(k)/(rho(k)*dz)
        pird1=dtday*am(k)/(rho(k)*dz*512.)
       sft(k)=sft(k)*pi(k)*dtday
       sfq(k)=-sfq(k)
cc     sfq(k)=-avc*sfq(k)*dtday
        d1(k)=(stf1(k)-stf1(k+1))*pird*a11
        d2(k)=-avc*(sqf1(k)-sqf1(k+1))*pird1*a11
         d3(k)=tls(k)*pi(k)*dtday
cc15     d4(k)=avc*qls(k)*dtday
         d4(k)=qls(k)
  15     d5(k)=qa(k)-taa(k)
      do 100 l=1,4
      do 100 k=1,nz
        pird=dtday/512.
       cond(k,l,m)=avc*coc(k,l,m)*a00
       evap(k,l,m)=-avc*coe(k,l,m)*a00
      bs1=avc*(coc(k,l,m)-coe(k,l,m))*a00
       dep(k,l,m)=asc*(gaut(k,l,m)+sdep(k,l,m)+gsub(k,l,m))*a00
       sub(k,l,m)=-asc*(pern(k,l,m))*a00
      bs2=asc*(gaut(k,l,m)-pern(k,l,m)+sdep(k,l,m)+gsub(k,l,m)-
     1    gwet(k,l,m)-gfr(k,l,m))*a00
       fmlt(k,l,m)=-afc*(smlt(k,l,m)+gmlt(k,l,m))*a00
      bs3=afc*(smlt(k,l,m)+gmlt(k,l,m))*a00
       otr(k,l,m)=afc*(sacw(k,l,m)+tacr(k,l,m)+sfw(k,l,m)+tmlt(k,l,m)+
     1            gacw(k,l,m)+gacr(k,l,m)+sacr(k,l,m)-d2rs(k,l,m))*a11
     2           -asc*(gwet(k,l,m)+gfr(k,l,m))*a00
      bs4=afc*(sacw(k,l,m)+tacr(k,l,m)+sfw(k,l,m)+tmlt(k,l,m)+
     1    gacw(k,l,m)+gacr(k,l,m)+sacr(k,l,m)-d2rs(k,l,m))*a11
       q1(k,l,m)=(bs1+bs2+bs4-bs3)*pird
cc     q2(k,l,m)=(bs1+bs2)*pird
       q2(k,l,m)=(bs1/avc+bs2/asc)/512.
       smf0(k,l,m)=smf0(k,l,m)*dtmas
        cond(k,l,m)=cond(k,l,m)*pird
        evap(k,l,m)=evap(k,l,m)*pird
        dep(k,l,m)=dep(k,l,m)*pird
        sub(k,l,m)=sub(k,l,m)*pird
        fmlt(k,l,m)=fmlt(k,l,m)*pird
        otr(k,l,m)=otr(k,l,m)*pird
  100  continue
       write(6,601) iaminf,m
      do 200 k1=1,nz
       k=nz+1-k1
  200 write(6,602) k,q1(k,1,m),q1(k,2,m),q1(k,3,m),q2(k,1,m),q2(k,2,m),
     1  q2(k,3,m),smf0(k,1,m),smf0(k,2,m),smf0(k,3,m)
       write(6,601) iaminf,m
      do 250 k1=1,nz
       k=nz+1-k1
  250 write(6,602) k,cond(k,1,m),cond(k,2,m),cond(k,3,m),evap(k,1,m),
     1  evap(k,2,m),evap(k,3,m),fmlt(k,1,m),fmlt(k,2,m),fmlt(k,3,m)
       write(6,601) iaminf,m
      do 270 k1=1,nz
       k=nz+1-k1
  270 write(6,602) k,dep(k,1,m),dep(k,2,m),dep(k,3,m),sub(k,1,m),
     1  sub(k,2,m),sub(k,3,m),otr(k,1,m),otr(k,2,m),otr(k,3,m)
      do 275 k1=1,nz
       k=nz+1-k1
       a1=-q2(k,1,4)-sfq(k)+d4(k)
  275 write(6,602) k,sfq(k),d4(k),a1,d5(k),d1(k),d2(k),d3(k),sft(k),
     1  sv(k)
  601 format(2x,'lev',4x,'q1(l=1)',4x,'q1(l=2)',4x,'q1(l=3)',4x,
     1 'q2(l=1)',4x,'q2(l=2)',4x,'q2(l=3)',4x,'mass(1)',4x,
     2 'mass(2)',4x,'mass(3)',10x,i4,'min',3x,i2,/)
  602 format(1x,i4,1x,9(1x,e10.3))
      return
      end
