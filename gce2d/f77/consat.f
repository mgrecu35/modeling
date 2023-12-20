
      subroutine consat (rho)
c    (r&h)  specify some constants in satice routine   ******
      parameter (NZ=43,NT=2880,nb=10*nz+5*nt)
      common/iceopt/ ice913,ilif
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/rterv/ zrc,zgc,zsc,vrc0,vrc1,vrc2,vrc3,vgc,vsc
      common/rsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,rn2,bnd2,rn3,rn4,
     $  rn5,rn50,rn51,rn52,rn53,rn6,rn60,rn61,rn62,rn63,rn7,rn8,rn9,
     $  rn10,rn101,rn102,rn10a,rn10b,rn10c,rn11,rn12,rn12a(31),
     $  rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn171,rn172,rn17a,rn17b,
     $  rn17c,rn18,rn18a,rn19,rn191,rn192,rn19a,rn20,rn20a,rn20b,rn30,
     $  rn30a,rn21,bnd21,rn22,rn23,rn231,rn232,rn25,rn25a(31),rn31,beta,
     $  rn32,rn33,rn331,rn332,rn34,rn35
      common/brh1/ srr(nz),qrr(nz),z1(nb)
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      real    a1(31),a2(31),rho(nz)
      data a1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5,
     1   .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5,
     2   .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6,
     3   .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6,
     4   .6513e-6,.5956e-6,.5333e-6,.4834e-6/
      data a2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      save
c     ******************************************************************
      cpi=4.*atan(1.)
      cpi2=cpi*cpi
c      grvt=980.
      c76=7.66
      c358=35.86
      c172=17.26939
      c409=4098.026
      c218=21.87456
      c580=5807.695
      c610=6.1078e3
      c149=1.496286e-5
      c879=8.794142
      c141=1.4144354e7
c     ***************
        tca=2.43e3
        dwv=.226
        dva=1.718e-4
        amw=18.016
        ars=8.314e7
      t0=273.16
      t00=238.16
      alv=2.5e10
      alf=3.336e9
      als=2.8336e10
      avc=alv/cp
      afc=alf/cp
      asc=als/cp
      rw=4.615e6
      cw=4.187e7
      ci=2.093e7
c***   define the density and size distribution of precipitation
      roqr=1.
      tnw=.08
        roqs=.1
        tns=.16
cfred   tns=1.
          roqg=.4
cfred          tng=.04
        tng=.08
c***   define the coefficients used in terminal velocity
       ag=351.2
       bg=.37
          as=78.63154
          bs=.11
c         AS=152.93
c         BS=.25
         aw=2115.
         bw=.8
c
       if (ice913 .eq. 1) then
          t00=238.16
          tns=0.08
          tng=.04
          ag=372.3
          as=78.63154
          bs=.11
       endif
c
       bgh=.5*bg
       bsh=.5*bs
       bwh=.5*bw
       bgq=.25*bg
       bsq=.25*bs
       bwq=.25*bw
      ga3=2.
      ga4=6.
      ga5=24.
      ga6=120.
      ga7=720.
      ga8=5040.
      ga9=40320.
        ga4g=11.63177
        ga3g=3.3233625
        ga5gh=1.608355
        if(bg.eq.0.37) ga4g=9.730877
        if(bg.eq.0.37) ga3g=2.8875
        if(bg.eq.0.37) ga5gh=1.526425
          ga3d=2.54925
          ga4d=8.285063
          ga5dh=1.456943
          if(bs.eq.0.57) ga3d=3.59304
          if(bs.eq.0.57) ga4d=12.82715
          if(bs.eq.0.57) ga5dh=1.655588
          if(bs.eq.0.11) ga3d=2.218906
          if(bs.eq.0.11) ga4d=6.900796
          if(bs.eq.0.11) ga5dh=1.382792

cccccc             rutledge and hobbs, 1984   cccccccccccccccccccccccccc
          ga6d=144.93124
        ac1=as
        ac2=ag
       zrc=(cpi*roqr*tnw)**0.25
       zsc=(cpi*roqs*tns)**0.25
       zgc=(cpi*roqg*tng)**0.25
       vrc0=-26.7
       vrc1=20600./zrc
       vrc2=-204500./(zrc*zrc)
       vrc3=906000./(zrc*zrc*zrc)
       vsc=ac1*ga4d/(6.*zsc**bs)
       vgc=ac2*ga4g/(6.*zgc**bg)
cs      cd1=6.e-1
cs      cd2=4.*grvt/(3.*cd1)
cs     vgc=ga4g*sqrt(cd2*roqg/zgc)/6.
c     ****************************
      rn1=9.4e-15
      rn2=1.e-3
       bnd2=2.0e-3
c
       if (ice913 .eq. 1)  bnd2=1.5e-3
c
       esi=.1
      rn3=.25*cpi*tns*ac1*esi*ga3d
       esc=1.
      rn4=.25*cpi*esc*tns*ac1*ga3d
       eri=1.
CLIN   ERI=1.
       ERI=0.1
c
       if (ice913 .eq. 1)  ERI=1.
c
      rn5=.25*cpi*eri*tnw
       rn50=-.267e2*ga3
       rn51=5.15e3*ga4
       rn52=-1.0225e4*ga5
       rn53=7.55e3*ga6
       ami=1./(24.*6.e-9)
      rn6=cpi2*eri*tnw*roqr*ami
       rn60=-.267e2*ga6
       rn61=5.15e3*ga7
       rn62=-1.0225e4*ga8
       rn63=7.55e3*ga9
       esr=.5
c
       if (ice913 .eq. 1)  esr=1.
c
      rn7=cpi2*esr*tnw*tns*roqs
       esr=1.
      rn8=cpi2*esr*tnw*tns*roqr
       egs=.1
      rn9=cpi2*egs*tns*tng*roqs
      rn10=4.*tns
       rn101=.65
       rn102=.44*sqrt(ac1/dva)*ga5dh
       rn10a=alv*als*amw/(tca*ars)
       rn10b=alv/tca
       rn10c=ars/(dwv*amw)
      rn11=2.*cpi*tns*tca/alf
c HFO  
c	ami50=4.8e-7
c KFL (1993)
       ami50=4.8e-7*(100./50.)**3
       ami40=2.46e-7
       ami40=2.46e-7*.5**3
c       ami40=3.84e-9
       eiw=1.
       ui50=100.
c HFO
c	ri50=5.e-3
c KFL (1993)
       ri50=2.*5.e-3
       cmn=1.05e-15
      rn12=cpi*eiw*ui50*ri50*ri50
      do 10 k=1,31
        y1=1.-a2(k)
       rn13(k)=a1(k)*y1/(ami50**y1-ami40**y1)
       rn12a(k)=rn13(k)/ami50
       rn12b(k)=a1(k)*ami50**a2(k)
       rn25a(k)=a1(k)*cmn**a2(k)
   10 continue
       egc=1.
      rn14=.25*cpi*egc*ac2*tng*ga3g
       egi=.1
      rn15=.25*cpi*egi*tng*ac2*ga3g
       egi=1.
      rn15a=.25*cpi*egi*tng*ac2*ga3g
       egr=1.
      rn16=cpi2*egr*tng*tnw*roqr
      rn171=2.*cpi*tng*alv*dwv
       rn172=2.*cpi*tng*tca
       rn17a=.31*ga5gh*sqrt(ac2/dva)
       rn17b=cw-ci
       rn17c=cw
       apri=.66
       bpri=1.e-4
CTAO
       bpri=.5*bpri
c
       if (ice913 .eq. 1)  bpri=1.e-4
c
      rn18=20.*cpi2*bpri*tnw*roqr
       rn18a=apri
      rn19=2.*cpi*tng*tca/alf
       rn191=.78
       rn192=.31*ga5gh*sqrt(ac2/dva)
       rn19a=cw/alf
      rn20=2.*cpi*tng
       rn20a=als*als*amw/(tca*ars)
       rn20b=als/tca
      rn30=2.*cpi*tng
       rn30a=alv*alv*amw/(tca*ars)
      rn21=1.e-3
       bnd21=1.5e-3
       erc=1.
      rn22=.25*cpi*erc*tnw
      rn23=2.*cpi*tnw
       rn231=.78
       rn232=.31*ga3*sqrt(3.e3/dva)
       cn0=1.e-8

c      rn25=cn0	!tao-new
       rn25=cn0

      rn31=1.e-17
       beta=-.6
      rn32=4.*51.545e-4
      rn33=4.*tns
       rn331=.65
       rn332=.44*sqrt(ac1/dva)*ga5dh
CTAO
       esc=0.5
c
       if (ice913 .eq. 1)  esc=1.
c
       amc=1./(24.*4.e-9)
      rn34=cpi2*esc*amc*ac1*roqs*tns*ga6d
      rn35=alv*alv/(cp*rw)
c     ****************************
      do 20 k=1,nz
       srr(k)=1./sqrt(rho(k))
   20  qrr(k)=sqrt(srr(k))
      return
      end
