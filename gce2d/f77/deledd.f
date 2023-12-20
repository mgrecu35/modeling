
c*********************************************************************

      subroutine deledd(tau,ssc,g0,csm,rr,tt,td)

c*********************************************************************
c
c-----uses the delta-eddington approximation to compute the
c     bulk scattering properties of a single layer
c     coded following King and Harshvardhan (JAS, 1986)
c
c  inputs:
c
c     tau: the effective optical thickness
c     ssc: the effective single scattering albedo
c     g0:  the effective asymmetry factor
c     csm: the effective secant of the zenith angle
c
c  outputs:
c
c     rr: the layer reflection of the direct beam
c     tt: the layer diffuse transmission of the direct beam
c     td: the layer direct transmission of the direct beam
c
c*********************************************************************

      implicit none

      real zero,one,two,three,four,fourth,seven,thresh
      parameter (one =1., three=3.)
      parameter (two =2., seven=7.)
      parameter (four=4., fourth=.25)
      parameter (zero=0., thresh=1.e-8)
      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----input parameters
      real tau,ssc,g0,csm

c-----output parameters
      real rr,tt,td

c-----temporary parameters

      real zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,
     *     all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4
      save

c---------------------------------------------------------------------

                zth = one / csm

c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  K & H eqs(27-29)

                ff  = g0*g0
                xx  = one-ff*ssc
                taup= tau*xx
                sscp= ssc*(one-ff)/xx
                gp  = g0/(one+g0)

c  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.

                xx  =  three*gp
                gm1 =  (seven - sscp*(four+xx))*fourth
                gm2 = -(one   - sscp*(four-xx))*fourth

c  akk is k as defined in eq(25) of K & H

                akk = sqrt((gm1+gm2)*(gm1-gm2))

                xx  = akk * zth
                st7 = one - xx
                st8 = one + xx
                st3 = st7 * st8

                if (abs(st3) .lt. thresh) then
                    zth = zth + 0.001
                    xx  = akk * zth
                    st7 = one - xx
                    st8 = one + xx
                    st3 = st7 * st8
                endif

c  extinction of the direct beam transmission

                td  = exp(-taup/zth)

c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H

                gm3  = (two - zth*three*gp)*fourth
                xx   = gm1 - gm2
                alf1 = gm1 - gm3 * xx
                alf2 = gm2 + gm3 * xx

c  all is last term in eq(21) of K & H
c  bll is last term in eq(22) of K & H

                xx  = akk * two
                all = (gm3 - alf2 * zth    )*xx*td
                bll = (one - gm3 + alf1*zth)*xx

                xx  = akk * gm3
                cll = (alf2 + xx) * st7
                dll = (alf2 - xx) * st8

                xx  = akk * (one-gm3)
                fll = (alf1 + xx) * st8
                ell = (alf1 - xx) * st7

                st2 = exp(-akk*taup)
                st4 = st2 * st2

                st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

c  rr is r-hat of eq(21) of K & H
c  tt is diffuse part of t-hat of eq(22) of K & H

                rr =   ( cll-dll*st4    -all*st2)*st1
                tt = - ((fll-ell*st4)*td-bll*st2)*st1

                rr = max(rr,zero)
                tt = max(tt,zero)

      return
      end
