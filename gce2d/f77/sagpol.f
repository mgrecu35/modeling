
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sagpol(tau,ssc,g0,rll,tll)

c*********************************************************************
c-----transmittance (tll) and reflectance (rll) of diffuse radiation
c     follows Sagan and Pollock (JGR, 1967).
c     also, eq.(31) of Lacis and Hansen (JAS, 1974).
c
c-----input parameters:
c
c      tau: the effective optical thickness
c      ssc: the effective single scattering albedo
c      g0:  the effective asymmetry factor
c
c-----output parameters:
c
c      rll: the layer reflection of diffuse radiation
c      tll: the layer transmission of diffuse radiation
c
c*********************************************************************

      implicit none

      real    one,three,four
      parameter (one=1., three=3., four=4.)

c-----output parameters:

      real    tau,ssc,g0,expmn

c-----output parameters:

      real    rll,tll

c-----temporary arrays

      real    xx,uuu,ttt,emt,up1,um1,st1

      save

          if(ssc .gt. 0.001) then

             xx  = one-ssc*g0
             uuu = sqrt( xx/(one-ssc))
             ttt = sqrt( xx*(one-ssc)*three )*tau
             emt = expmn(-ttt)
             up1 = uuu + one
             um1 = uuu - one
             xx  = um1*emt
             st1 = one / ((up1+xx) * (up1-xx))
             rll = up1*um1*(one-emt*emt)*st1
             tll = uuu*four*emt         *st1

          else

             rll = 0.0
             tll = expmn(-1.66*tau)

          endif

      return
      end
