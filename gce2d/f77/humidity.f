Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine humidity(T,P,Qsat)                                 
c
c     Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
c
      Qsat = (1.0007+3.46e-6*P)*6.1121*exp(17.502*T/(240.97+T))
      return
      end
