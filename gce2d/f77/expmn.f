
c*******************************************************************

      function expmn(fin)

      implicit none
c*******************************************************************
c compute exponential for arguments in the range 0> fin > -10.

      real    one,expmin,e1,e2,e3,e4
      parameter (one=1.0, expmin=-10.0)
      parameter (e1=1.0,        e2=-2.507213e-1)
      parameter (e3=2.92732e-2, e4=-3.827800e-3)

      real    fin,expmn

      if (fin .lt. expmin) fin = expmin
      expmn = ((e4*fin + e3)*fin+e2)*fin+e1
      expmn = expmn * expmn
      expmn = one / (expmn * expmn)

      return
      end
