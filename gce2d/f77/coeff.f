cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeff (n,x,f,w,iop,int,wk)
c-----------------------------------------------------------------------
c--- Compute the coefficients for the cubic spline interpolation     ---
c-----------------------------------------------------------------------
c--- Inputs:  x --- pressure at the model level                      ---
c---          f --- variables being interpolated                     ---
c--- Outputs: w --- coefficients                                     ---
c---          wk--- stores couple of coefficients at each point      ---
c-----------------------------------------------------------------------
      dimension x(n),f(n),w(n),wk(n,4),iop(2)
      save
      ii(i)=(i-1)*int+1
      j0=1
      do i=2,n
        jm=j0
        j0=j0+int
        wk(i,1)=x(i)-x(i-1)
        wk(i,2)=(f(j0)-f(jm))/wk(i,1)
        wk(i,3)=wk(i,1)/6.
        wk(i,1)=wk(i,1)/3.
      end do
      nn=n
      mk=iop(1)
      ml=iop(2)
      go to (102,103,104,105) ,mk
  102 continue
      wk(2,2)=wk(3,2)-wk(2,2)-wk(2,3)*w(1)
      wk(2,3)=0.
      wk(2,1)=wk(2,1)+wk(3,1)
      i1=2
      nn=nn-1
      go to 106
  103 continue
      wk(1,2)=wk(2,2)-w(1)
      wk(2,2)=wk(3,2)-wk(2,2)
      wk(1,3)=0.
      wk(1,1)=wk(2,1)
      wk(2,1)=wk(2,1)+wk(3,1)
      i1=1
      go to 106
  104 continue
      y2=wk(2,2)
      b2=wk(2,1)
      wk(2,2)=wk(3,2)-wk(2,2)
      wk(2,1)=wk(3,1)+wk(2,1)
      i1=2
      nn=nn-1
      go to 106
  105 continue
      a12=x(1)-x(2)
      a13=x(1)-x(3)
      a14=x(1)-x(4)
      a23=x(2)-x(3)
      a24=x(2)-x(4)
      a34=x(3)-x(4)
      j1=1
      j2=j1+int
      j3=j2+int
      j4=j3+int
      w(1)=(1./a12+1./a13+1./a14)*f(j1)-
     1     a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2     a12*a13/(a14*a24*a34)*f(j4)
      go to 103
  106 continue
      i2=n-2
      do i=3,i2
        wk(i,2)=wk(i+1,2)-wk(i,2)
        wk(i,1)=wk(i+1,1)+wk(i,1)
      end do
      in=ii(n)
      go to (108,109,110,111) ,ml
  108 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
      wk(n,3)=0.
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      nn=nn-1
      go to 112

  109 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)
      wk(n,2)=-wk(n,2)+w(in)
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      wk(1,4)=0.
      go to 112

  110 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)
      wk(n,2)=y2-wk(n,2)
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      wk(n,1)=wk(n,1)+b2
      wk(1,4)=wk(2,3)
      go to 112

  111 continue
      a12=x(n)-x(n-1)
      a13=x(n)-x(n-2)
      a14=x(n)-x(n-3)
      a23=x(n-1)-x(n-2)
      a24=x(n-1)-x(n-3)
      a34=x(n-2)-x(n-3)
      j1=in
      j2=j1-int
      j3=j2-int
      j4=j3-int
      w(in)=(1./a12+1./a13+1./a14)*f(j1)-
     1      a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2      a12*a13/(a14*a24*a34)*f(j4)
      go to 109

  112 continue
      ii1=ii(i1)
      call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(ii1),int)
      go to (114,114,113,114) ,mk

  113 continue
      w(1)=w(in)
  114 continue
      return
      end
