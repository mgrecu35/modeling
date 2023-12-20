cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp (n,x,f,w,y,i,int,tab,itab)
      dimension x(n),f(n),w(n),tab(3),itab(3)
      save
      ii(i)=(i-1)*int+1
      flk=x(i+1)-x(i)
      flp=x(i+1)-y
      fl0=y-x(i)
      i0=ii(i)
      ip=i0+int
      if(itab(1).ne.1) go to 102
      a=(w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)
      b=(f(ip)/flk-w(ip)*flk/6.)*fl0
      c=(f(i0)/flk-w(i0)*flk/6.)*flp
      tab(1)=a+b+c
  102 if(itab(2).ne.1) go to 104
      a=(w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
      b=(f(ip)-f(i0))/flk
      c=(w(i0)-w(ip))*flk/6.
      tab(2)=a+b+c
  104 if(itab(3).ne.1) go to 106
      tab(3)=(w(i0)*flp+w(ip)*fl0)/flk
  106 return
      end
