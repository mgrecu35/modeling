cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine terp1 (n,x,f,w,y,int,tab,itab)
      dimension x(n),f(n),w(n),tab(3),itab(3)
      save
      call search (y,x,n,i)
      call interp (n,x,f,w,y,i,int,tab,itab)
      return
      end
