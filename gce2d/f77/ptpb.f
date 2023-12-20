
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ptpb (m1)
      parameter (NX=514,NZ=43,NT=2880,nxi=nx-2,nzm1=nz-1)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1 st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b1c/ qcl(nx,nz)
      common/b1r/ qrn(nx,nz)
      common/b1i/ qci(nx,nz)
      common/b1s/ qcs(nx,nz)
      common/b1g/ qcg(nx,nz)
      common/bw1/ rst(nt,nxi),pcltop(nt,nxi),pclbot(nt,nxi)
      dimension tqe(nz)
      save
      do 100 i=2,nx-1
      ix=i-1
       do k=2,nz-1
         tqe(k)=qcl(i,k)+qrn(i,k)+qci(i,k)+qcs(i,k)+qcg(i,k)
       end do
         call ttop(tqe,ktop,kbase)
          pcltop(m1,ix)=1.e-3*p0(ktop+1)
            if ( ktop.eq.nzm1) pcltop(m1,ix)=1.e-3*p0(ktop)
          pclbot(m1,ix)=1.e-3*p0(kbase)
            if (kbase.eq. 2  ) pclbot(m1,ix)=1.e-3*p0(kbase)
           thin=pclbot(m1,ix)-pcltop(m1,ix)
            if(thin.lt.50.) go to 100
  100 continue
      return
      end          
