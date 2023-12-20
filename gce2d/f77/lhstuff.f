      subroutine initlh()
      parameter (NX=514,NZ=43,NT=2880,ITT=244)
      common /lhblock/ rlh(nx,nz)
      do i=1,nx
         do k=1,nz
            rlh(i,k)=0.
         enddo
      enddo
      end

      subroutine average(dtint)
      parameter (NX=514,NZ=43,NT=2880,ITT=244)
      common /lhblock/ rlh(nx,nz)
      do i=1,nx
         do k=1,nz
            rlh(i,k)=rlh(i,k)/dtint
         enddo
      enddo
      end
