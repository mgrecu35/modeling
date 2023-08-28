subroutine fscale(a,amean,astd,ascale,nx,ny)
    implicit none
    integer :: nx,ny, i,j,it
    real :: a(nx,72,ny),amean(72),astd(72)
    real, intent(out) :: ascale(nx*ny,72)
    it=1
    do i=1,nx
        do j=1,ny
            ascale(it,:)=(a(i,:,j)-amean(:))/astd(:)
            it=it+1
        end do
    end do

end subroutine fscale

subroutine fscale1(a,amean,astd,ascale,nx,ny,nz)
    implicit none
    integer :: nx,ny,nz,i,j,k,it
    real :: a(nx,nz,ny),amean,astd
    real, intent(out) :: ascale(nx*ny*nz)
    it=1
    do i=1,nx
        do j=1,ny
            do k=1,nz
                ascale(it)=(a(i,k,j)-amean)/astd
                it=it+1
            end do
        end do
    end do

end subroutine fscale1

subroutine fscale1d(a,amean,astd,ascale,nx,ny,nz)
    implicit none
    integer :: nx,ny,nz,i,j,k,it
    real :: a(nx,nz,ny),amean,astd
    real, intent(out) :: ascale(nx*ny,nz)
    it=1
    do i=1,nx
        do j=1,ny
            ascale(it,:)=(a(i,:,j)-amean)/astd
            it=it+1
        end do
    end do

end subroutine fscale1d

subroutine update_var(a_in,da,dt,amean,astd,a_out,nx,ny)
    implicit none
    integer :: nx,ny, i,j,it
    real :: a_in(nx,72,ny),amean(72),astd(72),da(nx*ny,72),dt
    real, intent(out) :: a_out(nx,72,ny)
    it=1
    do i=1,nx
        do j=1,ny
            a_out(i,:,j)=a_in(i,:,j)+(da(it,:)*astd(:)+amean(:))*dt
            it=it+1
        end do
    end do

end subroutine update_var


subroutine update_var1(a_in,da,dt,amean,astd,a_out,nx,ny,nz)
    implicit none
    integer :: nx,ny, nz, i,j,k,it
    real :: a_in(nx,nz,ny),amean,astd,da(nx*ny*nz),dt
    real, intent(out) :: a_out(nx,nz,ny)
    it=1
    do i=1,nx
        do j=1,ny
            do k=1,nz
                a_out(i,k,j)=a_in(i,k,j)+(da(it)*astd+amean)*dt
                it=it+1
            end do
        end do
    end do

end subroutine update_var1

subroutine update_var1d(a_in,da,dt,amean,astd,a_out,nx,ny,nz)
    implicit none
    integer :: nx,ny, nz, i,j,k,it
    real :: a_in(nx,nz,ny),amean,astd,da(nx*ny,nz),dt
    real, intent(out) :: a_out(nx,nz,ny)
    it=1
    do i=1,nx
        do j=1,ny
            a_out(i,:,j)=a_in(i,:,j)+(da(it,:)*astd+amean)*dt
            it=it+1
        end do
    end do

end subroutine update_var1d