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