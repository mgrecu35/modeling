subroutine read_data(fname, nx,ny, fields)
  implicit none
  integer, intent(in) :: nx,ny
  real, intent(out) :: fields(16,nx,72,ny)
  character(len=*), intent(in) :: fname
  integer :: i
  print*, fname
  open(10, file=fname, form='unformatted',convert='big_endian')
  do i=1,16
    read(10) fields(i,:,:,:)
  end do
  close(10)
end subroutine read_data

subroutine init_keras()
    use keras_def
    call net % load("/Users/mgrecu/cmbv7_T/FKB/mp_model_0d.txt")
    print*, 'init keras'
    
    input = [ -6.06659576E-02,  -4.95505370E-02, -0.110201314,&
              -0.138018325 ,     -8.25610980E-02,  0.351970166,&
               -0.651567221,     -0.712960541 ,      4.20613103E-02]
  
    ! run test input through network
    result = net % output(input)
    print *, result
  
end subroutine init_keras

subroutine updatefields(qv,qc,qr,qi,qs,qg,th,press,dz,&
    qv_new,qc_new,qr_new,qi_new,qs_new,qg_new,th_new,nx,ny,nz)
use keras_def
real :: xm(9), xs(9)
real :: ym(7), ys(7)
real :: qv(nx,nz,ny), qc(nx,nz,ny), qr(nx,nz,ny), qi(nx,nz,ny), qs(nx,nz,ny), qg(nx,nz,ny), &
    th(nx,nz,ny), press(nx,nz,ny), dz(nx,nz,ny)
real,intent(out) :: qv_new(nx,nz,ny), qc_new(nx,nz,ny), qr_new(nx,nz,ny), qi_new(nx,nz,ny), &
    qs_new(nx,nz,ny), qg_new(nx,nz,ny), th_new(nx,nz,ny)
integer, intent(in) :: nx,ny,nz
integer :: i,j,k
real :: temp
xm=[0.003315005218148738, 1.3991694213505624e-05, 1.460830579877713e-05, 2.4933441397318293e-05, &
4.743995609906075e-05, 0.0001070590156169587, 248.90307223405165, 42325.43806454194, 245.6880410989276]
xs=[0.005188289465698395, 0.00012867134045073828, 0.00016844543890582057, 0.00012484716493930284, &
0.00019047270694516342, 0.0007185409687322846, 31.769377267665877, 26040.11419624085, 52.57923401164694]

ym=[-9.800165435953876e-07, -1.5075518595215102e-07, -1.7199690209827199e-07, 5.8635545128302734e-08, &
1.323342222114488e-07, 1.299709084929269e-07, 0.002767458625174922]
ys=[1.6866934204044486e-05, 1.72771219177893e-05, 2.402103196794837e-05, 1.0752144470956599e-05, &
3.585489061213833e-06, 3.0273808704202674e-05, 0.04689521341170341]

!inputVars=["qv","qc","qr","qi","qs","qg","temp","press","dz"]
!outputVars=["qv_tend","qc_tend","qr_tend","qi_tend","qs_tend","qg_tend","temp_tend"]
print*, nx,ny,nz

do i=1,nx
    do k=1,nz
        do j=1,ny
            input(1)=(qv(i,k,j)-xm(1))/xs(1)
            input(2)=(qc(i,k,j)-xm(2))/xs(2)
            input(3)=(qr(i,k,j)-xm(3))/xs(3)
            input(4)=(qi(i,k,j)-xm(4))/xs(4)
            input(5)=(qs(i,k,j)-xm(5))/xs(5)
            input(6)=(qg(i,k,j)-xm(6))/xs(6)
            temp=th(i,k,j)*(press(i,k,j)/100000)**(287.0/1004.0)
            input(7)=(temp-xm(7))/xs(7)
            input(8)=(press(i,k,j)-xm(8))/xs(8)
            input(9)=(dz(i,k,j)-xm(9))/xs(9)
            result=net%output(input)
            qv_new(i,k,j)=qv(i,k,j)+result(1)*ys(1)+ym(1)
            qc_new(i,k,j)=qc(i,k,j)+result(2)*ys(2)+ym(2)
            qr_new(i,k,j)=qr(i,k,j)+result(3)*ys(3)+ym(3)
            qi_new(i,k,j)=qi(i,k,j)+result(4)*ys(4)+ym(4)
            qs_new(i,k,j)=qs(i,k,j)+result(5)*ys(5)+ym(5)
            qg_new(i,k,j)=qg(i,k,j)+result(6)*ys(6)+ym(6)
            th_new(i,k,j)=th(i,k,j)+(result(7)*ys(7)+ym(7))*(1e5/press(i,k,j))**(287.0/1004.0)
        enddo   
    enddo
enddo   
end subroutine updatefields
