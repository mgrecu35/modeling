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
    use  MODULE_MP_morr_two_moment2

    call morr_two_moment_init2
   
    call net % load("/Users/mgrecu/cmbv7_T/FKB/mp_model_0d_reset_log.txt")
    do i=1,10
        call net_tile(i)%load("/Users/mgrecu/cmbv7_T/FKB/mp_model_0d_reset_log.txt")
    enddo
    print*, 'init keras'
    
    !input = [ -6.06659576E-02,  -4.95505370E-02, -0.110201314,&
    !          -0.138018325 ,     -8.25610980E-02,  0.351970166,&
    !           -0.651567221,     -0.712960541 ,      4.20613103E-02]
    if(.not.allocated(input)) allocate(input(10,9))
    if(.not.allocated(result)) allocate(result(10,7))
    ! run test input through network
    !result = net % output(input)
    !print *, result
    print*, 'and out'
    !stop
end subroutine init_keras

subroutine updatefields(qv,qc,qr,qi,qs,qg,th,press,dz,&
    qv_new,qc_new,qr_new,qi_new,qs_new,qg_new,th_new,nx,ny,nz,itile)
use keras_def
real :: xm(9), xs(9)
real :: ym(7), ys(7)
real :: qv(nx,nz,ny), qc(nx,nz,ny), qr(nx,nz,ny), qi(nx,nz,ny), qs(nx,nz,ny), qg(nx,nz,ny), &
    th(nx,nz,ny), press(nx,nz,ny), dz(nx,nz,ny)
real,intent(out) :: qv_new(nx,nz,ny), qc_new(nx,nz,ny), qr_new(nx,nz,ny), qi_new(nx,nz,ny), &
    qs_new(nx,nz,ny), qg_new(nx,nz,ny), th_new(nx,nz,ny)
integer, intent(in) :: nx,ny,nz,itile
integer :: i,j,k
real :: p1
real :: temp
xm=[0.003315005218148738, 1.3991694213505624e-05, 1.460830579877713e-05, 2.4933441397318293e-05, &
4.743995609906075e-05, 0.0001070590156169587, 248.90307223405165, 42325.43806454194, 245.6880410989276]
xs=[0.005188289465698395, 0.00012867134045073828, 0.00016844543890582057, 0.00012484716493930284, &
0.00019047270694516342, 0.0007185409687322846, 31.769377267665877, 26040.11419624085, 52.57923401164694]

ym=[-9.800165435953876e-07, -1.5075518595215102e-07, -1.7199690209827199e-07, 5.8635545128302734e-08, &
1.323342222114488e-07, 1.299709084929269e-07, 0.002767458625174922]
ys=[1.6866934204044486e-05, 1.72771219177893e-05, 2.402103196794837e-05, 1.0752144470956599e-05, &
3.585489061213833e-06, 3.0273808704202674e-05, 0.04689521341170341]

xm=[0.002960509017827147, 1.0105936096235195e-05, 3.7480526540547464e-06, 1.547370552365313e-05, &
    1.2546235110746174e-05, 3.7969343019402076e-05, 216.59745474514978, 40254.284540731016, 245.9760871780734]
xs=[0.004918123213399891, 0.00012930930327912195, 9.479519751428191e-05, 0.00012595858643891413, &
    0.0001022809752772555, 0.00040436589821293107, 61.04199128710109, 26998.832965656176, 51.82192204734818]
ym=[-7.600808508934862e-07, -8.667786060777192e-08, -7.795199054679282e-10, 7.650110443948787e-08, &
    1.8685398142641414e-07, 3.1550781178184074e-07, 0.0021556757216830307]
ys=[1.2594798895898411e-05, 1.746612892371693e-05, 2.7133056827293538e-05, 9.670111399258711e-06, &
    3.4429783954797836e-06, 2.740044287612875e-05, 0.0351797748224097]



xm=[1.492011479793688, 0.025988618680400595, 0.011208313342917615, 0.05184656346743224, 0.04663025815558111, &
0.041444475973150575, 216.55349607150683, 40258.362056095444, 245.98054906652368]
xs=[1.0903266506687994, 0.20222921099576746, 0.11934684338003179, 0.2683834080121284, 0.25375568515900776, &
0.28448053420969976, 61.05797285364713, 26999.847063259014, 51.81974444364091]
ym=[-0.00022519540768959554, -0.0022948285191311444, -0.00022231164735151342, -0.00027946231809359835, &
0.00039113425244842015, 0.0002680951979741824, 0.0021518989187387623]
ys=[0.004173149648066692, 0.034305631040542114, 0.017506699620924254, 0.014724501474699588, 0.00760419049106299,&
 0.010852795422884587, 0.03514906199702389]

!inputVars=["qv","qc","qr","qi","qs","qg","temp","press","dz"]
!outputVars=["qv_tend","qc_tend","qr_tend","qi_tend","qs_tend","qg_tend","temp_tend"]
!print*, nx,ny,nz

do i=1,nx
    do k=1,nz
        do j=1,ny
            if(qv(i,k,j)<0)qv(i,k,j)=0
            if(qc(i,k,j)<0)qc(i,k,j)=0
            if(qr(i,k,j)<0)qr(i,k,j)=0
            if(qi(i,k,j)<0)qi(i,k,j)=0
            if(qs(i,k,j)<0)qs(i,k,j)=0
            if(qg(i,k,j)<0)qg(i,k,j)=0
            input(itile,1)=(log10(1+qv(i,k,j)/1e-5)-xm(1))/xs(1)
            input(itile,2)=(log10(1+qc(i,k,j)/1e-5)-xm(2))/xs(2)
            input(itile,3)=(log10(1+qr(i,k,j)/1e-5)-xm(3))/xs(3)
            input(itile,4)=(log10(1+qi(i,k,j)/1e-5)-xm(4))/xs(4)
            input(itile,5)=(log10(1+qs(i,k,j)/1e-5)-xm(5))/xs(5)
            input(itile,6)=(log10(1+qg(i,k,j)/1e-5)-xm(6))/xs(6)
            temp=th(i,k,j)*(press(i,k,j)/100000)**(287.0/1004.0)
            input(itile,7)=(temp-xm(7))/xs(7)
            input(itile,8)=(press(i,k,j)-xm(8))/xs(8)
            input(itile,9)=(dz(i,k,j)-xm(9))/xs(9)
            result(itile,:)=net_tile(itile)%output(input(itile,:))
            !qv_new(i,k,j)=qv(i,k,j)+result(itile,1)*ys(1)+ym(1)
            !qc_new(i,k,j)=qc(i,k,j)+result(itile,2)*ys(2)+ym(2)
            !qr_new(i,k,j)=qr(i,k,j)+result(itile,3)*ys(3)+ym(3)
            !qi_new(i,k,j)=qi(i,k,j)+result(itile,4)*ys(4)+ym(4)
            !qs_new(i,k,j)=qs(i,k,j)+result(itile,5)*ys(5)+ym(5)
            !qg_new(i,k,j)=qg(i,k,j)+result(itile,6)*ys(6)+ym(6)
            p1=log10(1+qv(i,k,j)/1e-5)+result(itile,1)*ys(1)+ym(1)
            if(p1<0) p1=0
            qv_new(i,k,j)=(10**p1-1)*1e-5
            p1=log10(1+qc(i,k,j)/1e-5)+result(itile,2)*ys(2)+ym(2)
            if(p1<0) p1=0
            qc_new(i,k,j)=(10**p1-1)*1e-5
            p1=log10(1+qr(i,k,j)/1e-5)+result(itile,3)*ys(3)+ym(3)
            if(p1<0) p1=0
            qr_new(i,k,j)=(10**p1-1)*1e-5
            p1=log10(1+qi(i,k,j)/1e-5)+result(itile,4)*ys(4)+ym(4)
            if(p1<0) p1=0
            qi_new(i,k,j)=(10**p1-1)*1e-5
            p1=log10(1+qs(i,k,j)/1e-5)+result(itile,5)*ys(5)+ym(5)
            if(p1<0) p1=0
            qs_new(i,k,j)=(10**p1-1)*1e-5
            p1=log10(1+qg(i,k,j)/1e-5)+result(itile,6)*ys(6)+ym(6)
            if(p1<0) p1=0
            qg_new(i,k,j)=(10**p1-1)*1e-5

            th_new(i,k,j)=th(i,k,j)+(result(itile,7)*ys(7)+ym(7))*(1e5/press(i,k,j))**(287.0/1004.0)
            if(qv_new(i,k,j)<0)qv_new(i,k,j)=0
            if(qc_new(i,k,j)<0)qc_new(i,k,j)=0
            if(qr_new(i,k,j)<0)qr_new(i,k,j)=0
            if(qi_new(i,k,j)<0)qi_new(i,k,j)=0
            if(qs_new(i,k,j)<0)qs_new(i,k,j)=0
            if(qg_new(i,k,j)<0)qg_new(i,k,j)=0
            
        enddo   
    enddo
enddo   
end subroutine updatefields
