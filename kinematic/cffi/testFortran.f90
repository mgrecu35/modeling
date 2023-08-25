program testCFFI
  real a(2,2)
  a(1,:)=(/1,2/)
  a(2,:)=(/3,4/)
  print*, a(1,*)
  call hello_world(a,2,2)
end program testCFFI
