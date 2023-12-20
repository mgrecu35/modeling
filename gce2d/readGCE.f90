
subroutine open_for_reading()
  open(14,FILE='data.3min',FORM='UNFORMATTED')
end subroutine open_for_reading

subroutine close_gce_file()
  close(14)
end subroutine close_gce_file

subroutine read_tape2(x,n1,n2)
  integer :: itape, n1,n2
  real,intent(out) ::  x(n1,n2)
  real*4 x4(n1,n2)
  itape=14
  read(itape) x4
  do k=1,n2
     do i=1,n1
        x(i,k)=x4(i,k)
     enddo
  enddo  
end subroutine read_tape2



subroutine wtap1(x,n)
  !c     ******   write 1-d array data to data file   ******
  integer :: n, itape
  
  real,intent(out)::x(n)
  real*4::x4(n)
  itape=14
  read (itape) x4
  do k=1,n
     x(k)=x4(k)
  enddo
end subroutine wtap1
      
