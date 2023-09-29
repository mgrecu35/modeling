

subroutine read_tabulated_mp
  use tabulated_mp
  use module_mp_morr_two_moment2
  integer :: nt1, np1, nf1
  call get_dims(nt1,np1,nf1)
  print*,nt1,np1,nf1
  nt=nt1
  np=np1
  nf=nf1
  if(.not.allocated(tdef)) allocate(tdef(nt1),pdef(np1),fract(nf1))
  if(.not.allocated(qvsat_map)) &
       allocate(qvsat_map(nt1,np1),qvsat_ice_map(nt1,np1),&
       fract_ice_map(nt1,np1,nf1),fract_map(nt1,np1,nf1))
  call read_tables(tdef,pdef,fract,qvsat_map,fract_map,qvsat_ice_map,fract_ice_map)
  print*, tdef(1), tdef(nt1)
  print*, pdef(1), pdef(np1)
  print*, fract(:)
  print*, qvsat_map(1,1), qvsat_map(nt1,np1)
  print*, fract_map(1,1,1), fract_map(nt1,np1,nf1)
  print*, '89,125',qvsat_map(89,125)
  do i=1,nt
     do j=1,np
        evs=polysvp(tdef(i),1)
        qvs = .622*EVS/(pdef(j)-evs)
        !print*, qvs, qvsat_ice_map(i,j), tdef(i), pdef(j)
     enddo
  enddo
        
  !stop
  !print*, fract
end subroutine read_tabulated_mp
