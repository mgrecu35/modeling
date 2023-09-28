module tabulated_mp
    integer :: nt, np, nf
    real,allocatable :: tdef(:), pdef(:), fract(:)
    real,allocatable :: qvsat_map(:,:)
    real,allocatable :: fract_map(:,:,:)
end module tabulated_mp

program test
    use tabulated_mp
    integer :: nt1, np1, nf1
    call get_dims(nt1,np1,nf1)
    print*,nt1,np1,nf1
    allocate(tdef(nt1),pdef(np1),fract(nf1))
    allocate(qvsat_map(nt1,np1),fract_map(nt1,np1,nf1))
    call read_tables(tdef,pdef,fract,qvsat_map,fract_map)
    print*, tdef(1), tdef(nt1)
    print*, pdef(1), pdef(np1)
    print*, fract(:)
    print*, qvsat_map(1,1), qvsat_map(nt1,np1)
    print*, fract_map(1,1,1), fract_map(nt1,np1,nf1)
    !print*, fract
end program test