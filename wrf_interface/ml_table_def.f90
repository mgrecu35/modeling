module tabulated_mp
    integer :: nt, np, nf
    real,allocatable :: tdef(:), pdef(:), fract(:)
    real,allocatable :: qvsat_map(:,:)
    real,allocatable :: fract_map(:,:,:)
    real,allocatable :: qvsat_ice_map(:,:)
    real,allocatable :: fract_ice_map(:,:,:)
end module tabulated_mp
