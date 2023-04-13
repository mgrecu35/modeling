module state_vars
  use init_mod
  use helpers_mod
  use boundary_mod
  use simulation_mod
  use output_mod
  !implicit none
  character(10) :: outname = "output"
  

  real :: t = 0, res = 0, del_t = 0.003, t_count = 0.0
  integer :: i, j, itersor = 0, ifluid = 0, ibound = 0, init_case, iters = 0
  integer :: output_counter = 0

  real u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1), p(0:imax+1, 0:jmax+1)
  real rhs(0:imax+1, 0:jmax+1), f(0:imax+1, 0:jmax+1), g(0:imax+1, 0:jmax+1)
  integer flag(0:imax+1, 0:jmax+1)

end module state_vars
