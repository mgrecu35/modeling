
  
subroutine init_nav(ix,jy,t_end_py,ui_py)
  

  ! by D. Orchard (2012) based on the classic code from:
  !
  ! Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
  ! Numerical Simulation in Fluid Dynamics,
  ! SIAM, 1998.
  !
  ! http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

  use state_vars
  integer, intent(out) :: ix, jy
  real, intent(out):: t_end_py
  real, intent(in) :: ui_py
  u = ui_py
  v = vi
  p = 0.0

  ix=imax
  jy=jmax
  t_end_py=t_end
  ! init flags
  call init_flag(delx, dely, flag, ibound)

  call apply_boundary_conditions(u, v, flag, 0.0)
  
end subroutine INIT_NAV

subroutine py_time_step(t_py,iters_py)
  use state_vars
  real, intent(out):: t_py
  integer, intent(out) :: iters_py
  !do while (t < t_end)
  del_t = set_timestep_interval(u, v, del_t)
  ifluid = (imax * jmax) - ibound

  call compute_tentative_velocity(u, v, f, g, flag, del_t)
  
  call compute_rhs(f, g, rhs, flag, del_t)
  
  if (ifluid > 0) then
     itersor = poisson(p, rhs, flag, res, ifluid)
  else
     itersor = 0
  end if
  
  print '(i0.1, " t:", f0.7, ", del_t:", f0.7, ", SOR iters:", i3.3, ", res:", f0.7, ", &
       bcells:", i0)', iters, t+del_t, del_t, itersor, res, ibound
  
  !("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
  !                iters, t+del_t, del_t, itersor, res, ibound);
  
  call update_velocity(u, v, f, g ,p, flag, del_t)
  
  call apply_boundary_conditions(u, v, flag, t)
  
  ! An alternate way to print output using the simulation time
  !if (toLogical(outputFlag) .and. (t >= t_count)) then
  !    call write_ppm(u, v, p, flag, output_counter, outname)
  !    output_counter = output_counter + 1
  !    t_count = t_count + output_freqt
  !end if
  
  if (toLogical(outputFlag) .and. (mod(iters, output_freq) == 0)) then
     call write_ppm(u, v, p, flag, iters, outname)
     output_counter = output_counter + 1
     t_count = t_count + output_freqt
  end if
  
  t = t + del_t
  iters = iters + 1
  !end do
  iters_py=iters
  t_py=t
end subroutine py_time_step

subroutine get_uvp(u_py,v_py,p_py,zeta_py,imax_py,jmax_py)
  use state_vars
  integer :: imax_py, jmax_py
  real, intent(out) ::  u_py(0:imax_py+1, 0:jmax_py+1), v_py(0:imax_py+1, 0:jmax_py+1), &
       p_py(0:imax_py+1, 0:jmax_py+1)
  real, intent(out) :: zeta_py(0:imax_py+1, 0:jmax_py+1)
  u_py=u
  v_py=v
  p_py=p
  call calc_zeta(u, v, flag, zeta_py)
end subroutine get_uvp


