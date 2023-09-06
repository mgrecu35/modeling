! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to hugh morrison's microphysics code 
!
module mphys_morr_two_moment

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, mom_names, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use module_mp_morr_two_moment

  Implicit None

  !Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains
  
  Subroutine mphys_morrison_interface
    
   


  end Subroutine mphys_morrison_interface

end module mphys_morr_two_moment
