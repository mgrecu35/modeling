
module keras_def
  use mod_kinds, only: ik, rk
  use mod_network, only: network_type

  implicit none

  type(network_type) :: net

  real(rk), allocatable :: result(:), input(:)
  !real(rk), allocatable :: result2(:), input2(:)

end module keras_def
