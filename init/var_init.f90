module var_init
! ---------------------------------------------------------
!
! interface for grid type, constant or stretching
!
! ---------------------------------------------------------
implicit none



contains

! #########################################################

subroutine temp_init (nz,temp)

  implicit none
  integer, intent(in)              :: nz
  real, dimension(nz), intent(inout) :: temp

  
  temp = 3.5  


end subroutine temp_init

! #########################################################


! #########################################################

subroutine omega_init (nz,temp)
  implicit none
  integer, intent(in)              :: nz
  real, dimension(nz), intent(inout) :: temp

end subroutine omega_init

! #########################################################

end module var_init
