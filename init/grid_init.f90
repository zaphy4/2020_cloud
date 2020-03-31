module grid_init
! ---------------------------------------------------------
!
! interface for grid type, constant or stretching
!
! ---------------------------------------------------------
implicit none
!use .. ~~

contains

! #########################################################

subroutine grid_const (zsfc,ztop,nz,zlev,szlev)
  implicit none
  integer, intent(in)             :: nz
  real, intent(in)                :: zsfc,ztop
  real, dimension(nz), intent(inout)   :: zlev
  real, dimension(nz+1), intent(inout) :: szlev

  ! private variable
  integer :: i
  real    :: dz

  dz = (ztop-zsfc)/nz

  szlev(1) = zsfc
  do i = 1, nz
    szlev(i+1) = dz + szlev(i) 
    zlev(i)    = (szlev(i)+szlev(i+1))/2.
  end do    

  !print*, szlev  !-- check
  !print*, zlev


end subroutine grid_const

! #########################################################

subroutine grid_stret (zsfc,ztop,nz,zlev,szlev)
  implicit none
  integer, intent(in)                :: nz
  real, intent(in)                   :: zsfc,ztop
  real, dimension(nz), intent(inout)   :: zlev
  real, dimension(nz+1), intent(inout) :: szlev

end subroutine grid_stret


! #########################################################




end module grid_init
