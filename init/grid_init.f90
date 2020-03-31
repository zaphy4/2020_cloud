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

subroutine grid_const (ps,ptop,nz,plev,zlev,splev,szlev)
  implicit none
  integer, intent(in)             :: nz
  real, intent(in)                :: ps,ptop
  real, dimension(nz), intent(inout)   ::  plev,  zlev
  real, dimension(nz+1), intent(inout) :: splev, szlev

  ! private variable
  integer :: i
  real    :: dz

  dz = (ptop-ps)/nz

  splev(1) = ps
  do i = 1, nz
    splev(i+1) = dz + splev(i) 
  end do    
  !szlev = 

  !print*, ps, splev(1), splev(nz+1)
  !print*, splev

end subroutine grid_const

! #########################################################

subroutine grid_stret (ps,ptop,nz,plev,zlev,splev,szlev)
  implicit none
  integer, intent(in)                :: nz
  real, intent(in)                   :: ps,ptop
  real, dimension(nz), intent(inout)   :: plev,zlev
  real, dimension(nz+1), intent(inout) :: splev, szlev

end subroutine grid_stret


! #########################################################




end module grid_init
