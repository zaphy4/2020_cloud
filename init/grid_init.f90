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
  integer, intent(in)                  :: nz
  real, intent(in)                     :: zsfc,ztop
  real, dimension(nz), intent(inout)   :: zlev
  real, dimension(nz+1), intent(inout) :: szlev

  ! private variable
  real, dimension(nz) :: array, z_weight, ratio
  integer :: i
  real    :: dz

  do i = 1, nz  
    array(i) = i
  end do

  z_weight = exp(-40./nz*(array-0)**2/(nz*1.5)) + exp(-40./nz*(array-(2*nz/3.))**2/(nz*1.5)) +0.1
  
  z_weight   = int( (1./z_weight)*100)/10 
  ratio      = real(z_weight)/sum(z_weight) 
  
  dz = ztop - zsfc

  szlev(1) = zsfc
  do i = 1, nz-1
    szlev(i+1) = szlev(i) + dz*ratio(i)
    zlev(i)    = (szlev(i)+szlev(i+1))/2.
  end do

  szlev(nz+1) = ztop
  zlev(nz)    = (szlev(nz)+szlev(nz+1))/2. 

  print*, szlev
  print*, zlev

end subroutine grid_stret


! #########################################################




end module grid_init
