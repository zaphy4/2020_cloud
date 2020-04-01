module adv
! ---------------------------------------------------------
!
! interface for grid type, constant or stretching
!
! ---------------------------------------------------------
implicit none

contains

! #########################################################

!subroutine adv_1d (temp, zlev, szlev, w, nz, nt)
subroutine adv_1d (temp, zlev, szlev, w, nz, dt)

  implicit none
  integer, intent(in)                    :: nz, dt 
  real, dimension(nz),   intent(inout)   :: temp, zlev
  real, dimension(nz+1), intent(inout)   :: w, szlev
  
  ! time var.
  integer                :: i
  real, dimension(nz+1)  :: z
  real, dimension(nz)    :: flux

  ! boundary condition
  w(nz+1) = 0
  w(1)    = 0

  
  do i = 1, nz-1               !-- i=1 calculates Cj  => i=nz calculates Cnz+1
  flux(i) = temp(i)*w(i+1)*dt  !-- FLUX start between Cj-1 and Cj

  if (i==1) then                !-- temp initial
  temp(1) = temp(1) - flux(1)*1./(szlev(2)-szlev(1))
  end if

  flux(i+1) = temp(i+1)*w(i+2)*dt
  temp(i+1) = temp(i+1) + (flux(i)-flux(i+1))/(szlev(i+2)-szlev(i+1))

  end do


!  flux(1) = temp(1)*w(2)*dt     ! -- staggered grid. FLUX start between Cj-1 and Cj
!  flux(2) = temp(2)*w(3)*dt
!  temp(2) = temp(2) + 1./(szlev(3)-szlev(2)) * (flux(1)-flux(2))  

!  flux(2) = temp(2)*w(3)*dt     ! -- staggered grid. FLUX start between Cj-1 and Cj
!  flux(3) = temp(3)*w(4)*dt
!  temp(3) = temp(3) + 1./(szlev(4)-szlev(3)) * (flux(2)-flux(3))  



  ! boundary condition
  !z(1,:)    = -u(2,:)  ; v(1,:)    =  v(2,:)   ; h(1,:)    = h(2,:)
  !z(nx+1,:) = -u(nx,:) ; v(nx+1,:) =  v(nx,:)  ; h(nx+1,:) = h(nx,:)
  !z(:,1)    =  u(:,2)  ; v(:,1)    = -v(:,2)   ; h(:,1)    = h(:,2)
  !z(:,ny+1) =  u(:,ny) ; v(:,ny+1) = -v(:,ny)  ; h(:,ny+1) = h(:,ny)
  !z(1)    = w(2)
  !z(ny+1) = w(ny)
 

end subroutine adv_1d

! #########################################################

subroutine intp_1d (z,Tz, nz,zlev,temp)
  implicit none
  integer, intent(in)                :: nz
  real, dimension(nz), intent(inout) :: temp, zlev
  real                               :: model, datas

  integer, parameter        :: xx=240, yy=121, zz=42 !-- MERRA2 dims
  real, dimension(zz)       :: Tz, z
  integer                   :: i
  integer, dimension(:), allocatable :: aa

  allocate(aa(zz))

  !-- z: merra2, zlev: intp
  do i = 1, nz
   aa = minloc(abs(z-zlev(i)))
   print*, "model:", zlev(i), "merra2:", z(aa(1)), Tz(aa(1))  !-- check nearest

  if (zlev(i) .eq. z(aa(1))) then
  temp(i) = Tz(aa(1))

  else if (zlev(i).lt.minval(z)) then !-- extrapolation
  temp(i) = Tz(aa(1)) + (zlev(i)-z(aa(1)))*(Tz(aa(1)+1)-Tz(aa(1)))/(z(aa(1)+1)-z(aa(1)))

  else if (zlev(i).gt.maxval(z)) then ! extrapolation, just adjust former gradient 
  temp(i) = Tz(aa(1)) + (zlev(i)-z(aa(1)))*(Tz(aa(1))-Tz(aa(1)-1))/(z(aa(1))-z(aa(1)-1))


  else if (zlev(i).gt.z(aa(1))) then ! interpolation
  temp(i) = Tz(aa(1)) + (zlev(i)-z(aa(1)))*(Tz(aa(1)+1)-Tz(aa(1)))/(z(aa(1)+1)-z(aa(1)))

  else if (zlev(i).lt.z(aa(1))) then ! interpolation
  temp(i) = Tz(aa(1)-1) + (zlev(i)-z(aa(1)-1))*(Tz(aa(1))-Tz(aa(1)-1))/(z(aa(1))-z(aa(1)-1))

  end if

  print*, temp(i) !-- check

  end do


end subroutine

! #########################################################

end module

