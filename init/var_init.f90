module var_init
! ---------------------------------------------------------
!
! interface for grid type, constant or stretching
!
! ---------------------------------------------------------
  use netcdf
implicit none


contains

! #########################################################

subroutine temp_init (nz,zlev,temp,q)
  use netcdf

  implicit none
  integer, intent(in)                :: nz
  real, dimension(nz), intent(inout) :: q,temp, zlev

  integer, parameter        :: xx=240, yy=121, zz=42 !-- MERRA2 dims
  real, dimension(xx,yy,zz) :: T, QV, TV
  real, dimension(zz)       :: lev, Tz, QVz, TVz, z
  real, dimension(yy)       :: lat
  real, dimension(xx)       :: lon
  integer                   :: ncid, varid, k
  integer, parameter        :: llat=62 , llon=107
  real, parameter           :: Rd=287., g=9.81 
  real                      :: dz
  character(len=200) :: INFILE
  

  write(INFILE,'(a)') '../data/MERRA2_400.inst3_3d_asm_Np.201708.SUB.nc4'

  call check( nf90_open(trim(INFILE), NF90_NOWRITE, ncid) ) !-- NF90_NOWRITE :
  call check( nf90_inq_varid(ncid, "T", varid) )          !-- Get the varID of
  call check( nf90_get_var(ncid, varid, T) )              !-- Read the data

  call check( nf90_inq_varid(ncid, "QV", varid) )    !-- specific humidity
  call check( nf90_get_var(ncid, varid, QV) )        !-- kg/kg

  call check( nf90_inq_varid(ncid, "lev", varid) )    
  call check( nf90_get_var(ncid, varid, lev) )      
 
  call check( nf90_inq_varid(ncid, "lat", varid) )    
  call check( nf90_get_var(ncid, varid, lat) )      

  call check( nf90_inq_varid(ncid, "lon", varid) )    
  call check( nf90_get_var(ncid, varid, lon) )      
  call check( nf90_close(ncid) )

  print*, "lat:", lat(llat), "lon:", lon(llon)
  !print*,T(llon,llat,:) ! -- check
  Tz  = T(llon,llat,:)
  QVz = QV(llon,llat,:)



  ! convert from p->z
  TVz = Tz*(1+0.61*QVz) !-- specfic humidity is almost same as mixing ratio
  z(1) = 0.             !-- Assume, 1000 hPa = 0 m

  do k = 2, zz          
    dz = -1*(Rd/g)*0.5*(TVz(k)+TVz(k-1))*(log(lev(k)*100) - log(lev(k-1)*100))
    z(k) = z(k-1) + dz

  !  print*, z(k), lev(k) !-- check
  end do


  ! interpolate T and QV to model z-level
  print*, "temp interpolate .."
  call intp_1d(z, Tz, nz,zlev,temp)
  print*, "q interpolate .."
  call intp_1d(z,QVz, nz,zlev,   q)

end subroutine

! #########################################################

subroutine temp_ideal (nz,zlev,temp,q)
  implicit none
  integer, intent(in)                :: nz
  real, dimension(nz), intent(inout) :: q,temp, zlev

  temp = 0
  temp(5) = 10

  q    = 0
  q(5) = 0.01

end subroutine

! #########################################################

subroutine w_init (nz, szlev, w)
  implicit none
  integer, intent(in)                  :: nz
  real, dimension(nz+1), intent(inout) :: w, szlev
  integer                              :: i

  do i = 1, nz+1
  w(i) = i
  end do

  ! method 1: constant
  w = 1 !-- m/s

  ! method 2: sin(-90) ~ sin(90)
  !w = (w-(nz+1)/2.)/real(nz+1) *3.14
  !w = sin(w)
  !print*, w

  ! method 3: sin^2
  !w = w/real(nz+1) * 3.14 * 2 
  !w = sin(w)**2
  !print*, w


end subroutine

! #########################################################

subroutine w_ideal (nz,szlev, w)
  implicit none
  integer, intent(in)                  :: nz
  real, dimension(nz+1), intent(inout) :: w, szlev




end subroutine

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

  subroutine check(status)
    implicit none
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

! #########################################################


end module

