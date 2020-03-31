!-----------------------------------------------------
!
! 2020 Cloud Modeling Class
! 1-D advection model 
! 
! Key structure : namelist 
!                 stretching 
!                 staggered grid
!                 finite volume method
!
! 20.3.27 : the simplest form of the namelist, subroutine, module works well.
!
!-----------------------------------------------------
program driver

  use netcdf
  use   file_io, only: file_exist
  use  var_init, only: temp_init, w_init, temp_ideal, w_ideal 
  use grid_init, only: grid_const, grid_stret

  implicit none
  ! [namelist input]
  integer          :: nz,   & ! The number of levels  
                      nt      ! The number of time

  real             :: zsfc,   & ! Surface pressure 
                      ztop      ! model top


  character(len=24) :: stat    ! al or real case
  character(len=24) :: grid    ! const or stret
  namelist /variable/nz, nt, zsfc, ztop, stat, grid



  ! [private variable]
  integer :: it
  real, dimension(:), allocatable :: q, w, temp, time, lev, zlev, szlev 
                                   !-- q [kg/kg], w[m/s], temp[K]


  ! [read namelist]
  if (file_exist('namelist')) then
     
  open(unit=101, file="namelist", status="old")
  read(101, nml=variable) 
  allocate(time(nz))
  allocate(temp(nz))
  allocate(  q(nz))
  allocate(zlev(nz))
  allocate(szlev(nz+1))
  allocate(    w(nz+1))
 
  else

  print*, "Namelist file doesn't exist"
  call exit ! program end

  end if



  ! [initialization]
  ! Make level first
  if (grid=='const') then
    print*, 'grid: constant'
    call grid_const(zsfc, ztop, nz, zlev, szlev)
  else if (grid=='stret') then
    print*, 'grid: stretching (NOT YET)'
    call grid_stret(zsfc, ztop, nz, zlev, szlev)
  end if

  ! read NetCDF file and convert from p->z and interpolation (temp, qv)
  if (stat=='ideal') then
    print*, 'variable: ideal case (NOT YET)'
    call exit                     !-- exit program
    call temp_ideal(nz,zlev,temp) !-- not staggered
    call    w_ideal(nz,szlev,w)   !-- staggered grid
  else if (stat=='real') then
    print*, 'variable: real case'
    call temp_init(nz,zlev,temp,q)
    call    w_init(nz,szlev,w)
  end if



  ! [time integral]
  do it = 1, nt
  end do
  

end program driver
