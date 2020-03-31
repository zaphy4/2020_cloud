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

  use   file_io, only: file_exist
  use  var_init, only: temp_init, omega_init
  use grid_init, only: grid_const, grid_stret

  implicit none
  ! [namelist input]
  integer          :: nz,   & ! The number of levels  
                      nt      ! The number of time

  real             :: ps,   & ! Surface pressure 
                      ptop    ! model top


  character(len=24) :: stat    ! al or real case
  character(len=24) :: grid    ! const or stret
  namelist /variable/nz, nt, ps, ptop, stat, grid



  ! [private variable]
  integer :: it
  real, dimension(:), allocatable :: temp, time, lev, plev, zlev, splev, szlev



  ! [read namelist]
  if (file_exist('namelist')) then
     
  open(unit=101, file="namelist", status="old")
  read(101, nml=variable) 
  allocate(time(nz))
  allocate(temp(nz))
  allocate(zlev(nz))
  allocate(plev(nz))
  allocate(szlev(nz+1))
  allocate(splev(nz+1))
 

  else

  print*, "Namelist file doesn't exist"
  call exit ! program end

  end if



  ! [initialization]
  

  if (grid=='const') then
    print*, 'grid: constant'
    call grid_const(ps, ptop, nz, plev, zlev, splev, szlev)
  else if (grid=='stret') then
    print*, 'grid: stretching'
    call grid_stret(ps, ptop, nz, plev, zlev, splev, szlev)
  end if



  if (stat=='ideal') then
    print*, 'variable: ideal case'
    call temp_init(nz,temp)
  else if (stat=='real') then
    print*, 'variable: real case'
  end if


  ! [time integral]
  do it = 1, nt
  end do
  write(*,*) temp
  

end program driver
