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
  use       adv, only: adv_1d
  use   file_io, only: file_exist
  use  var_init, only: temp_init, w_init, temp_ideal
  use grid_init, only: grid_const, grid_stret

  implicit none
  ! [namelist input]
  integer            :: nz,   & ! The number of levels  
                        nt,   & ! The number of time
                        dt      ! time interval

  integer            :: wflag   ! 1. constant (1m/s) 
                                ! 2. sin(-90) - sin(90)
                                ! 3. sin(0)^2 - sin(360)^2
                                ! 4. sfc and mid-level has strong wind

  real               :: zsfc,   & ! Surface pressure 
                        ztop      ! model top


  character(len=24) :: stat    ! al or real case
  character(len=24) :: grid    ! const or stret
  namelist /variable/nz, nt, zsfc, ztop, stat, grid, dt, wflag



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
    print*, 'grid: stretching '
    call grid_stret(zsfc, ztop, nz, zlev, szlev)
  end if

  ! read NetCDF file and convert from p->z and interpolation (temp, qv)
  if (stat=='ideal') then
    print*, 'variable: ideal case !!!!!'
    call temp_ideal(nz,zlev,temp,q)
    call    w_init(nz,szlev,w, wflag)
  else if (stat=='real') then
    print*, 'variable: real case'
    call temp_init(nz,zlev,temp,q)
    call    w_init(nz,szlev,w, wflag)
  end if



  ! [time integral]
write(*,'(a,i4,a,20f10.4,a,1f10.4)') "temp ",1, ":", q*1000, "  sum:", sum(q*1000)
  do it = 2, nt
    !call adv_1d(temp, zlev, szlev, w, nz,dt)
   !write(*,'(a,i4,a,20f10.4,a,1f10.4)') "temp ",it, ":", temp, "  sum:", sum(temp)
    call adv_1d(  q, zlev, szlev, w, nz,dt)     
   write(*,'(a,i4,a,20f10.4,a,1f10.4)') "temp ",it, ":", q*1000, "  sum:", sum(q*1000)
  end do
  

end program driver
