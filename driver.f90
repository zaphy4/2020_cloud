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

  use file_io, only: file_exist

  implicit none

  ! [private variable]
  integer :: it
  real, dimension(:), allocatable :: lev, temp

  ! [namelist input]
  integer          :: nz,   & ! The number of levels  
                      nt      ! The number of time
  real             :: ps,   & ! Surface pressure 
                      ptop    ! model top
  character(len=24) :: stat    ! ideal or real case

  namelist /variable/nz, nt, ps, ptop, stat


  ! [read namelist]
  if (file_exist('namelist')) then
     open(unit=101, file="namelist", status="old")
     read(101, nml=variable) 
  else
     print*, "Namelist file doesn't exist"
     call exit ! program end
  end if

  ! [allocatable]
  if (allocated(lev))  deallocate(lev)
  if (allocated(temp)) deallocate(temp)
  allocate(lev(nz))
  allocate(temp(nz))

  ! [initialization]
  if (stat=='ideal') then
     print*, 'make subroutine of ideal case'
     call temp_ideal(nz,temp)
  else if (stat=='real') then
     print*, 'make subroutine of real case'
  end if

  ! [time integral]
  do it = 1, nt
  end do
  write(*,*) temp

  if (allocated(lev))  deallocate(lev)
  if (allocated(temp)) deallocate(temp)

end program driver
