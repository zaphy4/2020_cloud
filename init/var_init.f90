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

subroutine w_init (nz, szlev, w, wflag)
  implicit none
  integer, intent(in)                  :: nz, wflag
  real, dimension(nz+1), intent(inout) :: w, szlev
  integer                              :: i

  do i = 1, nz+1
  w(i) = i
  end do

  ! method 1: constant
  if (wflag.eq.1) then
  w = 1 !-- m/s
  end if

  ! method 2: sin(-90) ~ sin(90)
  if (wflag.eq.2) then
  w = (w-(nz+1)/2.)/real(nz+1) *3.14
  w = sin(w)
  end if

  ! method 3: sin^2
  if (wflag.eq.3) then
  w = w/real(nz+1) * 3.14 * 2 
  w = sin(w)**2
  end if

  ! method 4: surface and cloud have strong w-wind
  if (wflag.eq.4) then
  w = exp(-1*(w-0)**2/(nz*5)) + exp(-1*(w-(nz/2))**2/(nz*10)) 
  end if

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

SUBROUTINE initial_sounding(nz, zlev, z, t, th, qv, p)

  IMPLICIT NONE

  INTEGER             :: nz
  REAL, DIMENSION(nz) :: zlev
  REAL, DIMENSION(nz) :: z, t, th, qv, p

  REAL,    PARAMETER :: kappa  = 0.288      !! check
  REAL,    PARAMETER :: p1000  = 1008.7
  REAL,    PARAMETER :: grav   = 9.81
  REAL,    PARAMETER :: rd     = 287.
  INTEGER, PARAMETER :: max_nl = 10000

  REAL :: pp, zz, tt, tdd, kt, pa
  REAL :: zl1, wgt, pp1

  INTEGER :: readstus, lnum, k, kk, ik

  REAL, DIMENSION(max_nl) :: z_in, t_in, th_in, qv_in, p_in

  CHARACTER(LEN=150) :: filename

  filename = '../data/OBS_SONDE.csv'

  OPEN(10, FILE=TRIM(filename))

  lnum = 0
  DO WHILE(.TRUE.)
     READ(10,*,IOSTAT=readstus)
     IF (readstus /= 0) EXIT
     lnum = lnum + 1
  ENDDO
!print*, lnum

  REWIND(10)

  DO k = lnum, 1, -1
!print *, k

     READ(10,*) pp, zz, tt, tdd

     kt = tt + 273.15
     pa = pp*100 ! hPa -> Pa

     p_in(k) = pa ! Pa
     p_in(k) = pp**(1./kappa)
     t_in(k) = kt ! K

     IF (int(zz) == 0) THEN
        z_in(k) = 44330.*(1.-(pp/p1000)**(1./5.255))
     ELSE
        z_in(k) = zz ! gpm
     ENDIF

     th_in(k) = (kt)*(p1000/pp)**kappa

     pp1 = exp(20.386 - (5133./kt))
     qv_in(k) = 0.622*(pp1/(pp-(0.378*pp1)))  ! kg/kg
     qv_in(k) = 0.622*(pp1/pp)

  ENDDO

  DO k = 1, nz

     zl1 = zlev(k)

     DO kk = 1, lnum
        IF (zl1 >= z_in(kk)) ik = kk
     ENDDO

     IF (ik == lnum) THEN
        th(k) = th_in(ik)
        qv(k) = qv_in(ik)
     ELSE
        wgt = (zl1 - z_in(ik)) / (z_in(ik+1) - z_in(ik))
        th(k) = (1.-wgt)*th_in(ik) + wgt*th_in(ik+1)
        qv(k) = (1.-wgt)*qv_in(ik) + wgt*qv_in(ik+1)
     ENDIF

     p(k) = p_in(ik)**kappa                 &
           - kappa * p1000**kappa * grav    &
           / rd / ((th(k) + th_in(ik))/2.0) &
           * (zl1 - z_in(ik))

     p(k) = p(k)**(1./kappa)

  ENDDO

  t = th*(p/p1000)**kappa

  CLOSE(10)

END SUBROUTINE initial_sounding

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

