module file_io
  implicit none

  contains
  function file_exist(filename) result(res)
    implicit none
    character(len=*), intent(in) :: filename
    logical                      :: res

    ! Check if the file exist
    inquire( file=trim(filename), exist=res)
  end function

end module
