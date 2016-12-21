module func_xml_var

  use global
  implicit none
  private

  type, public :: var
     private
     character(len=1), allocatable :: data(:)
     character(30), public :: name
     character(10), public :: type
     integer, public :: length
   contains
     procedure :: encode_generic
     procedure :: decode_integer
     procedure :: decode_real
     procedure :: print
     generic   :: encode => encode_generic
     generic   :: decode => decode_integer, decode_real
  end type var

contains

  subroutine encode_generic(self, data, name)
    class(*)     :: data(:)
    class(var)   :: self
    character(*) :: name
    integer      :: length
    ! prepare
    self % name = name
    self % length = size(data)
    ! encode
    select type (data)
    type is (integer)
       self % type = "integer"
       if (allocated(self % data)) deallocate(self % data)    
       length = size(transfer(data, self % data))
       allocate(self % data(length))
       self % data = transfer(data, self % data)
    type is (real(dp))
       self % type = "real(dp)"
       if (allocated(self % data)) deallocate(self % data)    
       length = size(transfer(data, self % data))
       allocate(self % data(length))
       self % data = transfer(data, self % data)
    class default
       stop "ERROR unknow data type"
    end select
    return
  end subroutine encode_generic

  subroutine decode_integer(self, data)
    integer, allocatable :: data(:)
    class(var) :: self
    ! safty check
    if (self % type /= "integer") stop "ERROR (DECODE) data type mismatch"
    ! clean memory
    if (allocated(data)) deallocate(data)
    ! decode
    allocate(data(self % length))
    data = transfer(self % data, data)
    return
  end subroutine decode_integer

  subroutine decode_real(self, data)
    real(dp), allocatable :: data(:)
    class(var) :: self
    ! safty check
    if (self % type /= "real(dp)") stop "ERROR (DECODE) data type mismatch"
    ! clean memory
    if (allocated(data)) deallocate(data)
    ! decode
    allocate(data(self % length))
    data = transfer(self % data, data)
    return
  end subroutine decode_real

  subroutine print(self)
    class(var) :: self
    integer , allocatable :: vI(:)
    real(dp), allocatable :: vF(:)
    if (self % type == "integer") then
       call self % decode(vI)
       print *, "var (INTEGER) ", trim(self % name), vI
       deallocate(vI)
    elseif (self % type == "real(dp)") then
       call self % decode(vF)
       print *, "var (REAL) ", trim(self % name), vF
       deallocate(vF)
    end if
    return
  end subroutine print
  
end module func_xml_var
