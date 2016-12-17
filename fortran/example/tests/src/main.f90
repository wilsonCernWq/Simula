module test
  type :: aaaa
     integer :: value
  end type aaaa
contains
  function get(a) result (r)
    type(aaaa), target :: a
    type(aaaa), pointer :: r
    r => a
    return
  end function get

end module test

program main
  use test
  implicit none
  type(aaaa) :: ssss, d
  integer, pointer :: a, b
  integer :: c

  allocate(a)
  allocate(b)
  a = 10
  b = a
  b = 20

  print *, a, b


  stop
end program main
