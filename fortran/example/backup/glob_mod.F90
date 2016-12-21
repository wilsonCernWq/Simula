module glob_mod
implicit none

!% constants
integer, parameter :: SIDE = 30
integer, parameter :: SYMM = 1
integer, parameter :: HOPNUM = 1E9
integer, parameter :: HOPSTP = 1E6

integer, parameter :: OUTNUM = 100
integer, parameter :: OUTCNT = HOPNUM / OUTNUM
integer, parameter :: TPYP = 1, LEAD = 2
real(8), parameter :: TEMPERATURE = 300, V0 = 1E10, KB = 8.617E-5
real(8), parameter :: COBH = 0.2, COBF = 0.4, EME = 0.0, HYB = 0.0, DEM = 0.0, DEP = 0.0

!% substrate
integer, dimension (SIDE, SIDE), save :: sub
contains

subroutine sub_init()
  sub = 0
  return
end subroutine sub_init

elemental function get_sub(x, y) result (v)
  integer, intent(in) :: x, y
  integer :: v
    v = sub(modulo(x-1,SIDE)+1, modulo(y-1,SIDE)+1)
  return
end function get_sub

subroutine set_sub(x, y, v)
  integer, intent(in) :: x, y, v
    !% increase the value of (x, y) by v
    sub(modulo(x-1,SIDE)+1, modulo(y-1,SIDE)+1) = sub(modulo(x-1,SIDE)+1, modulo(y-1,SIDE)+1) + v
  return
end subroutine set_sub

end module glob_mod