!------------------------------------------------------------------------------
! Department of Physics, The Hong Kong University of Science and Technology
!------------------------------------------------------------------------------
!
! Module: rcd_mod
!
!> @author
!> Wilson
!
! DESCRIPTION: 
!> Construct rcd list
!
!------------------------------------------------------------------------------
module rcd_mod
use glob_mod
use mono_mod
implicit none

integer, allocatable, save :: rcd(:,:)
integer, parameter :: dim = 4 !< x, y, d, name 
integer, save :: wid = 0

contains

subroutine rcd_init()
  integer :: ntot, i
  !=> compute total number of monomers
  wid = max_sur + dim
  ntot = 0
  do i = 1, no_of_comp
    ntot = ntot + comps(i)%num
  end do
  !=> initialize rcd list
  if (allocated(rcd)) deallocate(rcd)
  allocate(rcd(wid, ntot))
  rcd = 0
  return
end subroutine rcd_init

subroutine rcd_set_pos(id, name, pos, dir)
  integer :: id, name, pos(2), dir
    rcd(1, id) = modulo(pos(1)-1, SIDE)+1 !< position
    rcd(2, id) = modulo(pos(2)-1, SIDE)+1 !< position
    rcd(3, id) = dir !< direction
    rcd(4, id) = name !< type name
  return
end subroutine rcd_set_pos

elemental function rcd_name (id) result (name)
  integer, intent(in) :: id
  integer :: name
    if (id == 0) then
      name = 0
    else 
      name = rcd(4, id)
    end if
  return
end function rcd_name

end module rcd_mod