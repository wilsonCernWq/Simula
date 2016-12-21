subroutine move_test (pos, d, name, rate, sid)
use sim_mod
use mono_mod
implicit none
  integer, intent(in) :: pos(2), d, name, sid
  real(8), intent(out) :: rate
    !=> initialization
    if (test_occupied(pos, d, name, sid)) then
      rate = 0.0
    else
      !=> select monomer
      if (name == TPYP) then
        rate = RDEM
      !=> select monomer
      else if (name == LEAD) then
        rate = RDEP
      end if
    end if
  return
end subroutine move_test