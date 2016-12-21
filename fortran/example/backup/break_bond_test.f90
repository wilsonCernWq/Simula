subroutine break_bond_test (sid, tid, pos, rate)
use sim_mod
use rcd_mod
implicit none
  integer, intent(in) :: sid, tid, pos
  real(8), intent(out) :: rate
  integer :: name
    !=> initialization
    rate = 0.0
    !=> self information
    name = rcd_name(sid)
    !=> select monomer
    if (name == TPYP) then
      !=> hydrogen bond
      if (pos <= 4) then 
        rate = RBHY
      !=> coordination bond
      else if (pos > 4) then
        if (rcd(pos, tid)==0) then 
          rate = RBCA
        else
          rate = RBCB
        end if
      end if
    else if (name == LEAD) then 
      !=> metal bond
      if (pos <= 4) then
        rate = RBME
      !=> coordination bond
      else if (pos > 4) then
        if (rcd(modulo(pos-3,4)+5, sid)==0) then !=> target_pos
          rate = RBCA
        else
          rate = RBCB
        end if
      end if
    end if
  return
end subroutine break_bond_test