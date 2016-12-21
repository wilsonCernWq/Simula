subroutine form_bond_test (sid, pos, sur, rate, tid)
use glob_mod
use sim_mod
use rcd_mod
implicit none
  integer, intent(in) :: sid, pos, sur(2)
  integer, intent(out) :: tid
  real(8), intent(out) :: rate
  integer :: self_cor(2), target_cor(2)
  integer :: self_name, target_name, dis
  integer :: check_pos(8), check_result(8)

    !% initialization
    check_result = 0
    check_pos = 0
    rate = 0.0

    !% self information
    self_name = rcd_name(sid)
    self_cor = rcd(1:2,sid)
    !% target information
    tid = get_sub(sur(1), sur(2)) !% get target index
    if (tid /= 0) then
      target_name = rcd_name(tid) !% target name
      target_cor = rcd(1:2,tid) !% target coordinates
      dis = (target_cor(2)-self_cor(2))**2 + (target_cor(1)-self_cor(1))**2 !% distance
    else 
      rate = 0.0
      return !% (---> short circuit)
    end if

    !% select monomer
    !=>TPYP (5,6,7,8)--> hydrogen; (9,10,11,12)--> coordination
    !=>LEAD (5,6,7,8)--> metal; (9,10,11,12)--> coordination
    if (self_name == TPYP) then
      !% hydrogen bond
      if (pos <= 4) then !=> pos: 1~4
        if (target_name == TPYP .and. dis == 9) then
          !% get check position
          check_pos(1) = 8 + pos              !% self --> TPYP --> check no coordination bond
          check_pos(2) = 9 + modulo(pos-2,4)  !% self --> TPYP --> check no coordination bond
          check_pos(3) = 9 + modulo(pos,4)    !% target --> TPYP --> check no coordination bond
          check_pos(4) = 9 + modulo(pos+1,4)  !% target --> TPYP --> check no coordination bond
          !% check result
          check_result(1) = rcd(check_pos(1), sid)
          check_result(2) = rcd(check_pos(2), sid)
          check_result(3) = rcd(check_pos(3), tid)
          check_result(4) = rcd(check_pos(4), tid)
          if (sum(check_result) == 0) rate = RHYB
        end if
      !% coordination bond
      else if (pos > 4) then !=> pos: 5~8
        if (target_name == LEAD .and. dis == 8) then
          !% get check position
          check_pos(1) = pos                !% self  --> TPYP --> check no hydrogen bond
          check_pos(2) = 5 + modulo(pos,4)  !% self --> TPYP --> check no hydrogen bond
          check_pos(3) = pos                  !% target --> LEAD --> check no metal bond
          check_pos(4) = 5 + modulo(pos,4)    !% target --> LEAD --> check no metal bond 
          check_pos(5) = 5 + modulo(pos-2,4)  !% target --> LEAD --> check no metal bond
          check_pos(6) = 5 + modulo(pos-3,4)  !% target --> LEAD --> check no metal bond
          check_pos(7) = 9 + modulo(pos,4)    !% target --> LEAD --> check no coordination bond
          check_pos(8) = 9 + modulo(pos-2,4)  !% target --> LEAD --> check no coordination bond
          !% check result
          check_result(1) = rcd(check_pos(1), sid)
          check_result(2) = rcd(check_pos(2), sid)
          check_result(3) = rcd(check_pos(3), tid)
          check_result(4) = rcd(check_pos(4), tid)
          check_result(5) = rcd(check_pos(5), tid)
          check_result(6) = rcd(check_pos(6), tid)
          check_result(7) = rcd(check_pos(7), tid)
          check_result(8) = rcd(check_pos(8), tid)
          if (sum(check_result) == 0) then
            if (rcd(pos+4, tid)==0) then !% target_pos
              rate = RCOA
            else
              rate = RCOC
            end if
          end if
        end if
      end if
    !% select monomer
    else if (self_name == LEAD) then
      !% metal bond
      if (pos <= 4) then !=> pos: 1~4
        if (target_name == LEAD .and. dis == 1) then
          !% get check position
          check_pos(1) = 8 + pos              !% self --> LEAD --> check no coordiantion
          check_pos(2) = 9 + modulo(pos,4)    !% self --> LEAD --> check no coordiantion
          check_pos(3) = 9 + modulo(pos-2,4)  !% self --> LEAD --> check no coordiantion
          check_pos(4) = 9 + modulo(pos-3,4)  !% self --> LEAD --> check no coordiantion
          check_pos(5) = 8 + pos              !% target --> LEAD --> check no coordiantion
          check_pos(6) = 9 + modulo(pos,4)    !% target --> LEAD --> check no coordiantion
          check_pos(7) = 9 + modulo(pos-2,4)  !% target --> LEAD --> check no coordiantion
          check_pos(8) = 9 + modulo(pos-3,4)  !% target --> LEAD --> check no coordiantion
          !% check result
          check_result(1) = rcd(check_pos(1), sid)
          check_result(2) = rcd(check_pos(2), sid)
          check_result(3) = rcd(check_pos(3), sid)
          check_result(4) = rcd(check_pos(4), sid)
          check_result(5) = rcd(check_pos(5), tid)
          check_result(6) = rcd(check_pos(6), tid)
          check_result(7) = rcd(check_pos(7), tid)
          check_result(8) = rcd(check_pos(8), tid)
          if (sum(check_result) == 0) rate = REME
        end if
      !% coordination bond
      else if (pos > 4) then !=> pos: 5~8
        if (target_name == TPYP .and. dis == 8) then
          !% get check position
          check_pos(1) = pos                  !% self --> LEAD --> check no metal bond
          check_pos(2) = 5 + modulo(pos,4)    !% self --> LEAD --> check no metal bond 
          check_pos(3) = 5 + modulo(pos-2,4)  !% self --> LEAD --> check no metal bond
          check_pos(4) = 5 + modulo(pos-3,4)  !% self --> LEAD --> check no metal bond
          check_pos(5) = 4 + pos              !% self --> LEAD --> check no coordination bond
          check_pos(6) = 9 + modulo(pos-2,4)  !% self --> LEAD --> check no coordination bond
          check_pos(7) = 5 + modulo(pos-3,4)  !% target --> TPYP --> check no hydrogen
          check_pos(8) = 5 + modulo(pos-2,4)  !% target --> TPYP --> check no hydrogen
          !% check result
          check_result(1) = rcd(check_pos(1), sid)
          check_result(2) = rcd(check_pos(2), sid)
          check_result(3) = rcd(check_pos(3), sid)
          check_result(4) = rcd(check_pos(4), sid)
          check_result(5) = rcd(check_pos(5), sid)
          check_result(6) = rcd(check_pos(6), sid)
          check_result(7) = rcd(check_pos(7), tid)
          check_result(8) = rcd(check_pos(8), tid)
          if (sum(check_result) == 0) then
            if (rcd(modulo(pos-3,4)+9, sid)==0) then !% target_pos
              rate = RCOA
            else
              rate = RCOB
            end if
          end if
        end if
      end if
    end if
    if (rate == 0.0) tid = 0
  return
end subroutine form_bond_test