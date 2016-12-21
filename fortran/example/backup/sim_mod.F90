module sim_mod
use glob_mod
use rand_mod
use rate_mod
use mono_mod
use rcd_mod
use sys_mod
use output_mod
implicit none

integer, parameter :: moves(2,4) = reshape([1,0, 0,1, -1,0, 0,-1], [2,4])
integer, save :: record_iter_list(18) = 0
integer, save :: record_iter_len = 0

real(8), save :: RCOA, RCOB, RCOC
real(8), save :: RBCA, RBCB, RBCC
real(8), save :: REME, RBME
real(8), save :: RHYB, RBHY 
real(8), save :: RDEM, RDEP

contains

subroutine sim_single(temp)
  integer :: ntot
  real(8) :: time, temp
    ntot = 0
    time = 0.0
    call new_out_folder([0,0,0])
    call sub_init()
    call mono_init()
    call rcd_init()
    call rate_init()
    call init_param(temp)
    call evaporate(TPYP, ntot)
    call evaporate(LEAD, ntot) 
    call print_sub()
    call hopping(no_of_mono, time, temp)
    print *, "delta time =", time
  return
end subroutine sim_single

subroutine init_param(temp)
  real(8), intent(in) :: temp
    RCOA = V0 * exp(COBH/2/KB/temp)
    RCOB = V0 * exp((COBF-COBH)/2/KB/temp)
    RCOC = V0 * exp(COBF/2/KB/temp)
    RBCA = V0 * exp(-COBH/2/KB/temp)
    RBCB = V0 * exp(-(COBF-COBH)/2/KB/temp)
    RBCC = V0 * exp(-COBF/2/KB/temp)
    REME = V0 * exp(EME/2/KB/temp)
    RBME = V0 * exp(-EME/2/KB/temp)
    RHYB = V0 * exp(HYB/2/KB/temp)
    RBHY = V0 * exp(-HYB/2/KB/temp)
    RDEM = V0 * exp(-DEM/KB/temp)
    RDEP = V0 * exp(-DEP/KB/temp)
  return
end subroutine init_param

subroutine hopping(ntot, time, temp)
  integer, intent(in) :: ntot
  integer, pointer :: input_list(:), output_list(:)
  integer :: i
  real(8) :: time, temp
    !=> initialization
    if (associated(input_list)) deallocate(input_list)
    allocate(input_list(ntot))
    input_list = [(i, i=1, ntot)]
    !=> single hop
    do i = 1, HOPSTP
      output_list => hop_single(input_list, time, temp)
      deallocate(input_list)
      input_list => output_list
    end do
  return
end subroutine hopping

function hop_single(list, time, temp) result(next_list)
  integer, pointer :: list(:) !< input iteration array
  integer, pointer :: next_list(:) !< output iteration array for the next hop
  integer, allocatable :: sur(:,:) !< surrounding points array
  integer, allocatable :: bor(:,:) !< border points array
  integer :: info(12) !< state info => #1(xcor) #2(ycor) #3(dir) #4(name) #5-12(states)
  integer :: datas(6) !< datas => #1(index) #2(xcor) #3(ycor) #4(dir) #5(state-pos) #6(target index)
  integer :: i, j, sid, tid, name
  real(8) :: rate, temp, time

    !=> update rates
    do i = 1, size(list)
      !=> get basic information
      sid = list(i) !< self index
      info(:) = rcd(:, sid)
      datas(1) = sid
      datas(2:4) = info(1:3)
      !=> load surrounding array
      call get_sur(info(1:2), info(3), info(4), sur)
      !==> iterate stats
      do j = 1, 8
        !=> form bond test --> change 0 state to be a index -> occupied rcd 4-8
        if (info(j+4) == 0) then
          call form_bond_test(sid, j, sur(j,:), rate, tid)
          datas(5) = j+4
          datas(6) = tid
        !=> break bond test --> change non-0 state to be 0 -> occupied rcd 9-12
        else if (info(j+4) /= 0) then
          call break_bond_test(sid, info(j+4), j, rate)
          datas(5) = j+4
          datas(6) = 0
        end if
        call rate_set(rate_wid*(sid-1)+j, rate, datas) !< push rate
      end do
      !==> iterate moves -> occupied rcd 13-16
      do j = 1, 4
        if (sum(info(5:12)) == 0) then
          datas(2) = modulo(moves(1,j)+info(1)-1, SIDE)+1
          datas(3) = modulo(moves(2,j)+info(2)-1, SIDE)+1
          call move_test(datas(2:3), info(3), info(4), rate, sid)
        else
          rate = 0.0
        end if
        datas(5) = j+12
        datas(6) = 0 !=> set all target index to be 0
        call rate_set(rate_wid*(sid-1)+j+8, rate, datas) !< push rate
      end do
    end do
    deallocate(sur)

    !=> select event
    datas = rate_select() !< select event
    sid = datas(1) !< selected index
    call clear_iter()
    call push_iter(sid)
    !==> move monomer
    if (datas(5) > 12) then
      info(:) = rcd(:, sid)
      name = info(4) !=> get name
      !=> update list at original place (use target list)
      call get_bor(info(1:2), info(3), name, bor)
      do i = 1, max_bor(name)
        call push_iter(get_sub(bor(i,1), bor(i,2)))
      end do
      !=> update list at the new place (use self list)
      call get_bor(datas(2:3), datas(4), name, bor)
      do i = 1, max_bor(name)
        call push_iter(get_sub(bor(i,1), bor(i,2)))
      end do
      !=> move place
      call deland(info(1:2), info(3), name, sid)
      call land(datas(2:3), datas(4), name, sid)
      call rcd_set_pos(sid, name, datas(2:3), datas(4)) !=> update record
    !==> change bonding states
    else
      tid = datas(6) !=> target index
      !=> form bonds
      if (tid /= 0) then
        rcd(datas(5), sid) = tid
        if (datas(5) > 8) then
          rcd(modulo(datas(5)-3,4)+9, tid) = sid
        else
          rcd(modulo(datas(5)-3,4)+5, tid) = sid
        end if
      !=> break bonds
      else
        tid = rcd(datas(5), sid)
        rcd(datas(5), sid) = 0
        if (datas(5) > 8) then
          rcd(modulo(datas(5)-3,4)+9, tid) = 0
        else
          rcd(modulo(datas(5)-3,4)+5, tid) = 0
        end if
      end if
      !=> update list at target
      info(:) = rcd(:, tid)
      name = info(4) !=> get name
      call get_bor(info(1:2), info(3), name, bor)
      do i = 1, max_bor(name)
        call push_iter(get_sub(bor(i,1), bor(i,2)))
      end do
      !=> update list at self
      info(:) = rcd(:, sid)
      name = info(4) !=> get name
      call get_bor(info(1:2), info(3), name, bor)
      do i = 1, max_bor(name)
        call push_iter(get_sub(bor(i,1), bor(i,2)))
      end do
    end if
    deallocate(bor)
    !=> update next list
    next_list => get_iter()
    time = rate_update_time(time)
  return
end function hop_single

subroutine clear_iter()
  record_iter_len = 0
  record_iter_list = 0
  return
end subroutine clear_iter

subroutine push_iter(n)
  integer :: n, i
  logical :: repeat
    repeat = .false.
    if (n > 5) then
      continue
    end if
    if (n > 5) then
      continue
    end if
    do i = 1, record_iter_len
      if (record_iter_list(i) == n .or. n == 0) repeat = .true.
    end do
    if (.not. repeat) then
      record_iter_len = record_iter_len + 1
      record_iter_list(record_iter_len) = n
    end if
  return
end subroutine push_iter

function get_iter()
  integer, pointer :: get_iter(:)
    allocate(get_iter(record_iter_len))
    get_iter = record_iter_list(1:record_iter_len)
  return
end function get_iter

subroutine evaporate(name, ntot)
  integer, intent(in) :: name
  integer :: ntot !< total number of monomer
  integer :: i, n !< total number of evaporation
    !% evaporation
    n = min(comps(name)%flu, comps(name)%num-eva_log(name))
    no_of_mono = no_of_mono - ntot
    do i = 1, n
      ntot = ntot + 1 !% increase end index
      eva_log(name) = eva_log(name) + 1 !% increase evaporation log by 1
      call eva_single(name, ntot)
    end do
    no_of_mono = no_of_mono + ntot
  return
end subroutine evaporate

subroutine eva_single (name, id)
  integer, intent(in) :: name, id
  integer :: pos(2), d
    EVA_LOOP: do while (.true.)
      pos(1) = rand_int(1, SIDE) !% random x coordinate
      pos(2) = rand_int(2, SIDE) !% random y coordinate
      d = rand_int(0, SYMM-1) !% direction
      if (.not. test_occupied(pos, d, name)) exit EVA_LOOP
    end do EVA_LOOP
    call land(pos, d, name, id)
    call rcd_set_pos(id, name, pos, d)
  return
end subroutine eva_single

! --> debug subroutine
subroutine eva_at (name, id, pos)
  integer, intent(in) :: name, id, pos(2)
  integer :: d
    d = 0
    call land(pos, d, name, id)
    call rcd_set_pos(id, name, pos, d)
  return
end subroutine eva_at

end module sim_mod
