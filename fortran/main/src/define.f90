!---------------------------------------------------------------------------  
! DESCRIPTION
!> this is a module should be handle by user for defining all constant
!  variables, molecule structures and chemical reactions
!---------------------------------------------------------------------------  
module define
  
  use class_mtype
  use class_molecule
  use func_helper
  use func_substrate
  use func_rate_kmc
  implicit none

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> define molecule type here
  !--------------------------------------------------------------------------- 
  type (mtype), save, pointer :: tpyp
  type (mtype), save, pointer :: lead
  integer, save :: tpyp_idx = 2000
  integer, save :: lead_idx = 1000

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> define the hopping time
  integer, parameter :: HOP_TIME = 1E4
  !--------------------------------------------------------------------------- 

contains

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> This function can define free movement (left/right/up/down/rotation)
  !  relatively easily
  !--------------------------------------------------------------------------- 
  subroutine def_free_move(obj, r, move, energy)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< class instance
    integer             , intent (in) :: r       !< reaction index
    integer             , intent (in) :: move(3) !< movement vector [x,y,d]
    real(dp)            , intent (in) :: energy  !< movement energy
    !-- local variables
    integer, allocatable :: array(:,:)
    integer, allocatable :: tarpos(:) !< target position
    integer, allocatable :: astate(:) !< summary of all state changes
    integer, pointer     :: istate(:) !< initial states for each position
    integer              :: i, j, n
    !
    !-- function definition
    n = obj % comp_num() !< number of molecule size
    !< define a temporary array to store information 5xN where N is the size
    !  of the molecule
    !< the five elements are
    !  [new-x, new-y, flag, ori-x, ori-y]
    call alloc(array, 5, n)
    array = 1 !< initialize array value with one
    !< now we need to find all locations that 
    do i = 1, n
       array(1:2,i) = obj%translate(obj%comp(i), move)
       !< check overlapping
       do j = 1, n
          array(4:5,j) = obj%xy(j)
          if (all(array(1:2,i) == array(4:5,j))) array(3,i) = 0
       end do
    end do
    !
    !< the array to store all target positions
    call alloc_I1(tarpos, sum(array(3,:)) * 3)
    !< the array to store all state information
    call alloc_I1(astate, n * 2)
    !
    !< fill the arraies 
    j = 1 !< this is a counter
    istate => obj%state_initial() !< initial states
    do i = 1, n
       if (array(3,i) == 1) then
          tarpos(3*j-2) = array(1,i)
          tarpos(3*j-1) = array(2,i)
          tarpos(3*j  ) = 0
          j = j + 1
       end if
       astate(2*i-1) = istate(i) !< donot need to change state for free move
       astate(2*i  ) = istate(i) !< donot need to change state for free move
    end do
    !
    ! basic information
    call obj%reac(r)%set_energy (energy)
    call obj%reac(r)%set_move   (move)
    ! two conditions
    call obj%reac(r)%alloc_cond (2)
    ! condition for molecule itself
    call obj%reac(r)%cond(1)%set_tp (obj%idx_def())
    call obj%reac(r)%cond(1)%alloc_opt (1)
    call obj%reac(r)%cond(1)%opt(1)%set_pos ([0,0,0])
    call obj%reac(r)%cond(1)%opt(1)%set_state (astate)
    ! condition for background checking (empty checking)
    call obj%reac(r)%cond(2)%set_tp (0) !< this is the index for background
    call obj%reac(r)%cond(2)%alloc_opt (1)
    call obj%reac(r)%cond(2)%opt(1)%set_pos (tarpos)
    call obj%reac(r)%cond(2)%opt(1)%set_state ([0,0]) 
    return
  end subroutine def_free_move

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> this function defines the hydrogen bond between tpyps
  !--------------------------------------------------------------------------- 
  subroutine def_tpyp_1tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< TPYP class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir           !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: k, r, s, c
    !
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  2,3,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  posN(k),posN(k),  2,3])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (32)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,0])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,2])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s),posN(s+1)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,0])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,2])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s),posN(s+1)])
          !* position #3
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s),posN(s+1),  posN(k),posN(k)])
          !* position #4
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s),posN(s+1),  posN(k),posN(k)])
       end do
    end do
    !
    !---------------------------    
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  3,2,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  posN(k),posN(k),  3,2])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (32)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,0])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,2])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s+1),posN(s)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,0])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,2])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s+1),posN(s)])
          !* position #3
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s+1),posN(s),  posN(k),posN(k)])
          !* position #4
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s+1),posN(s),  posN(k),posN(k)])
       end do
    end do
    return
  end subroutine def_tpyp_1tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> this function defines the hydrogen bond between tpyps
  !--------------------------------------------------------------------------- 
  subroutine def_tpyp_2tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< TPYP class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir           !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: k, r, c, s
    !
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (3)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  3,4,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  posN(k),posN(k),  3,4])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (16)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s),posN(s+1)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,+1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,+1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s),posN(s+1),  posN(k),posN(k)])
       end do
    end do
    !-- #3 tptp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (16)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,-1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,-1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,+1,2])
          else
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,+1,0])
          end if
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s+1),posN(s+1)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,+1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,+1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,-1,2])
          else
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,-1,0])
          end if
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1,  posN(s+1),posN(s+1),  posN(k),posN(k)])
       end do
    end do
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (3)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  4,3,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  posN(k),posN(k),  4,3])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (16)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s+1),posN(s)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,-1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,-1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([+2,+1,2])
          else
             call obj % reac (r) % cond (2) % opt(c) % set_pos ([-2,+1,0])
          end if
          call obj % reac (r) % cond (2) % opt(c) % set_state &
               ([1,1,  posN(s+1),posN(s),  posN(k),posN(k)])
       end do
    end do
    !-- #3 tptp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (16)
    !* possible positions
    c = 0
    do k = 1, 4
       do s = 1, 2
          !* position #1
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,-1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,-1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,+1,2])
          else
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,+1,0])
          end if
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1,  posN(k),posN(k),  posN(s+1),posN(s+1)])
          !* position #2
          c = c + 1
          if (dir == 1) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,+1,0])
          else if (dir == 2) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,+1,2])
          else if (dir == 3) then
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,-1,2])
          else
             call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,-1,0])
          end if
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1,  posN(s+1),posN(s+1),  posN(k),posN(k)])
       end do
    end do
    return
  end subroutine def_tpyp_2tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> reaction tpyp+lead
  !---------------------------------------------------------------------------
  subroutine def_tpyp_lead (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir           !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: r, k
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  2,5,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([2,2,  posN(k),posN(k),  2,5])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (lead_idx)
    call obj % reac (r) % cond (2) % alloc_opt (4)
    !* possible positions
    !* position #1
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(1) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(1) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(1) % set_state ([6,7])
    !* position #2
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(2) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(2) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(2) % set_state ([7,8])
    !* position #3
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(3) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(3) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(3) % set_state ([8,9])
    !* position #4
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(4) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(4) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(4) % set_state ([9,10])
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (1) % alloc_opt (4)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  5,2,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(k) % set_state &
               ([1,1,  posN(k),posN(k),  5,2])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (lead_idx)
    call obj % reac (r) % cond (2) % alloc_opt (4)
    !* possible positions
    !* position #1
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(1) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(1) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(1) % set_state ([7,6])
    !* position #2
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(2) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(2) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(2) % set_state ([8,7])
    !* position #3
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(3) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(3) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(3) % set_state ([9,8])
    !* position #4
    if (dir == 1) then
       call obj % reac (r) % cond (2) % opt(4) % set_pos ([+2,0,0])
    else
       call obj % reac (r) % cond (2) % opt(4) % set_pos ([-2,0,0])
    end if
    call obj % reac (r) % cond (2) % opt(4) % set_state ([10,9])
    return 
  end subroutine def_tpyp_lead

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> reaction lead+tpyp
  !---------------------------------------------------------------------------
  subroutine def_lead_1tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< LEAD class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir     !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: r, k
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([6,7])
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 2,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 2,5])
    end do
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (2)
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([7,6])
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 5,2, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,2])
    end do
    return
  end subroutine def_lead_1tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> reaction lead with 2 tpyp
  !---------------------------------------------------------------------------
  subroutine def_lead_2tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< LEAD class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir     !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: r, k, c
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (3) ! one lead and two tpyps
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([7,8])
    !-- #2 tpyp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 2,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 2,5])
    end do
    !-- #3 tpyp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (24)
    !* possible positions
    do k = 1, 4
       !* position #1
       c = 6*k-5
       if (dir == 1) then      !< [-2,0,0]          
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3])
       else if (dir == 2) then !< [0,-2,1]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0])
       else if (dir == 3) then !< [+2,0,2]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1])
       else                    !< [0,+2,3]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2])
       end if
       do c = 6*k-5, 6*k-3
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1, 5,5, posN(k),posN(k)])
       end do
       !* position #2
       c = 6*k-2
       if (dir == 1) then      !< [+2,0,0]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3])
       else if (dir == 2) then !< [0,+2,1]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0])
       else if (dir == 3) then !< [-2,0,2]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1])
       else                    !< [0,-2,3]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2])
       end if
       do c = 6*k-2, 6*k
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1, posN(k),posN(k), 5,5])
       end do
    end do
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (3)
    !-- #1 tpyp
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([8,7])
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 5,2, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,2])
    end do
    !-- #3 tpyp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (24)
    !* possible positions
    do k = 1, 4
       !* position #1
       c = 6*k-5
       if (dir == 1) then      !< [-2,0,0]          
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3])
       else if (dir == 2) then !< [0,-2,1]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0])
       else if (dir == 3) then !< [+2,0,2]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1])
       else                    !< [0,+2,3]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,2])
       end if
       do c = 6*k-5, 6*k-3
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1, 5,5, posN(k),posN(k)])
       end do
       !* position #2
       c = 6*k-2
       if (dir == 1) then      !< [+2,0,0]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3])
       else if (dir == 2) then !< [0,+2,1]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0])
       else if (dir == 3) then !< [-2,0,2]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,-2,3]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1])
       else                    !< [0,-2,3]
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([+2,0,0]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([0,+2,1]); c=c+1
          call obj % reac (r) % cond (3) % opt(c) % set_pos ([-2,0,2])
       end if
       do c = 6*k-2, 6*k
          call obj % reac (r) % cond (3) % opt(c) % set_state &
               ([1,1, posN(k),posN(k), 5,5])
       end do
    end do
    return
  end subroutine def_lead_2tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> reaction lead with 3 tpyp
  !---------------------------------------------------------------------------
  subroutine def_lead_3tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< LEAD class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir     !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: r, k, c
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (4) ! one lead and two tpyps
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([8,9])
    !-- #2 tpyp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 2,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 2,5])
    end do
    !-- #3 tpyp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,+2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([-2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,-2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([+2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #4 tpyp
    call obj % reac (r) % cond (4) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (4) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([-2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,-2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([+2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,+2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (4) ! one lead and two tpyps
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([9,9])
    !-- #2 tpyp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 5,2, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,2])
    end do
    !-- #3 tpyp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,+2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([-2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,-2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([+2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #4 tpyp
    call obj % reac (r) % cond (4) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (4) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([-2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,-2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([+2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,+2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    return
  end subroutine def_lead_3tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> reaction lea 4 tpyp
  !---------------------------------------------------------------------------
  subroutine def_lead_4tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< LEAD class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir     !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: r, k, c
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    r = rid
    call obj % reac (r) % set_energy (Eform)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (5)
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([9,10])
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 2,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 2,5])
    end do
    !-- #3 tptp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,+2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([-2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,-2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([+2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #4 tptp
    call obj % reac (r) % cond (4) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (4) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([-2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,-2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([+2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,+2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #5 tptp
    call obj % reac (r) % cond (5) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (5) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 2) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 3) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([0,-2,1])
       else
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([+2,0,2])
       end if
       call obj % reac (r) % cond (5) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 2) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 3) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([0,+2,1])
       else
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([-2,0,2])
       end if
       call obj % reac (r) % cond (5) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !
    !---------------------------
    !> bond break
    r = rid + 1
    call obj % reac (r) % set_energy (Ebreak)
    call obj % reac (r) % set_move   ([0,0,0])
    call obj % reac (r) % alloc_cond (5)
    !-- #1 lead
    call obj % reac (r) % cond (1) % set_tp (lead_idx)
    call obj % reac (r) % cond (1) % alloc_opt (1)
    call obj % reac (r) % cond (1) % opt(1) % set_pos ([0,0,0])
    call obj % reac (r) % cond (1) % opt(1) % set_state ([10,9])
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([+2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k-1) % set_pos ([0,+2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k-1) % set_state &
            ([1,1, 5,2, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 2) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 3) then
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([-2,0,2])
       else
          call obj % reac (r) % cond (2) % opt(2*k) % set_pos ([0,-2,3])
       end if
       call obj % reac (r) % cond (2) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,2])
    end do
    !-- #3 tptp
    call obj % reac (r) % cond (3) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (3) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,-2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([0,+2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k-1) % set_pos ([-2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,+2,1])
       else if (dir == 2) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 3) then
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([0,-2,3])
       else
          call obj % reac (r) % cond (3) % opt(2*k) % set_pos ([+2,0,0])
       end if
       call obj % reac (r) % cond (3) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #4 tptp
    call obj % reac (r) % cond (4) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (4) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([+2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([-2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k-1) % set_pos ([0,-2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([-2,0,2])
       else if (dir == 2) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 3) then
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([+2,0,0])
       else
          call obj % reac (r) % cond (4) % opt(2*k) % set_pos ([0,+2,1])
       end if
       call obj % reac (r) % cond (4) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do
    !-- #5 tptp
    call obj % reac (r) % cond (5) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (5) % alloc_opt (8)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([0,+2,3])
       else if (dir == 2) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([-2,0,0])
       else if (dir == 3) then
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([0,-2,1])
       else
          call obj % reac (r) % cond (5) % opt(2*k-1) % set_pos ([+2,0,2])
       end if
       call obj % reac (r) % cond (5) % opt(2*k-1) % set_state &
            ([1,1, 5,5, posN(k),posN(k)])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([0,-2,3])
       else if (dir == 2) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([+2,0,0])
       else if (dir == 3) then
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([0,+2,1])
       else
          call obj % reac (r) % cond (5) % opt(2*k) % set_pos ([-2,0,2])
       end if
       call obj % reac (r) % cond (5) % opt(2*k) % set_state &
            ([1,1, posN(k),posN(k), 5,5])
    end do

    return
  end subroutine def_lead_4tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> initialization
  !---------------------------------------------------------------------------
  subroutine init ()
    type(mtype), pointer :: obj
    real(dp)             :: Eform, Ebreak

    !-------------------------------------------------------------------
    ! handled by user
    !
    !> adding types into the record
    !> DOC
    !  Here you need to define how many kinds of molecule you need
    !  for example I defined two kinds of molecules here: TPYP and LEAD
    call tlist_init (2)
    tpyp => tlist_new ()
    lead => tlist_new ()
    !-------------------------------------------------------------------
    write (*,*) 'TPYP'
    !-------------------------------------------------------------------
    !
    !> define molecule type TPyP
    !
    ! how to check one condition
    ! 1) check if one of the relative positions matches: cond->opt->pos
    ! 2) check if one of the relative direction matches: cond->opt->dir
    ! 3) check if ALL components' initial states match : cond->state
    ! pass the checking and do movement and update component state
    !
    call tpyp % set_symm    (4) !> define rotation angle (180/n)
    call tpyp % set_idx_def (tpyp_idx) !> type id should be within [1000, 9999]
    call tpyp % set_eva_num (200) !> evaporation number
    call tpyp % alloc_comp (3)    !> number of components
    !
    !> define molecule structure
    call tpyp % set_comp(1, [ 0, 0, 1]) !> x, y, initial-state
    call tpyp % set_comp(2, [ 1, 0, 2]) !> x, y, initial-state
    call tpyp % set_comp(3, [-1, 0, 2]) !> x, y, initial-state
    !
    !> define reactions
    call tpyp % alloc_reac (23)
    ! free movements
    call def_free_move(tpyp, 1, [ 1, 0, 0], 0.68_dp)
    call def_free_move(tpyp, 2, [ 0, 1, 0], 0.68_dp)
    call def_free_move(tpyp, 3, [-1, 0, 0], 0.68_dp)
    call def_free_move(tpyp, 4, [ 0,-1, 0], 0.68_dp)
    call def_free_move(tpyp, 5, [ 0, 0, 1], 0.68_dp) ! rotate by  90 degrees
    call def_free_move(tpyp, 6, [ 0, 0, 2], 0.68_dp) ! rotate by 180 degrees
    call def_free_move(tpyp, 7, [ 0, 0, 3], 0.68_dp) ! rotate by 270 degrees
    ! reaction: tpyp+tpyp
    Eform  = 0.68_dp
    Ebreak = 2.68_dp
    write (*,*) 'tpyp+1tpyp'
    call def_tpyp_1tpyp (tpyp,  8, 1, Eform, Ebreak)
    call def_tpyp_1tpyp (tpyp, 10, 2, Eform, Ebreak)
    write (*,*) 'tpyp+2tpyp'
    call def_tpyp_2tpyp (tpyp, 12, 1, Eform, Ebreak)
    call def_tpyp_2tpyp (tpyp, 14, 2, Eform, Ebreak)
    call def_tpyp_2tpyp (tpyp, 16, 3, Eform, Ebreak)
    call def_tpyp_2tpyp (tpyp, 18, 4, Eform, Ebreak)
    ! reaction: tpyp+lead
    Eform  = 0.68_dp
    Ebreak = 2.68_dp
    write (*,*) 'tpyp+lead'
    call def_tpyp_lead (tpyp, 20, 1, Eform, Ebreak)
    call def_tpyp_lead (tpyp, 22, 2, Eform, Ebreak)
    !-------------------------------------------------------------------
    write (*,*) 'LEAD'
    !-------------------------------------------------------------------
    !
    !> define type Lead
    !
    call lead % set_symm (1) 
    call lead % set_idx_def (lead_idx)
    call lead % set_eva_num (100)   
    call lead % alloc_comp (1)
    call lead % set_comp (1, [0, 0, 6]) 
    call lead % alloc_reac (36)
    call def_free_move(lead, 1, [ 1, 0,0], 0.68_dp)
    call def_free_move(lead, 2, [ 0, 1,0], 0.68_dp)
    call def_free_move(lead, 3, [-1, 0,0], 0.68_dp)
    call def_free_move(lead, 4, [ 0,-1,0], 0.68_dp)
    ! reaction: lead+tpyp
    Eform  = 0.00_dp
    Ebreak = 2.68_dp
    write (*,*) 'lead+1tpyp'
    call def_lead_1tpyp (lead, 5, 1, Eform, Ebreak)
    call def_lead_1tpyp (lead, 7, 2, Eform, Ebreak)
    call def_lead_1tpyp (lead, 9, 3, Eform, Ebreak)
    call def_lead_1tpyp (lead,11, 4, Eform, Ebreak)
    write (*,*) 'lead+2tpyp'
    call def_lead_2tpyp (lead,13, 1, Eform, Ebreak)
    call def_lead_2tpyp (lead,15, 2, Eform, Ebreak)
    call def_lead_2tpyp (lead,17, 3, Eform, Ebreak)
    call def_lead_2tpyp (lead,19, 4, Eform, Ebreak)
    write (*,*) 'lead+3tpyp'
    call def_lead_3tpyp (lead,21, 1, Eform, Ebreak)
    call def_lead_3tpyp (lead,23, 2, Eform, Ebreak)
    call def_lead_3tpyp (lead,25, 3, Eform, Ebreak)
    call def_lead_3tpyp (lead,27, 4, Eform, Ebreak)
    write (*,*) 'lead+4tpyp'
    call def_lead_4tpyp (lead,29, 1, Eform, Ebreak)
    call def_lead_4tpyp (lead,31, 2, Eform, Ebreak)
    call def_lead_4tpyp (lead,33, 3, Eform, Ebreak)
    call def_lead_4tpyp (lead,35, 4, Eform, Ebreak)
    !--------------------------------------------------

    !--------------------------------------------------
    !> Auto initialization
    call init_random_seed()    !< initialize random seed
    call init_substrate(40,40) !< substrate size (20x20 for window)
    call init_rates()          !< this should be placed after
    call mlist_init()          !  type definitions 
    !--------------------------------------------------
    
    return
  end subroutine init

end module define
