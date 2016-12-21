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
  subroutine def_tpyp_tpyp (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< TPYP class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir           !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: posN(4)=[2,3,4,5] !< 1HYB, 2HYB, 3HYB, COB
    integer :: k, r
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
    call obj % reac (r) % cond (1) % alloc_opt (8)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(2*k-1) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
               ([1,1,  2,3,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
               ([1,1,  posN(k),posN(k),  2,3])
       end if
       !* position #2
       call obj % reac (r) % cond (1) % opt(2*k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(2*k) % set_state &
               ([1,1,  3,4,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(2*k) % set_state &
               ([1,1,  posN(k),posN(k),  3,4])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (24)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-5) % set_pos ([+2,1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-5) % set_pos ([-2,1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-5) % set_state &
            ([1,1,  posN(k),posN(k),  2,3])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-4) % set_pos ([+2,-1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-4) % set_pos ([-2,-1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-4) % set_state &
            ([1,1,  posN(k),posN(k),  2,3])
       !* position #3
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-3) % set_pos &
               ([+2,1,0,  +2,-1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-3) % set_pos &
               ([-2,1,2,  -2,-1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-3) % set_state &
            ([1,1,  posN(k),posN(k),  3,4])
       !* position #4
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-2) % set_pos ([+2,1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k-2) % set_pos ([-2,1,0])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-2) % set_state &
            ([1,1,  2,3,  posN(k),posN(k)])
       !* position #5
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-1) % set_pos ([+2,-1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k-1) % set_pos ([-2,-1,0])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-1) % set_state &
            ([1,1,  2,3,  posN(k),posN(k)])
       !* position #6
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k) % set_pos &
               ([+2,1,2,  +2,-1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k) % set_pos &
               ([-2,1,0,  -2,-1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k) % set_state &
            ([1,1,  3,4,  posN(k),posN(k)])
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
    call obj % reac (r) % cond (1) % alloc_opt (8)
    do k = 1, 4
       !* position #1
       call obj % reac (r) % cond (1) % opt(2*k-1) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
               ([1,1,  3,2,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
               ([1,1,  posN(k),posN(k),  3,2])
       end if
       !* position #2
       call obj % reac (r) % cond (1) % opt(2*k) % set_pos ([0,0,0])
       if (dir == 1) then
          call obj % reac (r) % cond (1) % opt(2*k) % set_state &
               ([1,1,  4,3,  posN(k),posN(k)])
       else
          call obj % reac (r) % cond (1) % opt(2*k) % set_state &
               ([1,1,  posN(k),posN(k),  4,3])
       end if
    end do
    !-- #2 tptp
    call obj % reac (r) % cond (2) % set_tp (tpyp_idx)
    call obj % reac (r) % cond (2) % alloc_opt (24)
    !* possible positions
    do k = 1, 4
       !* position #1
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-5) % set_pos ([+2,1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-5) % set_pos ([-2,1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-5) % set_state &
            ([1,1,  posN(k),posN(k),  3,2])
       !* position #2
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-4) % set_pos ([+2,-1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-4) % set_pos ([-2,-1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-4) % set_state &
            ([1,1,  posN(k),posN(k),  3,2])
       !* position #3
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-3) % set_pos &
               ([+2,1,0,  +2,-1,0])
       else
          call obj % reac (r) % cond (2) % opt(6*k-3) % set_pos &
               ([-2,1,2,  -2,-1,2])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-3) % set_state &
            ([1,1,  posN(k),posN(k),  4,3])
       !* position #4
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-2) % set_pos ([+2,1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k-2) % set_pos ([-2,1,0])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-2) % set_state &
            ([1,1,  3,2,  posN(k),posN(k)])
       !* position #5
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k-1) % set_pos ([+2,-1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k-1) % set_pos ([-2,-1,0])
       end if
       call obj % reac (r) % cond (2) % opt(6*k-1) % set_state &
            ([1,1,  3,2,  posN(k),posN(k)])
       !* position #6
       if (dir == 1) then
          call obj % reac (r) % cond (2) % opt(6*k) % set_pos &
               ([+2,1,2,  +2,-1,2])
       else
          call obj % reac (r) % cond (2) % opt(6*k) % set_pos &
               ([-2,1,0,  -2,-1,0])
       end if
       call obj % reac (r) % cond (2) % opt(6*k) % set_state &
            ([1,1,  4,3,  posN(k),posN(k)])
    end do
    return
  end subroutine def_tpyp_tpyp

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> initialization
  !---------------------------------------------------------------------------
  subroutine def_tpyp_lead (obj, rid, dir, Eform, Ebreak)
    !-- arguments
    type(mtype), pointer, intent (inout) :: obj  !< class instance
    integer             , intent (in) :: rid     !< reaction index
    integer             , intent (in) :: dir           !< direction
    real(dp)            , intent (in) :: Eform, Ebreak !< movement energy
    !-- local variables
    integer :: r
    !---------------------------
    !-- function definition
    !
    !---------------------------
    !> bond formation
    ! 6(789) 121
    r = rid
    ! call obj % reac (r) % set_energy (Eform)
    ! call obj % reac (r) % set_move   ([0,0,0])
    ! call obj % reac (r) % alloc_cond (2)
    ! !-- #1 tpyp
    ! call obj % reac (r) % cond (1) % set_tp (tpyp_idx)
    ! call obj % reac (r) % cond (1) % alloc_opt (8)
    ! do k = 1, 4
    !    !* position #1
    !    call obj % reac (r) % cond (1) % opt(2*k-1) % set_pos ([0,0,0])
    !    if (dir == 1) then
    !       call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
    !            ([1,1,  2,5,  posN(k),posN(k)])
    !    else
    !       call obj % reac (r) % cond (1) % opt(2*k-1) % set_state &
    !            ([2,2,  posN(k),posN(k),  1,4])
    !    end if
    !    !* position #2
    !    call obj % reac (r) % cond (1) % opt(2*k) % set_pos ([0,0,0])
    !    if (dir == 1) then
    !       call obj % reac (r) % cond (1) % opt(2*k) % set_state &
    !            ([2,2,  4,5,  posN(k),posN(k)])
    !    else
    !       call obj % reac (r) % cond (1) % opt(2*k) % set_state &
    !            ([2,2,  posN(k),posN(k),  4,5])
    !    end if
    ! end do
    
    return 
  end subroutine def_tpyp_lead

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
    call tpyp % alloc_reac (15)
    ! free movements
    call def_free_move(tpyp, 1, [ 1, 0, 0], 0.68_dp)
    call def_free_move(tpyp, 2, [ 0, 1, 0], 0.68_dp)
    call def_free_move(tpyp, 3, [-1, 0, 0], 0.68_dp)
    call def_free_move(tpyp, 4, [ 0,-1, 0], 0.68_dp)
    call def_free_move(tpyp, 5, [ 0, 0, 1], 0.68_dp) ! rotate by  90 degrees
    call def_free_move(tpyp, 6, [ 0, 0, 2], 0.68_dp) ! rotate by 180 degrees
    call def_free_move(tpyp, 7, [ 0, 0, 3], 0.68_dp) ! rotate by 270 degrees
    ! reaction: tpyp+tpyp
    Eform  = 0.00_dp
    Ebreak = 2.68_dp
    call def_tpyp_tpyp (tpyp,  8, 1, Eform, Ebreak)
    call def_tpyp_tpyp (tpyp, 10, 2, Eform, Ebreak)
    call def_tpyp_tpyp (tpyp, 12, 1, Eform, Ebreak)
    call def_tpyp_tpyp (tpyp, 14, 2, Eform, Ebreak)
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    !
    !> define type Lead
    !
    call lead % set_symm (1) 
    call lead % set_idx_def (lead_idx)
    call lead % set_eva_num (0)   
    call lead % alloc_comp (1)
    call lead % set_comp (1, [0, 0, 6]) 
    call lead % alloc_reac (4)
    call def_free_move(lead, 1, [ 1, 0,0], 0.5_dp)
    call def_free_move(lead, 2, [ 0, 1,0], 0.5_dp)
    call def_free_move(lead, 3, [-1, 0,0], 0.5_dp)
    call def_free_move(lead, 4, [ 0,-1,0], 0.5_dp)
    !--------------------------------------------------

    !--------------------------------------------------
    !> Auto initialization
    call init_random_seed()    !< initialize random seed
    call init_substrate(20,20) !< substrate size (20x20 for window)
    call init_rates()          !< this should be placed after
    call mlist_init()          !  type definitions 
    !--------------------------------------------------
    
    return
  end subroutine init

end module define
