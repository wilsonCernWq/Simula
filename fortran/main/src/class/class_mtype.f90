!-----------------------------------------------------------------------------
!
! DESCRIPTION
!> @brief This class defines all the common properties shared by the same 
!         molecule type
!
! FIELD
!  
!-----------------------------------------------------------------------------
module class_mtype

  ! used modules
  use global        , only: SYS_SYMM
  use func_helper   , only: alloc, rotate
  use class_reaction, only: reaction
  implicit none
  private
  
  ! moledule type class
  type, public :: mtype
     logical, private :: m_init_status (7)
     ! part [1]
     integer, public  :: symm         ! number of symmetry (rotation)
     integer, private :: m_rmat (2,2) ! rotation matrix (auto)
     ! part [2]
     integer, private :: m_eva_num    ! evaporation number
     ! part [3]
     integer, private :: m_idx_def    ! user defined index for reference
     ! part [4]
     integer, private :: m_idx_off    ! molecule index offset (auto)
     ! part [5]
     integer, private :: m_idx_gen    ! data index in storage (auto)
     ! part [6]
     integer, private, pointer :: m_comp (:,:) => null() ! 3xN {x, y, si}
     integer, private          :: m_comp_num             ! component number
     ! part [7]
     type (reaction), public , pointer ::   reac (:) => null() ! reaction list
     integer        , private          :: m_reac_num           ! reaction number
   contains
     procedure :: check    => m_check
     ! getters
     procedure :: eva_num  => m_get_eva_num
     procedure :: idx_def  => m_get_idx_def
     procedure :: idx_off  => m_get_idx_off
     procedure :: idx_gen  => m_get_idx_gen
     procedure :: state_initial  => m_get_state_initial
     procedure :: x  => m_get_x
     procedure :: y  => m_get_y
     procedure :: xy => m_get_xy
     procedure :: comp     => m_get_comp
     procedure :: comp_num => m_get_comp_num
     procedure :: reac_num => m_get_reac_num
     ! setters / allocators
     procedure :: set_symm     => m_set_symm
     procedure :: set_eva_num  => m_set_eva_num
     procedure :: set_idx_def  => m_set_idx_def
     procedure :: set_idx_gen  => m_set_idx_gen
     procedure :: set_idx_off  => m_set_idx_off
     procedure :: set_comp     => m_set_comp
     procedure :: alloc_comp   => m_alloc_comp
     procedure :: alloc_reac   => m_alloc_reac
     ! calculator
     procedure :: abs_id       => m_abs_idx
     procedure :: rel_id       => m_rel_idx
     procedure :: rotate       => m_rotate
     procedure :: translate    => m_translate
  end type mtype

  !---------------------------------------------------------------------------  
  !> global type list
  type (mtype), public, allocatable, target :: tlist (:)

  !---------------------------------------------------------------------------  
  !> global functions
  public  :: tlist_num, tlist_init, tlist_new

  !---------------------------------------------------------------------------  
  !> private functions (just a reminder)
  private :: tlist_build_background

contains

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief check if type definition is completed
  !---------------------------------------------------------------------------  
  subroutine m_check(this)
    class(mtype), intent (in) :: this
    if (.not. all(this % m_init_status)) then
       print *, this % m_init_status
       stop "ERROR: Molecule type definition imcomplete"
    end if
    return
  end subroutine m_check
  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for molecule rotational symmetry (range for direction, or in
  !         the other words, number of 90/30 degrees it needs for completing 
  !         one full rotation) and compute corresponding rotational matrix
  !> @param n: symmetry number
  !---------------------------------------------------------------------------  
  subroutine m_set_symm (this, n)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: n
    this % symm   = n
    this % m_rmat = 0
    if (SYS_SYMM == 4) then
       if (n == 1) then
          this % m_rmat (1,1) =  1
          this % m_rmat (2,2) =  1
       else if (n == 2) then
          this % m_rmat (1,1) = -1
          this % m_rmat (2,2) = -1
       else if (n == 4) then
          this % m_rmat (1,2) = -1
          this % m_rmat (2,1) =  1
       else
          stop "ERROR: Invalid symmetry"
       end if
    else 
       stop "ERROR: unknown system symmetry"
    end if
    this % m_init_status(1) = .true.
    return
  end subroutine m_set_symm

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter type evaporation amount
  !> @param n: amount
  !---------------------------------------------------------------------------  
  subroutine m_set_eva_num (this, n)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: n
    this % m_eva_num = n
    this % m_init_status(2) = .true.
    return
  end subroutine m_set_eva_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for molecule evaporation number
  !---------------------------------------------------------------------------  
  function m_get_eva_num (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_eva_num
    return 
  end function m_get_eva_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for user defined type index
  !> @param i: index
  !---------------------------------------------------------------------------  
  subroutine m_set_idx_def (this, i)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: i
    this % m_idx_def = i
    this % m_init_status(3) = .true.
    return
  end subroutine m_set_idx_def

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for defined type index
  !---------------------------------------------------------------------------  
  function m_get_idx_def (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_idx_def
    return 
  end function m_get_idx_def

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for program generated type index
  !> @param i: index
  !---------------------------------------------------------------------------  
  subroutine m_set_idx_gen (this, i)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: i
    this % m_idx_gen = i
    this % m_init_status(4) = .true.
    return
  end subroutine m_set_idx_gen
  
  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for generated type index
  !---------------------------------------------------------------------------  
  function m_get_idx_gen (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_idx_gen
    return 
  end function m_get_idx_gen

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for molecule index offset
  !> @param i: index
  !---------------------------------------------------------------------------  
  subroutine m_set_idx_off (this, i)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: i
    this % m_idx_off = i
    this % m_init_status(5) = .true.
    return
  end subroutine m_set_idx_off

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for molecule index offset
  !---------------------------------------------------------------------------  
  function m_get_idx_off (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_idx_off
    return 
  end function m_get_idx_off

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Set component number and allocate comp array
  !> @param n: maximum component (dot) number
  !---------------------------------------------------------------------------  
  subroutine m_alloc_comp (this, n)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: n
    call alloc (this % m_comp, 3, n)
    this % m_comp_num = n
    this % m_comp     = 0
    this % m_init_status(6) = .true.
    return
  end subroutine m_alloc_comp

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for component
  !> @param i: index
  !> @param v: value [x, y, initial-state]
  !---------------------------------------------------------------------------  
  subroutine m_set_comp (this, i, v)
    class(mtype), intent (inout) :: this
    integer     , intent (in) :: i, v(3)
    this % m_comp(:,i) = v
    return 
  end subroutine m_set_comp

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for component information
  !> @param i: index
  !---------------------------------------------------------------------------  
  function m_get_comp (this, i) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: i
    integer                   :: r(3)
    r = this % m_comp(:,i)
    return 
  end function m_get_comp

  function m_get_x (this, i) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: i
    integer                   :: r
    r = this % m_comp(1,i)
    return 
  end function m_get_x

  function m_get_y (this, i) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: i
    integer                   :: r
    r = this % m_comp(2,i)
    return 
  end function m_get_y

  function m_get_xy (this, i) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: i
    integer                   :: r(2)
    r = this % m_comp(1:2,i)
    return 
  end function m_get_xy

  function m_get_state_initial (this) result (r)
    class(mtype), intent (in) :: this
    integer, pointer          :: r(:)
    r => this % m_comp(3,:)
    return 
  end function m_get_state_initial

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for component number
  !---------------------------------------------------------------------------  
  function m_get_comp_num (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_comp_num
    return 
  end function m_get_comp_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Set total number of reactions and allocate reactions array
  !> @param n: reaction number
  !---------------------------------------------------------------------------  
  subroutine m_alloc_reac (this, n)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: n
    integer                      :: i, status
    ! multiple allocation check
    if (associated(this % reac)) stop "ERROR: (reac) multiple definitions"
    ! assign value
    this % m_reac_num = n
    ! not a basic type, allocate manually
    allocate (this % reac (n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    ! initialization
    do i = 1, n
       call this % reac (i) % set_rid (i)
    end do
    this % m_init_status(7) = .true.
    return
  end subroutine m_alloc_reac

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for reaction number
  !---------------------------------------------------------------------------  
  function m_get_reac_num (this) result (r)
    class(mtype), intent (in) :: this
    integer                   :: r
    r = this % m_reac_num
    return 
  end function m_get_reac_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for tlist size
  !---------------------------------------------------------------------------
  function tlist_num () result (r)
    integer :: r
    r = size(tlist) - 1
    return
  end function tlist_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Allocate tlist
  !> @param n: number
  !> @remark allocate tlist from 0 to N, where #0 is prepared for the 
  !          pre-defiend background type
  !---------------------------------------------------------------------------  
  subroutine tlist_init (n)
    integer, intent (in) :: n
    integer              :: status
    ! check multiple allocation
    if (allocated(tlist)) stop "ERROR: (alloc tlist) multiple definitions"
    ! not a basic type, allocate manually
    allocate (tlist (0:n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    ! insert background
    call tlist_build_background()
    return
  end subroutine tlist_init

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief add a new type (by reference) to the tlist
  !> @param tp: the new molecule type to be added
  !> @remark intert type from index 1, index zero should be allocate 
  !          seperately using other function
  !---------------------------------------------------------------------------  
  function tlist_new () result (r)
    type (mtype), pointer  :: r
    integer     , save     :: i = 0
    i = i + 1
    ! generate type indices
    tlist(i) % m_idx_gen = i
    ! point pointer to the data
    r => tlist(i)
    return
  end function tlist_new

  !---------------------------------------------------------------------------  
  ! (private)
  ! DESCRIPTION
  !> @brief add a new background type (by reference) to the tlist
  !> @remark intert background type as index zero
  !---------------------------------------------------------------------------  
  subroutine tlist_build_background ()
    ! setup background type and insert it into tlist
    ! @remark probably you should not change anything here !
    call tlist (0) % set_symm (1)      ! it has to be 1
    call tlist (0) % set_eva_num (0)   ! it can never be evaporated
    call tlist (0) % set_idx_def (0)   ! zero is the reserved index
    call tlist (0) % set_idx_gen (0)   ! zero is the reserved index
    call tlist (0) % set_idx_off (0)   ! zero is the reserved index
    call tlist (0) % alloc_comp  (1)   ! one component only
    call tlist (0) % alloc_reac  (0)   ! no reaction (no movement)
    tlist (0) % m_comp(:,1) = [0,0]    ! the only component is at origin
    call tlist(0) % check()
    return
  end subroutine tlist_build_background

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Rotate vector v for N times using the rotational matrix
  !> @param v: the input vector
  !> @param n: number of rotations
  !---------------------------------------------------------------------------  
  function m_rotate (this, v, n) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: v(2), n
    integer                   :: r(2), np
    np = modulo(n, this % symm)
    r = rotate (this % m_rmat, v, np)
    return
  end function m_rotate
  
  function m_translate (this, v, m) result (r)
    class(mtype), intent (in) :: this
    integer     , intent (in) :: v(2), m(3)
    integer                   :: r(2)
    r = this % rotate(v, m(3)) + m(1:2)
    return
  end function m_translate

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief function to convert relative index (index with respect to the same
  !         type) to absolute index (index with respect to all molecules)
  !---------------------------------------------------------------------------  
  function m_abs_idx (this, i) result (r)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: i
    integer                      :: r
    r = i + this % m_idx_off
    return
  end function m_abs_idx

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief function to convert absolute index (index with respect to all 
  !         molecules) to relative index (index with respect to the same
  !         type)
  !---------------------------------------------------------------------------  
  function m_rel_idx (this, i) result (r)
    class(mtype), intent (inout) :: this
    integer     , intent (in)    :: i
    integer                      :: r
    r = i - this % m_idx_off
    return
  end function m_rel_idx

end module class_mtype
