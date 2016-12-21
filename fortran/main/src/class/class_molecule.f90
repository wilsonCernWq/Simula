!-----------------------------------------------------------------------------  
!
! DESCRIPTION
!> This class defines the behaviors of one molecule
!
! FIELD
!
!-----------------------------------------------------------------------------  
module class_molecule

  !---------------------------------------------------------------------------  
  !> used module
  use func_helper, only: alloc
  use class_mtype, only: mtype, tlist
  implicit none
  private
  
  !---------------------------------------------------------------------------  
  !> each molecule 
  type, public :: molecule
     integer, public , pointer :: tp
     integer, public , pointer :: id
     integer, public , pointer ::   pos   (:) => null()
     integer, private          :: m_pos_num
     integer, public , pointer ::   state (:) => null()
     integer, private          :: m_state_num
     integer, private, pointer :: m_data(:)   => null()
     ! => [t,i,x,y,d,s1,s2,s3,...]
   contains
     procedure :: init      => m_init
     procedure :: pos_num   => m_get_pos_num
     procedure :: state_num => m_get_state_num
  end type molecule

  !---------------------------------------------------------------------------  
  !> Global variables 
  type (molecule), public, allocatable, target :: mlist (:) ! molecule list

  !---------------------------------------------------------------------------  
  !> global functions
  public :: mlist_init, mlist_num

contains

  subroutine m_init (this, n)
    class(molecule), intent (inout) :: this
    integer        , intent (in)    :: n
    call alloc (this % m_data, 5+n)
    this % tp  => this % m_data(1)
    this % id  => this % m_data(2)
    this % pos   => this % m_data(3:5)
    this % state => this % m_data(6:5+n)
    return
  end subroutine m_init

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for state list size
  !---------------------------------------------------------------------------  
  function m_get_state_num (this) result (r)
    class(molecule), intent (in) :: this
    integer                      :: r
    r = this % m_state_num
    return
  end function m_get_state_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for pos list size
  !---------------------------------------------------------------------------  
  function m_get_pos_num (this) result (r)
    class(molecule), intent (in) :: this
    integer                      :: r
    r = this % m_pos_num
    return
  end function m_get_pos_num

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  !> @brief Getter for tlist size
  !---------------------------------------------------------------------------
  function mlist_num () result (r)
    integer :: r
    r = size(mlist)-1
    return
  end function mlist_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief initialize molecule list (mlist)
  !---------------------------------------------------------------------------  
  subroutine mlist_init()
    integer :: n         ! total number of molecules
    integer :: s         ! index step sum 
    integer :: t         ! molecule type
    integer :: i, status

    type(mtype)   , pointer :: t_obj
    type(molecule), pointer :: m_obj

    integer, pointer :: p(:)

    ! calculate total number of molecules
    n = 0
    EACH_TYPE: do t = 1, size(tlist)-1
       n = n + tlist(t) % eva_num()
    end do EACH_TYPE

    ! Printing total number of molecules
    write (*,'(" No. of molecules", I6)') n

    ! not a basic type, allocate manually
    ! @remark id = 0 is researved for background type
    if (allocated(mlist)) stop "ERROR: (alloc mlist) multiple definitions"
    allocate (mlist (0:n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"

    ! initialize molecule
    t = 0
    s = 0 
    EACH_MOLECULE: do i = 0, n
       ! @remark update current molecule type. this check should always be
       !         done before other things. it will ignore the first case
       if (i > s) then 
          t = t + 1
          call tlist(t) % set_idx_off (s) ! index offset for molecule
          s = s + tlist(t) % eva_num ()   ! set next index step
       end if
       ! allocate state number
       t_obj => tlist(t)
       m_obj => mlist(i)
       call m_obj % init( t_obj % comp_num() )
       ! set molecule properties
       m_obj % tp  = t
       m_obj % id  = i
       m_obj % pos = 0
       if (i == 0) then
          m_obj % state = 0 ! initial state for background
       else          
          m_obj % state =  t_obj % state_initial() ! initial state
       end if
    end do EACH_MOLECULE

    return
  end subroutine mlist_init

end module class_molecule
