!-----------------------------------------------------------------------------
! 
! DESCRIPTION
!> @brief This is a class for dealing with reaction condition checking
!
! FIELDS
!> @var mtp: target type index (generated index)
!
!> @var opt: option list. it defines all the information for checking in 
!            details. the condition is determined to be ture if at least one
!            of the options is observed to be true
!
!-----------------------------------------------------------------------------
module class_condition

  ! used modules
  use func_helper , only: alloc
  use class_option, only: option
  implicit none
  private

  ! changes for one reactant within one reaction (bonding)  
  type, public :: condition
     type(option), public , pointer :: opt (:) => null()
     integer     , private          :: m_opt_num
     integer     , private          :: m_tp
   contains
     ! functions for type index
     procedure :: tp_eq_to  => m_tp_eq_to
     procedure :: tp        => m_get_tp
     procedure :: set_tp    => m_set_tp
     ! functions for option list
     procedure :: opt_num   => m_get_opt_num
     procedure :: alloc_opt => m_alloc_opt
  end type condition

contains

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for condition target type
  !> @param t: target type of the condition instance
  !---------------------------------------------------------------------------  
  subroutine m_set_tp (this, t)
    class(condition), intent (inout) :: this
    integer         , intent (in)    :: t
    this % m_tp = t
    return
  end subroutine m_set_tp

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for condition target type
  !---------------------------------------------------------------------------  
  function m_get_tp (this) result (r)
    class(condition), intent (in) :: this
    integer :: r
    r = this % m_tp
    return
  end function m_get_tp

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Comparator for target type
  !> @param t: target type
  !---------------------------------------------------------------------------  
  function m_tp_eq_to (this, t) result (r)
    class(condition), intent (in) :: this
    integer         , intent (in) :: t
    logical                       :: r
    r = (this % m_tp == t)
    return
  end function m_tp_eq_to

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for option number
  !---------------------------------------------------------------------------  
  function m_get_opt_num (this) result (r)
    class(condition), intent (in) :: this
    integer                       :: r
    r = this % m_opt_num
    return
  end function m_get_opt_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Allocator for option list
  !> @param n: number of options 
  !---------------------------------------------------------------------------  
  subroutine m_alloc_opt (this, n)
    class(condition), intent (inout) :: this
    integer         , intent (in)    :: n
    integer                          :: status
    ! check allocation
    if (associated(this % opt)) stop "ERROR: (alloc opt) multiple definitions"
    ! allocate memory space
    allocate (this % opt (n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory"
    ! assign array size
    this % m_opt_num = n
    return
  end subroutine m_alloc_opt

end module class_condition
