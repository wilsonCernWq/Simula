!-----------------------------------------------------------------------------
!
! DESCRIPTION
!> @brief This class defines the behaviors of one single general reaction. The
!         word general means it contains not only chemical bondings, but also
!         all free movements, bond breakings, etc. All those behaviors should
!         be explicitly defined either by hard coding or a NAMELIST.
!
! FIELD
!> @var idx: reaction index
!> @var ene: reaction energy
!> @var mov: reaction moving (rotation) direction / specification
!> @var cond_num: length of condition list
!> @var conds   : condition list, the reaction can be executed if and only if 
!                 all the conditions listed here are passed
!
!-----------------------------------------------------------------------------
module class_reaction

  ! used modules
  use func_helper     , only: dp
  use class_condition , only: condition
  implicit none
  private

  ! descriotion for one reaction
  type, public :: reaction
     real(dp), public  :: energy   ! reaction energy
     integer , public  :: rid      ! reaction id for further reference
     integer , public  :: move (3) ! action specification [x, y, d]
     type (condition), public, pointer :: cond (:) => null()
     integer         , private         :: m_cond_num
   contains
     procedure :: set_energy => m_set_energy
     procedure :: set_rid    => m_set_rid
     procedure :: set_move   => m_set_move
     procedure :: alloc_cond => m_alloc_cond
     procedure :: cond_num   => m_get_cond_num
  end type reaction

contains

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for reaction index
  !> @param i : reaction index (automatically generated)
  !---------------------------------------------------------------------------  
  subroutine m_set_rid (this, i)
    class(reaction), intent (inout) :: this
    integer        , intent (in)    :: i
    this % rid  = i
    return
  end subroutine m_set_rid

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for reaction energy
  !> @param e : energy value 
  !             since the default fortran real precision is short precision, a
  !             type conversion from reak(4) to real(dp) is performed here
  !---------------------------------------------------------------------------  
  subroutine m_set_energy (this, e)
    class(reaction), intent (inout) :: this
    real(dp)       , intent (in)    :: e 
    this % energy = e
    return
  end subroutine m_set_energy

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for reaction action
  !> @param m : moving specification
  !---------------------------------------------------------------------------  
  subroutine m_set_move (this, m)
    class(reaction), intent (inout) :: this
    integer        , intent (in)    :: m (3)
    this % move = m
    return
  end subroutine m_set_move
  
  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Allocator for cond list
  !> @param n : number of conditions
  !---------------------------------------------------------------------------  
  subroutine m_alloc_cond (this, n)
    class(reaction), intent (inout) :: this
    integer        , intent (in)    :: n
    integer                         :: status
    ! allocation check
    if (associated(this % cond)) stop "ERROR: (alloc cond) multiple definitions"
    ! assign number value
    this % m_cond_num = n
    ! not a basic type, allocate manually
    allocate(this % cond(n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory! (class reaction)"
    return
  end subroutine m_alloc_cond

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for condition number
  !---------------------------------------------------------------------------  
  function m_get_cond_num (this) result (r)
    class(reaction), intent (inout) :: this
    integer :: r
    r = this % m_cond_num
    return
  end function m_get_cond_num
  
end module class_reaction
