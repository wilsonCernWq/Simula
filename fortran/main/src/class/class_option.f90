!-----------------------------------------------------------------------------
! 
! DESCRIPTION
!> @brief This is a class for store detail information for condition checkng
!
! FIELDS
!> @var state: it defines the required initial states and final states for
!              executing the reaction. the value pairs are stored following 
!              the component defining order
!
!> @var pos: relative xyd-coordinate list for targets with respect to current
!            object. It requires all the points should be filled with identical
!            objects, with identical corresponding relative direction. The 
!            direction value should follow the molecule symmetry definition, 
!            which means d < symm
!
!-----------------------------------------------------------------------------
module class_option

  ! used modules
  use func_helper, only: alloc
  implicit none
  private

  ! option information for condition
  type, public :: option
     integer, private, pointer :: m_state (:,:) => null() ! (2,M)
     integer, private          :: m_state_num             ! number M
     integer, private, pointer :: m_pos (:,:)   => null() ! (3,N) {[x,y,d]...}
     integer, private          :: m_pos_num               ! number N
   contains
     procedure :: pos_eq_to   => m_pos_eq_to
     procedure :: state_eq_to => m_state_eq_to
     ! getters
     procedure :: pos     => m_get_pos
     procedure :: pos_num => m_get_pos_num
     procedure :: x   => m_get_x
     procedure :: y   => m_get_y
     procedure :: d   => m_get_d
     procedure :: xy  => m_get_xy
     procedure :: yx  => m_get_yx
     procedure :: xd  => m_get_xd
     procedure :: dx  => m_get_dx
     procedure :: yd  => m_get_yd
     procedure :: dy  => m_get_dy
     procedure :: xyd => m_get_xyd
     procedure :: xdy => m_get_xdy
     procedure :: yxd => m_get_yxd
     procedure :: ydx => m_get_ydx
     procedure :: dxy => m_get_dxy
     procedure :: dyx => m_get_dyx
     procedure :: state_num  => m_get_state_num
     procedure :: state_vec  => m_get_state_vec
     procedure :: state_iptr => m_get_state_iptr
     procedure :: state_fptr => m_get_state_fptr
     ! setters
     procedure :: set_all   => m_set_all
     procedure :: set_pos   => m_set_pos
     procedure :: set_state => m_set_state
  end type option

contains

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  ! Setters
  !---------------------------------------------------------------------------  
  !> @brief Setter for option checking position list, a (x,y,d)-coordinate
  !         list. if all the positions listed here are valid, the option is
  !         considered to be true
  !> @param p: values for position array in 1d form {x1,y1,d1, x2,y2,d2 ...}
  !---------------------------------------------------------------------------
  subroutine m_set_pos (this, p)
    class(option), intent (inout) :: this
    integer      , intent (in)    :: p(:)
    integer                       :: n, i
    ! allocate array
    n = size (p) / 3
    call alloc (this % m_pos, 3, n)
    ! assign array
    do i = 1, n
       this % m_pos (1,i) = p (i * 3 - 2)
       this % m_pos (2,i) = p (i * 3 - 1)
       this % m_pos (3,i) = p (i * 3    )
    end do
    ! assign array size
    this % m_pos_num = n
    return
  end subroutine m_set_pos

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Setter for condition state array. All positions listed in position
  !         array should fulfill the same state specification
  !> @param s: values for state array in 1D form of {state_i,state_f, ...}
  !---------------------------------------------------------------------------  
  subroutine m_set_state (this, s)
    class(option), intent (inout) :: this
    integer      , intent (in)    :: s(:)
    integer                       :: n, i
    ! allocate array
    n = size(s) / 2
    call alloc (this % m_state, 2, n)
    ! assign matrix
    do i = 1, n
       this % m_state (1,i) = s (2 * i - 1)
       this % m_state (2,i) = s (2 * i    )
    end do
    ! assign array size
    this % m_state_num = n
    return
  end subroutine m_set_state

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief General setter for position and states
  !> @param p: same as set_pos
  !> @param d: same as set_state
  !---------------------------------------------------------------------------  
  subroutine m_set_all (this, p, s)
    class(option), intent (inout) :: this
    integer      , intent (in)    :: p(:)
    integer      , intent (in)    :: s(:)
    call this % set_pos   (p)
    call this % set_state (s)
    return
  end subroutine m_set_all

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  ! Getters
  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for i-th position vector
  !> @param i: index
  !---------------------------------------------------------------------------  
  function m_get_pos (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r = this % m_pos (:,i)
    return
  end function m_get_pos

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for position array size
  !---------------------------------------------------------------------------  
  function m_get_pos_num (this) result (r)
    class(option), intent (in) :: this
    integer                    :: r
    r = this % m_pos_num
    return
  end function m_get_pos_num

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for i-th position with different forms 
  !> @param i: index
  !---------------------------------------------------------------------------  
  function m_get_x (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r
    r = this % m_pos (1,i)
    return
  end function m_get_x

  function m_get_y (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r
    r = this % m_pos (2,i)
    return
  end function m_get_y

  function m_get_d (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r
    r = this % m_pos (3,i)
    return
  end function m_get_d

  function m_get_xy (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r = this % m_pos (1:2,i)
    return
  end function m_get_xy

  function m_get_yx (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r(1) = this % m_pos (2,i)
    r(2) = this % m_pos (1,i)
    return
  end function m_get_yx

  function m_get_xd (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r(1) = this % m_pos (1,i)
    r(2) = this % m_pos (3,i)
    return
  end function m_get_xd

  function m_get_dx (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r(1) = this % m_pos (3,i)
    r(2) = this % m_pos (1,i)
    return
  end function m_get_dx

  function m_get_yd (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r(1) = this % m_pos (2,i)
    r(2) = this % m_pos (3,i)
    return
  end function m_get_yd

  function m_get_dy (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r(1) = this % m_pos (3,i)
    r(2) = this % m_pos (2,i)
    return
  end function m_get_dy

  function m_get_xyd (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r = this % m_pos (:,i)
    return
  end function m_get_xyd

  function m_get_xdy (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r(1) = this % m_pos (1,i)
    r(2) = this % m_pos (3,i)
    r(3) = this % m_pos (2,i)
    return
  end function m_get_xdy

  function m_get_yxd (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r(1) = this % m_pos (2,i)
    r(2) = this % m_pos (1,i)
    r(3) = this % m_pos (3,i)
    return
  end function m_get_yxd

  function m_get_ydx (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r(1) = this % m_pos (2,i)
    r(2) = this % m_pos (3,i)
    r(3) = this % m_pos (1,i)
    return
  end function m_get_ydx

  function m_get_dxy (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r(1) = this % m_pos (3,i)
    r(2) = this % m_pos (1,i)
    r(3) = this % m_pos (2,i)
    return
  end function m_get_dxy

  function m_get_dyx (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(3)
    r(1) = this % m_pos (3,i)
    r(2) = this % m_pos (2,i)
    r(3) = this % m_pos (1,i)
    return
  end function m_get_dyx

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for option state array i-th vector
  !> @param i: index
  !---------------------------------------------------------------------------  
  function m_get_state_vec (this, i) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i
    integer                    :: r(2)
    r = this % m_state(:, i)
    return
  end function m_get_state_vec

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for option state array's all initial states
  !> @return array pointer to the original data !!! be careful !!!
  !---------------------------------------------------------------------------  
  function m_get_state_iptr (this) result (r)
    class(option), intent (in) :: this
    integer, pointer           :: r(:)
    r => this % m_state(1, :)
    return
  end function m_get_state_iptr

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for option state array's all final states
  !> @return array pointer to the original data !!! be careful !!!
  !---------------------------------------------------------------------------  
  function m_get_state_fptr (this) result (r)
    class(option), intent (in) :: this
    integer, pointer           :: r(:)
    r => this % m_state(2, :)
    return
  end function m_get_state_fptr

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief Getter for condition state array size
  !---------------------------------------------------------------------------  
  function m_get_state_num (this) result (r)
    class(option), intent (in) :: this
    integer                    :: r
    r = this % m_state_num
    return
  end function m_get_state_num

  function m_pos_eq_to (this, i, pos) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: i, pos(3)
    logical                    :: r
    r = all(this % m_pos(:,i) == pos)
    return
  end function m_pos_eq_to

  function m_state_eq_to (this, s) result (r)
    class(option), intent (in) :: this
    integer      , intent (in) :: s(:)
    logical                    :: r
    r = all(this % m_state(1,:) == s)
    return
  end function m_state_eq_to

end module class_option
