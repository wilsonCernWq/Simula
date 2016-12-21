!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!> Define all the helper functions independent to the simulation project
!
!---------------------------------------------------------------------------  
module func_helper

  implicit none

  !---------------------------------------------------------------------------  
  ! define double precision here
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
  integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))

  !---------------------------------------------------------------------------  
  ! function overload
  interface alloc
     module procedure &
          alloc_I1, alloc_I1_ptr, &
          alloc_I2, alloc_I2_ptr, &
          alloc_F1, alloc_F2, alloc_F1_range
  end interface alloc

  interface binary_search
     module procedure binary_search_F1, binary_search_I1
  end interface binary_search

contains
  
  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Function to generate random integer (include boundaries)
  !> @param[in] sint: lower bound
  !> @param[in] eint: upper bound
  !> @return integer random number
  !--------------------------------------------------------------------------- 
  function rand_int(sint, eint)
    integer,intent(IN) :: sint, eint
    integer :: rand_int
    real(8) :: r
    call random_number(r)
    rand_int = floor(r*(eint-sint+1))+sint
    return
  end function rand_int

  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Function to generate random float
  !> @param[in] spnt: lower bound
  !> @param[in] epnt: upper bound
  !> @return real random number
  !--------------------------------------------------------------------------- 
  function rand_uniform(spnt, epnt)
    real(8),intent(IN) :: spnt, epnt
    real(8) :: r, rand_uniform
    call random_number(r)
    rand_uniform = r*(epnt-spnt) + spnt
    return
  end function rand_uniform

  !---------------------------------------------------------------------------  
  ! Description:
  !> @reference https://gcc.gnu.org/onlinedocs/gfortran/RAND.html
  !
  !   Restarts or queries the state of the pseudorandom number generator used
  !   by RANDOM_NUMBER.
  !
  !   If RANDOM_SEED is called without arguments, it is initialized to a default
  !   state. The example below shows how to initialize the random seed based on
  !   the system's time.
  !
  ! Standard : F95 and later
  !
  ! Class    : Subroutine
  !
  ! Syntax   : CALL RANDOM_SEED(SIZE, PUT, GET)
  !
  ! Arguments:
  !     SIZE (Optional) Shall be a scalar and of type default INTEGER, with
  !                     INTENT(OUT). It specifies the minimum size of the arrays
  !                     used with the PUT and GET arguments.
  !     PUT (Optional)  Shall be an array of type default INTEGER and rank one.
  !                     It is INTENT(IN) and the size of the array must be 
  !                     larger than or equal to the number returned by the SIZE
  !                     argument.
  !     GET (Optional)  Shall be an array of type default INTEGER and rank one. 
  !                     It is INTENT(OUT) and the size of the array must be
  !                     larger than or equal to the number returned by the SIZE
  !                     argument.
  !---------------------------------------------------------------------------  
  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine init_random_seed

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief rotate vector once by rotational matrix
  !> @param mat: matrix
  !> @param vec: vector
  !
  ! MATRIX INDEX
  ! (1,1) (1,2)
  ! (2,1) (2,2)
  !---------------------------------------------------------------------------  
  function rotate_once(mat, vec)
    integer, intent(in) :: mat(2,2), vec(2)
    integer :: rotate_once(2)
    rotate_once(1) = vec(1) * mat(1,1) + vec(2) * mat(1,2)
    rotate_once(2) = vec(1) * mat(2,1) + vec(2) * mat(2,2)
    return
  end function rotate_once

  !---------------------------------------------------------------------------  
  ! DESCRIPTION
  !> @brief rotate vector n times by rotational matrix
  !> @param mat: matrix
  !> @param vec: vector
  !> @param n  : number of rotation
  !---------------------------------------------------------------------------  
  function rotate(mat, vec, n)
    integer, intent(in) :: mat(2,2), vec(2), n
    integer :: tmp(2), rotate(2), i
    tmp = vec
    do i = 1,n
       tmp = rotate_once(mat, tmp)
    end do
    rotate = tmp
    return
  end function rotate
  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief allocate allocatable array for basic types
  !> @param array: the allocatable target array
  !> @param {n}  : dimension specifications
  !---------------------------------------------------------------------------  
  subroutine alloc_I1 (array, n)
    integer, allocatable, intent(out) :: array(:)
    integer             , intent(in)  :: n
    integer :: status
    if (allocated(array)) stop "ERROR: (I1) multiple definitions"
    allocate(array(n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_I1

  subroutine alloc_I1_ptr (array, n)
    integer, pointer, intent(out) :: array(:)
    integer         , intent(in)  :: n
    integer :: status
    if (associated(array)) stop "ERROR: (I1 ptr) multiple definitions"
    allocate(array(n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_I1_ptr

  subroutine alloc_I2 (array, n1, n2)
    integer, allocatable, intent(out) :: array(:,:)
    integer             , intent(in)  :: n1, n2
    integer :: status
    if (allocated(array)) stop "ERROR: (I2) multiple definitions"
    allocate(array(n1,n2), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_I2

  subroutine alloc_I2_ptr (array, n1, n2)
    integer, pointer, intent(out) :: array(:,:)
    integer         , intent(in)  :: n1, n2
    integer :: status
    if (associated(array)) stop "ERROR: (I2 ptr) multiple definitions"
    allocate(array(n1,n2), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_I2_ptr

  subroutine alloc_F1 (array, n)
    real(8), allocatable, intent(out) :: array(:)
    integer             , intent(in)  :: n
    integer :: status
    if (allocated(array)) stop "ERROR: (F1) multiple definitions"
    allocate(array(n), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_F1

  subroutine alloc_F1_range (array, n1, n2)
    real(8), allocatable, intent(out) :: array(:)
    integer             , intent(in)  :: n1, n2
    integer :: status
    if (allocated(array)) stop "ERROR: (F1 range) multiple definitions"
    allocate(array(n1:n2), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_F1_range

  subroutine alloc_F2 (array, n1, n2)
    real(8), allocatable, intent(out) :: array(:,:)
    integer             , intent(in)  :: n1, n2
    integer :: status
    if (allocated(array)) stop "ERROR: (F2) multiple definitions"
    allocate(array(n1, n2), STAT = status)
    if (status /= 0) stop "ERROR: Not enough memory!"
    return
  end subroutine alloc_F2

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief binary search algorithm for finding bin index for a given number
  !> @param  array: bin edge array, which should have N+1 values
  !> @param  p    : number to select
  !> @return r    : integer within [1,N] so that array(r) < p < array(r+1)
  !---------------------------------------------------------------------------  
  function binary_search_F1(array, p) result (r)
    real(dp), intent (in) :: array (:)
    real(dp), intent (in) :: p
    integer               :: r
    integer :: nl, nr, nm
    nl = 1
    nr = size(array)
    nm = (nl + nr) / 2
    do while (nr - nl > 1) 
       if (p <= array(nm)) then
          nr = nm
       else
          nl = nm
       end if
       nm = (nl + nr) / 2
    end do
    r = nl
    return
  end function binary_search_F1

  function binary_search_I1(array, p) result (r)
    integer, intent (in) :: array (:)
    integer, intent (in) :: p
    integer              :: r
    integer :: nl, nr, nm
    nl = 1
    nr = size(array)
    nm = (nl + nr) / 2
    do while (nr - nl > 1) 
       if (p <= array(nm)) then
          nr = nm
       else
          nl = nm
       end if
       nm = (nl + nr) / 2
    end do
    r = nl
    return
  end function binary_search_I1

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief function to convert integer to string
  !> @param  i: input integer
  !---------------------------------------------------------------------------  
  function int2str(i) result(c)
    integer, intent(in)  :: i
    character (len = int(log(real(i))/log(10.0_dp))+1) :: c
    character (len = 2)                                :: f    
    write (f,'(I2)') int(log(real(i))/log(10.0_dp))+1
    write (c,'(I'//trim(f)//')') i
    return
  end function int2str

end module func_helper
