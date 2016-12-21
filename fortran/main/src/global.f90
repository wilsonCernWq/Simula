!-----------------------------------------------------------------------------
!
! DESCRIPTION
!> @brief definitions for global variables and their setters
!
!-----------------------------------------------------------------------------
module global

  use func_helper
  implicit none
  public

  !---------------------------------------------------------------------------
  !> @var it defines the substrate symmetry, which will be related to 
  !       rotational generation and substrate indexing
  integer , save :: SYS_SYMM = 4
  real(dp), save :: tot_time = 0.0_dp
  real(dp), save :: env_temp = 400.0_dp

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  !> @brief system symmetry setter
  !---------------------------------------------------------------------------
  subroutine set_sys_symm(s)
    integer, intent (in) :: s
    SYS_SYMM = s
    return
  end subroutine set_sys_symm

  function rate_compute(e) result(r)
    real(dp) :: e, r
    real(dp), parameter :: KB = 8.625E-5
    r = 10E12  * exp(e / KB / env_temp)   
    return
  end function

end module global
