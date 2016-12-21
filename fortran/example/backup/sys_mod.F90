!------------------------------------------------------------------------------
! Department of Physics, The Hong Kong University of Science and Technology
!------------------------------------------------------------------------------
!
! Module: sys_mod
!
!> @author
!> Wilson
!
! DESCRIPTION: 
!> Module for system error output control
!
!------------------------------------------------------------------------------
module sys_mod
implicit none
contains
!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief Subroutine to print out error and stop the program
!> @param outputs: output string indicating the error type
!> @return None
!--------------------------------------------------------------------------- 
subroutine error(outputs)
	implicit none
	character*(*), optional, intent (IN) :: outputs
		if(present(outputs))then
			write(*,*) "ERROR: "//outputs
		else
			write(*,*) "ERROR"
		end if
	stop
end subroutine error

end module sys_mod