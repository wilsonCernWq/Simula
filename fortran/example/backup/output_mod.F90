!------------------------------------------------------------------------------
! Department of Physics, The Hong Kong University of Science and Technology
!------------------------------------------------------------------------------
!
! Module: output_mod
!
!> @author
!> Wilson
!
! DESCRIPTION: 
!> Construct component
!
!------------------------------------------------------------------------------
module output_mod
use glob_mod
use rcd_mod
implicit none

character(len=20), save :: folder = ""
character(len=10), save :: fbase = ""

contains
!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief Subroutine to set base folder
!> @param name: base folder name
!> @return none
!--------------------------------------------------------------------------- 
subroutine set_base_folder(name)
	character*(*) :: name
		fbase = trim(adjustl(name))
	return
end subroutine set_base_folder
!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief Subroutine to set current folder name
!> @param name: folder name index
!> @return none
!--------------------------------------------------------------------------- 
subroutine new_out_folder(name)
	integer :: name(:), i
	character(len=4) :: part
  folder = ""
	do i = 1, size(name)
		write(part, '(I4)') name(i)
		folder = trim(adjustl(folder))//trim(adjustl(part))
		if (i /= size(name)) folder = trim(adjustl(folder))//"-"
	end do
	call system('mkdir .\'//trim(adjustl(fbase))//'\'//trim(adjustl(folder)))
	return
end subroutine new_out_folder
!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> @brief Subroutine to print substrate
!> @param sub: substrate
!> @return none
!--------------------------------------------------------------------------- 
subroutine print_sub()
  integer :: i, j
  integer :: islandp(SIDE,SIDE)
  integer, save :: out_log = 0
  character (len=4) :: fname !< file name
	character (len=4) :: fsize !< matrix size
		!% write information
		write(fname,'(I4)') out_log
		write(fsize,'(I4)') SIDE
    out_log = out_log + 1
		! open file
		open(unit=90, file='.\'//trim(adjustl(fbase))//'\'//trim(adjustl(folder))//'\'//trim(adjustl(fname))//'.txt')
		! copy substrate into integer array
		islandp = 0
		forall (i=1:SIDE:1, j=1:SIDE:1)
			islandp(i,j) = rcd_name(get_sub(i,j))
		end forall
		!\ print substrate
		do i = 1, SIDE
			write(90, '('//trim(adjustl(fsize))//'I6)') islandp(i, 1:SIDE)
		end do
		close(unit=90, status='KEEP')
	return
end subroutine print_sub

end module output_mod