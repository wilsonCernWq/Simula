program main

  !use FoX_dom
  use global
  use func_substrate
  use func_rate_kmc
  use func_xml_reader
  use define
  implicit none

  integer :: i
  integer, parameter :: HOP_STP=1E2
  integer, parameter :: OUT_STP=HOP_STP/10

  call xml_read()

! #ifdef _WIN32
!   call set_root_dir("..\..\..\output")
! #else
!   call set_root_dir("/home/qiwu/work/dev/simula/fortran/output")
! #endif

!   call set_proj_dir("test-qiwu")

!   print *, ">>> Simula"
!   print *, ">>> Initialization"
  
!   call init()
!   print *, ">>> Evaporation"


!   call activate_new(1)
!   print *, land_one(1, 5, 1,0)
!   call activate_new(1)
!   print *, land_one(3, 5, 1,0)
!   call activate_new(1)
!   print *, land_one(5, 5, 1,0)


!   call evaporate(tpyp, 10)
!   call evaporate(lead, 10)
!   call print_to(6, 4)

!   print *, ">>> Calculate rate"
!   print *, ""
!   do i = 1, HOP_STP
!      call compute_rates(verbose = .false.)
!      if (modulo(i, OUT_STP)==0) then
!         call start_file(90)
!         call print_to(90,4)
!         call close_file(90)
!         print *, "total time", tot_time
!      end if
!   end do
!   print *, ""
!   call print_to(6, 4)
  
  stop
end program main
