program main
  !-------------------------------------------------
  !> Load Modules
  !use FoX_dom
  use global
  use func_substrate
  use func_rate_kmc
  use func_xml_reader
  use define
  implicit none
  !-------------------------------------------------
  !> Local Variables
  integer, parameter :: OUT_TIME=HOP_TIME/100
  integer :: time
  !-------------------------------------------------
  !> XML reader
  ! call xml_read()
  !-------------------------------------------------
  !> Set Pathes
#ifdef _WIN32
  call set_root_dir("..\..\..\output")
#else
  call set_root_dir("../output")
#endif
  call set_proj_dir("debug")
  !-------------------------------------------------
  !> Print Header
  print *, ">>> Simula"
  print *, ">>> Initialization"  
  !> Initialization
  call init()
  !-------------------------------------------------
  !> Evaluation
  print *, ">>> Evaporation"
  call evaporate(tpyp, 40)
  call evaporate(lead, 10)
  call print_to(6, 4)
  !-------------------------------------------------
  !> KMC
  print *, ">>> Calculate rate"
  print *, ""
  do time = 1, HOP_TIME
     call compute_rates(verbose = .false.)
     if (modulo(time, OUT_TIME)==0) then
        !> print to file
        call start_file(90)
        call print_to(90,4)
        call close_file(90)
        !> print to screen
        call print_to(6, 4)
        print *, "total time", tot_time
     end if
  end do
  print *, ""
  call print_to(6, 4)
  !-------------------------------------------------
  stop
end program main
