program main
use sim_mod
use output_mod
  real(8) :: t1, t2

  call random_seed()
  call set_base_folder("Output")

  call cpu_time(t1)
  call sim_single(TEMPERATURE)
  call print_sub()
  call cpu_time(t2)
  print *, t2 - t1

  stop
end program