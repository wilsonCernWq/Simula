!------------------------------------------------------------------------------
! Department of Physics, The Hong Kong University of Science and Technology
!------------------------------------------------------------------------------
!
! Module: rate_mod
!
!> @author
!> Wilson
!
! DESCRIPTION: 
!> Construct rate list
!
!------------------------------------------------------------------------------
module rate_mod
use rand_mod
use mono_mod
implicit none

real(8), allocatable, save :: rates (:)
integer, allocatable, save :: infos (:,:)
real(8), save :: rate_sum = 0.0
integer, save :: rate_dat_sz = 6
integer, save :: rate_len = 0
integer, save :: rate_wid = 0
contains

subroutine rate_init()
  integer :: ntot, i
  !=> compute total number of monomers
  ntot = 0
  do i = 1, no_of_comp
    ntot = ntot + comps(i)%num
  end do
  !=> no of monte carlo states for each monomer
  rate_wid = (5 * SYMM-1) + 8
  rate_sum = 0.0
  rate_len = rate_wid * ntot
  !=> initialize rates list
  if (allocated(rates)) deallocate(rates)
  allocate(rates(rate_wid * ntot))
  rates = 0.0
  !=> initialize infos list
  if (allocated(infos)) deallocate(infos)
  allocate(infos(rate_dat_sz, rate_wid*ntot))
  infos = 0
  return
end subroutine rate_init

subroutine rate_set(n, v, info)
  integer, intent(in) :: n, info(rate_dat_sz)
  real(8), intent(in) :: v
    !=> update sum
    rate_sum = rate_sum - rates(n) + v
    !=> set value
    rates(n) = v
    infos(:,n) = info
  return
end subroutine rate_set

function rate_select() result (info)
  integer :: i, info(rate_dat_sz)
  real(8) :: re, rs
    rs = 0.0 !=> current sum
    !=> target rate
    re = rand_uniform(real(0.0,8), rate_sum)
    do i = 1, rate_len
      rs = rs + rates(i)
      if (rs >= re) exit
    end do
    !=> load information
    info = infos(:,i)
  return
end function rate_select

function rate_update_time(in_time) result (out_time)
  real(8), intent(in) :: in_time
  real(8) :: out_time, r
    call random_number(r)
    out_time = in_time + (log(1.0/r)/rate_sum)
  return
end function rate_update_time

end module rate_mod