module mono_mod
use glob_mod
implicit none

type::mono
  integer, allocatable, dimension(:,:) :: str
  integer, allocatable, dimension(:,:) :: sur
  integer, allocatable, dimension(:,:) :: bor
  integer :: flu = 0
  integer :: num = 0
end type

integer, save :: no_of_comp = 0
integer, save :: no_of_mono = 0
integer, save :: max_sur = 0
integer, allocatable, save :: max_bor(:)
integer, allocatable, save :: eva_log(:) !< evaporation log
type(mono), allocatable, save :: comps(:)

contains

subroutine mono_init()
  integer :: no, i
  !=> overall initialization
  no_of_comp = 2 !< no of components
  if (allocated(comps)) deallocate(comps)
  allocate(comps(no_of_comp))
  if (allocated(eva_log)) deallocate(eva_log)
  allocate(eva_log(no_of_comp))
  if (allocated(max_bor)) deallocate(max_bor)
  allocate(max_bor(no_of_comp))
  eva_log = 0
  !=> comp --> TPYP
  no = 1
  allocate(comps(no)%str(9, 2)) !=> structure
  comps(no)%str = transpose(reshape([0,0, 1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1], [2,9]))
  allocate(comps(no)%sur(8, 2)) !=> surrounding points
  comps(no)%sur = transpose(reshape([2,0, 0,2, -2,0, 0,-2, 2,2, -2,2, -2,-2, 2,-2], [2,8]))
  allocate(comps(no)%bor(16, 2)) !=> border points
  comps(no)%bor = transpose(reshape([2,0, 0,2, -2,0, 0,-2, 2,1, 2,-1, -1,2, 1,2, -2,1, -2,-1, -1,-2, 1,-2, 2,2, -2,2, -2,-2, 2,-2], [2,16]))
  max_bor(no) = 16
  comps(no)%flu = 5
  comps(no)%num = 50
  !=> comp --> LEAD
  no = 2
  allocate(comps(no)%str(1, 2)) !=> structure
  comps(no)%str = transpose(reshape([0,0], [2,1]))
  allocate(comps(no)%sur(8, 2)) !=> surrounding points
  comps(no)%sur = transpose(reshape([1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1], [2,8]))
  allocate(comps(no)%bor(8, 2)) !=> border points
  comps(no)%bor = transpose(reshape([1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1], [2,8]))
  max_bor(no) = 8
  comps(no)%flu = 20
  comps(no)%num = 100
  !=> total monomer number
  no_of_mono = 0
  !=> max surrounding number
  max_sur = 8
  return
end subroutine mono_init

!> @name rotation function
!> @{
elemental function rotate_x (x, y, d) result (rx)
  integer, intent(in) :: x, y, d
  integer :: rd, rx
    rd = modulo(d, 4)
    select case (rd)
    case (0)
      rx =  x
    case (1)
      rx = -y
    case (2)
      rx = -x
    case (3)
      rx =  y 
    end select
  return
end function rotate_x

elemental function rotate_y (x, y, d) result (ry)
  integer, intent(in) :: x, y, d
  integer :: rd, ry
    rd = modulo(d, 4)
    select case (rd)
    case (0)
      ry =  y
    case (1)
      ry =  x
    case (2)
      ry = -y
    case (3)
      ry = -x
    end select
  return
end function rotate_y
!> @}

!> @{
function test_occupied (pos, d, name, sid) result (stat)
  integer, intent(in) :: pos(2), name, d
  integer, optional :: sid !< self index
  logical :: stat
  integer :: i, rank
  integer, allocatable :: sta(:,:), cor(:,:), r(:)
    stat = .false.
    !=> copy structure
    allocate (sta, source = comps(name)%str)
    allocate (cor, source = sta)
    !=> shape of array
    rank = size(sta, 1) 
    !=> rotate + translation
    cor(:,1) = rotate_x(sta(:,1), sta(:,2), d) + pos(1)
    cor(:,2) = rotate_y(sta(:,1), sta(:,2), d) + pos(2)
    !=> check occupied
    allocate (r(rank), source = get_sub(cor(:,1), cor(:,2)))
    if (.not. present(sid)) then
      if (sum(r) /= 0) stat = .true.
    else
      do i = 1, rank
        if (r(i) /= sid .and. r(i) /= 0) stat = .true.
      end do
    end if
  return
end function test_occupied

subroutine land (pos, d, name, id)
  integer, intent(in) :: pos(2), d, name, id
  integer :: i, rank
  integer, allocatable :: sta(:,:), cor(:,:), r(:)
    !=> copy structure
    allocate (sta, source = comps(name)%str)
    allocate (cor, source = sta)
    rank = size(sta, 1) !< shape of array
    !=> rotate + translation
    cor(:,1) = rotate_x(sta(:,1), sta(:,2), d) + pos(1)
    cor(:,2) = rotate_y(sta(:,1), sta(:,2), d) + pos(2)
    !=> land new monomer
    do i = 1, rank
      call set_sub(cor(i,1), cor(i,2), id)
    end do
  return
end subroutine land

subroutine deland (pos, d, name, id)
  integer, intent(in) :: pos(2), d, name, id
  integer :: i, rank
  integer, allocatable :: sta(:,:), cor(:,:), r(:)
    !=> copy structure
    allocate (sta, source = comps(name)%str)
    allocate (cor, source = sta)
    rank = size(sta, 1) !< shape of array
    !=> rotate
    cor(:,1) = rotate_x(sta(:,1), sta(:,2), d)
    cor(:,2) = rotate_y(sta(:,1), sta(:,2), d)
    !=> translation
    cor(:,1) = cor(:,1) + pos(1)
    cor(:,2) = cor(:,2) + pos(2)
    !=> deland monomer
    do i = 1, rank
      call set_sub(cor(i,1), cor(i,2), -id)
    end do
  return
end subroutine deland
!> @}

subroutine get_sur (pos, d, name, r)
  integer, intent(in) :: pos(2), d, name
  integer, allocatable, intent(out) :: r(:,:) !< return array
  integer, allocatable :: sta(:,:)
  integer :: i, rank
    !=> copy surrounding
    allocate (sta, source = comps(name)%sur)
    rank = size(sta, 1) !< shape of array
    if (.not.allocated(r)) allocate (r, source = sta)
    !=> rotate + translation
    r(:,1) = rotate_x(sta(:,1), sta(:,2), d) + pos(1)
    r(:,2) = rotate_y(sta(:,1), sta(:,2), d) + pos(2)
  return
end subroutine get_sur

subroutine get_bor (pos, d, name, r)
  integer, intent(in) :: pos(2), d, name
  integer, allocatable, intent(out) :: r(:,:) !< return array
  integer, allocatable :: sta(:,:)
  integer :: i, rank
    !=> copy surrounding
    allocate (sta, source = comps(name)%bor)
    rank = size(sta, 1) !< shape of array
    if (.not.allocated(r)) allocate (r, source = sta)
    !=> rotate + translation
    r(:,1) = rotate_x(sta(:,1), sta(:,2), d) + pos(1)
    r(:,2) = rotate_y(sta(:,1), sta(:,2), d) + pos(2)
  return
end subroutine get_bor

end module mono_mod