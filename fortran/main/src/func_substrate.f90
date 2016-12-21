!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!> @brief Define actions related to substrate
!
!> @function init_substrate
!> @function set_sub
!> @function get_sub
!> ...
!
!--------------------------------------------------------------------------- 
#include "macro_slash.f90"
module func_substrate

  use func_helper   , only: alloc, rand_int
  use class_mtype   , only: tlist, tlist_num, mtype
  use class_molecule, only: mlist, mlist_num, molecule
  implicit none

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief substrate static variable
  integer, private :: m_xsize
  integer, private :: m_ysize
  integer, private, allocatable :: m_sub (:,:) ! substrate data
  integer, private, allocatable :: m_num (:)   ! number of molecule activated

  character(len=200), save :: m_root_dir = "."
  character(len=100), save :: m_curr_dir = "."

contains
  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Subroutine to set base folder
  !> @param name: base folder name
  !> @return none
  !--------------------------------------------------------------------------- 
  subroutine set_root_dir(name)
    character*(*) :: name
    m_root_dir = trim(adjustl(name))
    write (*, '(" Setting output root directory: ", A50)') m_root_dir
    return
  end subroutine set_root_dir

  subroutine set_proj_dir(name)
    character*(*) :: name
    character (len=8)  :: dates
    character (len=10) :: times
    ! write date information into project folder name
    call date_and_time(DATE=dates,TIME=times)
    ! create folder name and make folder
    m_curr_dir = trim(adjustl(name))//"-"//dates//"-"//times(1:6)
    write (*, '(" Setting project output directory: ", A50)') m_curr_dir
    call system('mkdir '//               &
         trim(adjustl(m_root_dir))//macro_SLASH//&
         trim(adjustl(m_curr_dir)))
    return
  end subroutine set_proj_dir

  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Subroutine to open and close file
  !> @param u: file unit idx
  !--------------------------------------------------------------------------- 
  subroutine start_file(u)
    integer  , intent(in)   :: u
    integer  , save         :: file_idx  = 0
    integer  , dimension(8) :: dates
    integer                 :: i
    character (len=4 )      :: fpart(8)
    character (len=99)      :: fname    ! file name
    ! write information
    call date_and_time(VALUES=dates)
    ! construct unique filename
    do i = 1, 8
       write (fpart(i), '(I4)') dates(i)
    end do
    write(fname,'(I4,"-")') file_idx
    do i = 1, 8
       fname = trim(adjustl(fname))//trim(adjustl(fpart(i)))
    end do
    fname = &
         trim(adjustl(m_root_dir))//macro_SLASH// &
         trim(adjustl(m_curr_dir))//macro_SLASH// &
         trim(adjustl(fname))//'.txt'
    write (*, '(" Saving data into ",A99)') fname
    ! openfile
    open(u, file=trim(adjustl(fname)))
    file_idx  = file_idx + 1
    return
  end subroutine start_file

  subroutine close_file (u)
    integer, intent(in) :: u
    close(u, status='KEEP')
    return
  end subroutine close_file

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief Get number of activated molecules of type i
  !> @param i: ith molecule type
  !--------------------------------------------------------------------------- 
  function activated_num (i) result (r)
    integer, intent (in) :: i
    integer              :: r
    if (i > tlist_num()) stop "ERROR: index out of bound (module-substrate)"
    r = m_num(i)
    return
  end function activated_num

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief Increase activated molecules of type i by one
  !> @param i: ith molecule type
  !--------------------------------------------------------------------------- 
  subroutine activate_new (i)
    integer, intent (in) :: i
    if (i > tlist_num()) stop "ERROR: index out of bound (module-substrate)"
    m_num (i) = m_num (i) + 1
    return
  end subroutine activate_new

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief Allocate and initialize substrate
  !> @param xlen
  !> @param ylen 
  !--------------------------------------------------------------------------- 
  subroutine init_substrate (xlen, ylen)
    integer, intent (in) :: xlen, ylen
    m_xsize = xlen
    m_ysize = ylen
    call alloc (m_sub, xlen, ylen)
    call alloc (m_num, tlist_num())
    m_sub = 0
    m_num = 0
    return
  end subroutine init_substrate

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief generate a random X coordinate
  !--------------------------------------------------------------------------- 
  function rand_subX ()
    integer :: rand_subX
    rand_subX = rand_int(1, m_xsize)
    return
  end function rand_subX

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief generate a random Y coordinate
  !--------------------------------------------------------------------------- 
  function rand_subY ()
    integer :: rand_subY
    rand_subY = rand_int(1, m_ysize)
    return
  end function rand_subY

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief compress dot information into one number
  !> @param mid  : molecule index
  !> @param tid  : type index
  !> @param comp : component
  !> @param state: dot state
  !--------------------------------------------------------------------------- 
  function convert_to_land (mid, tid, comp, state)
    integer :: mid, tid, comp, state, convert_to_land
    !> check if numbers are all valid
    if (tid   > 99) stop "ERROR: number of types cannot be larger than 99"
    if (comp  > 99) stop "ERROR: number of dots cannot be larger than 99"
    if (state > 99) stop "ERROR: number of states cannot be larger than 99"
    !> compute values
    convert_to_land = state + 100 * comp + 10000 * tid + 1000000 * mid 
    return
  end function convert_to_land

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief convert substrate value back to readable number
  !> @param v: substrate value
  !> @param k: 1 => molecule id, 2 => type id, 3 => component, 4 => state
  !--------------------------------------------------------------------------- 
  function convert_from_land(v, k)
    integer :: v, k, convert_from_land
    if (k == 1) then
       convert_from_land = v / 1000000
    else if (k == 2) then
       convert_from_land = mod(v / 10000, 100)
    else if (k == 3) then
       convert_from_land = mod(v / 100, 100)
    else if (k == 4) then
       convert_from_land = mod(v, 100)
    end if
    return
  end function convert_from_land

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief Access substrate point at (x,y)
  !> @param m_sub: substrate structure
  !> @param x: x coordinate
  !> @param y: y coordinate
  !--------------------------------------------------------------------------- 
  elemental function get_sub (x, y)
    integer             :: get_sub
    integer, intent(in) :: x, y
    get_sub = m_sub ( modulo(x-1,m_xsize) + 1, modulo(y-1,m_ysize) + 1 )
    return 
  end function get_sub

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief Assign substrate point at (x,y)
  !> @param m_sub: substrate structure
  !> @param x: x coordinate
  !> @param y: y coordinate
  !> @param v: value
  !--------------------------------------------------------------------------- 
  subroutine set_sub (x, y, v)
    integer, intent(in) :: x, y, v
    m_sub ( modulo(x-1,m_xsize) + 1, modulo(y-1,m_ysize) + 1 ) = v
    return
  end subroutine set_sub

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief to land one molecule at xc, yc, dc position
  !> @param m: molecule index
  !> @param xc
  !> @param yc
  !> @param dc
  !> @return true if landing 
  !--------------------------------------------------------------------------- 
  function land_one (m, xc, yc, dc)
    logical :: land_one
    logical :: empty
    integer, intent(in)  :: m, xc, yc, dc
    integer, allocatable :: vec(:,:)
    integer              :: i, x, y, v, t, s, c

    type(mtype)   , pointer :: t_obj
    type(molecule), pointer :: m_obj

    ! retrieve molecule type & molecule object
    m_obj => mlist(m)
    t_obj => tlist(m_obj % tp)

    ! calculate all positions {x,y,comp-id}
    call alloc (vec, 2, t_obj % comp_num())

    ! check position is empty by the way
    empty = .true.
    do i = 1, t_obj % comp_num()
       vec(:,i) = t_obj % rotate(t_obj % xy(i), dc)
       x = vec(1,i) + xc
       y = vec(2,i) + yc
       if (get_sub(x, y) /= 0) empty = .false.
    end do

    ! land molecule if the site is empty, do nothing otherwise
    if (empty) then

       m_obj % pos(1) = xc ! however this works if the variable is defined
       m_obj % pos(2) = yc ! as a pointer
       m_obj % pos(3) = modulo(dc, t_obj % symm) 

       ! land each dots
       do i = 1, t_obj % comp_num()
          ! part [1]
          t = t_obj % idx_gen() ! type index generated
          c = i                 ! comp-id
          s = m_obj % state(i)  ! comp-state
          ! part [2]
          x = vec(1,i) + xc
          y = vec(2,i) + yc
          v = convert_to_land (m, t, c, s)
          ! part 3
          call set_sub(x, y, v)
       end do
       land_one = .true.
    else
       print *, "  landing failed once"
       land_one = .false.
    end if

    return
  end function land_one

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief to move m-th molecule from its original position to (nx,ny,nd)
  !> @param m: molecule index
  !> @param rx
  !> @param ry
  !> @param rd
  !> @remark: this function does not perform empty test ! user has the duty to
  !           make sure the destination is completely empty
  !---------------------------------------------------------------------------
  subroutine move_one (m, rx, ry, rd)
    integer, intent(in)  :: m, rx, ry, rd ! relative position
    integer              :: nx, ny, nd
    integer              :: ox, oy, od
    integer              :: comp(2), opos(2), npos(2)
    integer              :: i, t, c, s
    type(mtype)   , pointer :: t_obj
    type(molecule), pointer :: m_obj 

    ! retrieve molecule type & molecule object
    m_obj => mlist(m)          ! pointing pointer to oject
    t_obj => tlist(m_obj % tp) ! this is pointer assignment !

    ox = m_obj % pos(1)
    oy = m_obj % pos(2)
    od = m_obj % pos(3)
    
    nx = rx + ox
    ny = ry + oy
    nd = rd + od
    !write (*, '("(",3I2,")=>(",3I2,")")') [ox,oy,od], [nx,ny,nd] 

    !mlist(m) % pos(1) = nx ! here I am directly accessing the data in order
    !mlist(m) % pos(2) = ny ! to make sure the changes are saved. Do changes
    !mlist(m) % pos(3) = nd ! on local variable does not gaurantee that

    m_obj % pos(1) = nx ! however this works if the variable is defined
    m_obj % pos(2) = ny ! as a pointer
    m_obj % pos(3) = modulo(nd, t_obj % symm) 

    ! reset old points
    do i = 1, t_obj % comp_num()
       opos = [ox,oy] + t_obj % rotate(t_obj % xy(i), od)       
       call set_sub(opos(1), opos(2), 0)
    end do

    ! assign new points
    do i = 1, t_obj % comp_num()
       npos = [nx,ny] + t_obj % rotate(t_obj % xy(i), nd)
       t = t_obj % idx_gen() ! type index generated
       c = i                 ! comp-id
       s = m_obj % state(i)  ! comp-state
       call set_sub(npos(1), npos(2), convert_to_land (m, t, c, s))
    end do

    return
  end subroutine move_one

  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Define subroutine for evaporation
  !> @param mtp: copy/pointer to molecule type
  !> @param num: number of molecules will be evaporate
  !--------------------------------------------------------------------------- 
  subroutine evaporate (t_obj, num)
    type(mtype), intent(in) :: t_obj
    integer    , intent(in) :: num
    integer :: i, t, k, x, y, d

    ! debug -------------------------------
    write (*,'(" plan to evaporate type  ",I6,"   for ",I6)') &
         t_obj % idx_def(), num
    ! debug -------------------------------

    ! evaporate new molecules
    t = t_obj % idx_gen() !> get molecule index
    LAND_LOOP: do i = 1, num
       ! check if next landing will exceed maximum molecule number
       if (activated_num(t) >= t_obj % eva_num()) exit LAND_LOOP
       ! land a new molecule
       call activate_new(t)
       SEARCH_LOOP: do while (.true.)
          k = activated_num(t) + t_obj % idx_off() ! offset molecule id
          x = rand_subX()
          y = rand_subY()
          d = rand_int(1, t_obj % symm) 

          ! land_one will return true if landing succeed, vice versa
          if (land_one(k, x, y, d)) exit SEARCH_LOOP
       end do SEARCH_LOOP
    end do LAND_LOOP

    ! debug -------------------------------
    write (*,'(" done evaporation of type",I6,"   for ",I6)') &
         t_obj % idx_def(),activated_num(t)
    ! debug -------------------------------

    return
  end subroutine evaporate

  !--------------------------------------------------------------------------- 
  ! DESCRIPTION
  !> @brief print subsrate to screen for debuging
  !--------------------------------------------------------------------------- 
  subroutine print_to(u, k)
    integer, intent(in) :: u, k
    integer :: i, j, v
    !> @remark values for k: 1 -> mid, 2 ->tid, 3 ->comp 4 -> state
    if (u /= 6) then
       ! printing to files
       do j = 1, m_ysize
          ! in case you want to print all substrate values
          write (u, "(100G2.5)") &
               (convert_from_land(m_sub(i,j), k), i=1,m_xsize)
       end do
    else
       ! printing to screen
       ! do j = 1, m_ysize
       !    ! in case you want to print all substrate values
       !    do i = 1, m_xsize
       !       v = convert_from_land(m_sub(i,j),k)
       !       if (v < 0) then
       !          print *, m_sub(i,j)
       !          stop
       !       end if
       !       write (u, "(I2)",advance="no") v
       !    end do
       !    write (u,*) ""
       ! end do
       do j = 0, m_ysize+1
          ! skip printing when substrate is too large
#ifdef _WIN32
          if (m_xsize > 25) return 
#else
          if (m_xsize > 45) return 
#endif
          ! draw left boundary
          write (u, "(A1)",advance="no") "|"
          do i = 1, m_xsize
             if (j == 0 .or. j == m_ysize + 1) then
                ! draw upper and lower boundary 
                write (u, "(A2)",advance="no") "=="
             else
                ! print values
                v = convert_from_land(m_sub(i,j),k) 
                if (v /= 0) then
                   write (u, "(I2)",advance="no") v
                else 
                   write (u, "(A2)",advance="no") "  "
                end if
             end if
          end do
          ! draw right boundary
          write (u, "(A1)",advance="no") "|"
          ! change a new line
          write (u,*) ""
       end do
    end if
    return
  end subroutine print_to

end module func_substrate
