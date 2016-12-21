module func_xml_reader

  use global
  use func_xml_var
  use FoX_dom
  implicit none

  interface getDataAttribute
     module procedure :: getDataAttribute_integer, getDataAttribute_real
  end interface getDataAttribute
contains

  subroutine getDataAttribute_integer(p, v, name)
    type(Node), pointer     :: p
    integer,    intent(out) :: v
    character(*) :: name
    if (hasAttribute(p, name)) then
       call extractDataAttribute(p, name, v)
    else
       print *, "ERROR", name, " no found"
       stop
    end if
    return
  end subroutine getDataAttribute_integer

  subroutine getDataAttribute_real(p, v, name)
    type(Node), pointer     :: p
    real(dp),   intent(out) :: v
    character(*) :: name
    if (hasAttribute(p, name)) then
       call extractDataAttribute(p, name, v)
    else
       print *, "ERROR", name, " no found"
       stop
    end if
    return
  end subroutine getDataAttribute_real


  subroutine xml_read()
    type(Node)    , pointer :: myDoc, d, r, m, c, p
    type(NodeList), pointer :: dList, rList, mList, cList, pList

    integer :: i,j, k, v
    integer :: vec(3)
    type(var), allocatable :: varList (:)

    print *, "--- reading xml input ---"

    ! Load in the document
    myDoc => parseFile("/home/qiwu/work/dev/simula/fortran/main/src/h2o.xml")

    ! Find all the reactions:
    dList => getElementsByTagName(myDoc, "defReaction")
    mList => getElementsByTagName(myDoc, "molecule"   )

    ! print
    print*, "Found ", getLength(dList), " reaction definitions."
    print*, "Found ", getLength(mList), " molecule definitions."

    ! Loop over the parameter list. Note that the DOM counts from zero, not
    ! from one.
    do i = 0, getLength(mList)-1
       m => item(mList, i)
       
       ! get molecule name
       call getDataAttribute(m, v, "id"    )
       print *, "read id      ", v
       call getDataAttribute(m, v, "symm"  )
       print *, "read symmetry", v
       call getDataAttribute(m, v, "amount")       
       print *, "read amount  ", v

       ! component
       cList => getElementsByTagName(m, "component")
       do j = 0, getLength(cList)-1
          c => item(cList, j)

          call getDataAttribute(c, vec(1), "x")       
          call getDataAttribute(c, vec(2), "y")       
          call getDataAttribute(c, vec(3), "init")       
          print *, "component", vec

       end do
       
       ! reaction
       rList => getElementsByTagName(m, "reaction")
       do j = 0, getLength(rList)-1
          r => item(rList, j)
          pList => getElementsByTagName(r, "var")
          call read_var(pList, varList)

          ! print
          do k = 1, size(varList)
             call varList(k) % print()
          end do

       end do

    end do
    
    print *, "-----------------------"

    ! Clear up all allocated memory
    call destroy(myDoc)

  end subroutine xml_read

  subroutine read_var(nlist, vlist)
    type(Node)    , pointer :: p
    type(NodeList), pointer :: nlist
    type(var), allocatable  :: vlist (:)
    integer :: i, length
    integer               :: sI (1)
    real(dp)              :: sF (1)
    integer , allocatable :: vI (:)
    real(dp), allocatable :: vF (:)
    character(len=20)     :: name

    ! allocate memory
    if (allocated(vlist)) stop "ERROR non-empty var list"
    allocate(vlist(getLength(nlist)))    
    ! loop over all variable
    do i = 0, getLength(nlist)-1
       p => item(nList, i)
       ! for each element
       if (hasAttribute(p, "name")) then
          if (hasAttribute(p, "type")) then
             if (getAttribute(p, "type")=="scalarInt") then

                call extractDataContent(p, sI)
                call vlist(i+1) % encode(sI, getAttribute(p, "name"))

             elseif (getAttribute(p, "type")=="scalarFloat") then

                call extractDataContent(p, sF)              
                call vlist(i+1) % encode(sF, getAttribute(p, "name"))

             elseif (getAttribute(p, "type")=="vectorInt") then

                if (hasAttribute(p, "length")) then                      
                   call extractDataAttribute(p, "length", length)
                   call alloc_I1(vI, length)
                   call extractDataContent(p, vI)
                   call vlist(i+1) % encode(vI, getAttribute(p, "name"))
                   deallocate(vI)
                else
                   stop "ERROR (XML) vector length missing"
                end if

             elseif (getAttribute(p, "type")=="vectorFloat") then

                if (hasAttribute(p, "length")) then                      
                   call extractDataAttribute(p, "length", length)
                   call alloc_F1(vF, length)
                   call extractDataContent(p, vF)
                   call vlist(i+1) % encode(vF, getAttribute(p, "name"))
                   deallocate(vF)
                else
                   stop "ERROR (XML) vector length missing"
                end if

             else
                stop "ERROR (XML) unknown data type"
             end if

          else
             stop "ERROR (XML) missing variable type"
          end if

       else
          stop "ERROR (XML) missing variable name"
       end if

    end do
    return
  end subroutine read_var

end module func_xml_reader
