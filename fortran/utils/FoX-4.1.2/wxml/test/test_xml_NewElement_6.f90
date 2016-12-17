program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement
  use FoX_common, only : FoX_set_fatal_errors
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call FoX_set_fatal_errors(.true.)
  call xml_OpenFile(filename, xf)
  call xml_NewElement(xf, "!root")
  call xml_Close(xf)

end program test
