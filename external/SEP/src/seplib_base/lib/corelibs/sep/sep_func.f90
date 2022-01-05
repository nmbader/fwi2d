module sep_func_mod
use, intrinsic :: iso_c_binding
implicit none



 
interface 

    subroutine sepwarn(err,string) bind(c,name="sepwarn")
    import
    integer(c_int), intent(in) :: err
    character(C_CHAR), dimension(*),intent(in) :: string
  end subroutine
    subroutine fseperr(string) bind(c,name="seperr")
    import
    character(C_CHAR), dimension(*),intent(in) :: string
  end subroutine

  subroutine sep_add_doc_line(string) bind(c,name="sep_add_doc_line")
    import
    character(C_CHAR), dimension(*) :: string
  end subroutine
  subroutine sep_prog(string) bind(c,name="sep_prog")
    import
    character(C_CHAR), dimension(*) :: string
  end subroutine
  subroutine sep_begin_prog() bind(c,name="sep_begin_prog")
    import
 end subroutine
   subroutine sep_end_prog() bind(c,name="sep_end_prog")
    import
end subroutine
integer function strlen(strin) bind(c,name="strlen")
  import
  character(C_CHAR),dimension(*) :: strin
end function
end interface

contains
  subroutine erexit(string) 
    character(len=*) :: string
    call seperr(string//C_NULL_CHAR)
  end subroutine

  subroutine seperr(string) 
    character(len=*) :: string
    call fseperr(string//C_NULL_CHAR)
    call exit(-1);
  end subroutine

subroutine c2forstr(str)
 character(len=*) :: str
 integer ::n,i
 
 n=strlen(str)
   do i=n+1,len(str)
     str(i:i)=" "
   end do
end subroutine


end module
