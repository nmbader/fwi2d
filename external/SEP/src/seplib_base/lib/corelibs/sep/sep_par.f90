module sep_par_mod
use, intrinsic :: iso_c_binding
use sep_func_mod
implicit none



 
interface 

 
 subroutine sep_prog_stat(nm,i1,i2,i3) bind(c,name="sep_prog_stat")
  import
    character(C_CHAR), dimension(*) :: nm
    integer(C_INT),value :: i2,i1,i3
  end subroutine
 
  subroutine doc1(string) bind(c,name="doc")
    import
    character(C_CHAR), dimension(*) :: string
  end subroutine 

  subroutine putlin(string) bind(c,name="putlin")
    import
    character(C_CHAR), dimension(*) :: string
  end subroutine

  subroutine auxputlin(tag,string) bind(c,name="sfauxputlin")
    import
    character(C_CHAR), dimension(*) :: string,tag 
  end subroutine
 
 
 !GETCH

   integer function getch_l_f(arg,typ,val) bind(c,name="getch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
    integer function getch_s_f(arg,typ,val) bind(c,name="getch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    character(C_CHAR),dimension(*) :: val
  end function
    integer function getch_f_f(arg,typ,val) bind(c,name="getch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT) :: val
  end function
    integer function getch_g_f(arg,typ,val) bind(c,name="getch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE) :: val
  end function
   integer function getch_i_f(arg,typ,val) bind(c,name="getch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
      integer function getch_f_f_a(arg,typ,val) bind(c,name="getch_f_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT),dimension(*) :: val
  end function
      integer function getch_g_f_a(arg,typ,val) bind(c,name="getch_g_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE),dimension(*) :: val
  end function
   integer function getch_i_f_a(arg,typ,val) bind(c,name="getch_i_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT),dimension(*) :: val
  end function
  
  
  !PUTCH
    integer function putch_i_f(arg,typ,val) bind(c,name="putch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
  
    integer function putch_s_f(arg,typ,val) bind(c,name="putch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    character(C_CHAR),dimension(*) :: val
  end function
    integer function putch_f_f(arg,typ,val) bind(c,name="putch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT) :: val
  end function
    integer function putch_l_f(arg,typ,val) bind(c,name="putch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    logical(C_BOOL) :: val
  end function
    integer function putch_g_f(arg,typ,val) bind(c,name="putch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE) :: val
  end function
  
    !FETCH
    integer function fetch_i_f(arg,typ,val) bind(c,name="fetch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
  
    integer function fetch_s_f(arg,typ,val) bind(c,name="fetch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    character(C_CHAR),dimension(*) :: val
  end function
    integer function fetch_f_f(arg,typ,val) bind(c,name="fetch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT) :: val
  end function
    integer function fetch_l_f(arg,typ,val) bind(c,name="fetch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    logical(C_BOOL) :: val
  end function
    integer function fetch_g_f(arg,typ,val) bind(c,name="fetch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE) :: val
  end function
  
      !TETCH
    integer function tetch_i_f(arg,typ,val) bind(c,name="tetch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
  
    integer function tetch_s_f(arg,typ,val) bind(c,name="tetch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    character(C_CHAR),dimension(*) :: val
  end function
    integer function tetch_f_f(arg,typ,val) bind(c,name="tetch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT) :: val
  end function
    integer function tetch_l_f(arg,typ,val) bind(c,name="tetch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    logical(C_BOOL) :: val
  end function
    integer function tetch_g_f(arg,typ,val) bind(c,name="tetch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE) :: val
  end function
  
  !HETCH
    integer function hetch_i_f(arg,typ,val) bind(c,name="hetch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT) :: val
  end function
  
    integer function hetch_s_f(arg,typ,val) bind(c,name="hetch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    character(C_CHAR),dimension(*) :: val
  end function
    integer function hetch_f_f(arg,typ,val) bind(c,name="hetch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT) :: val
  end function
    integer function hetch_l_f(arg,typ,val) bind(c,name="hetch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    logical(C_BOOL) :: val
  end function
    integer function hetch_g_f(arg,typ,val) bind(c,name="hetch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE) :: val
  end function
      integer function hetch_f_f_a(arg,typ,val) bind(c,name="hetch_f_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_FLOAT),dimension(*) :: val
  end function
      integer function hetch_i_f_a(arg,typ,val) bind(c,name="hetch_i_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    integer(C_INT),dimension(*) :: val
  end function
      integer function hetch_l_f_a(arg,typ,val) bind(c,name="hetch_l_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    logical(C_BOOL),dimension(*) :: val
  end function
      integer function hetch_g_f_a(arg,typ,val) bind(c,name="hetch_g_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ
    real(C_DOUBLE),dimension(*) :: val
  end function
  
  !AUXP
      integer function auxpar_l_f(arg,typ,val,tag) bind(c,name="auxpar_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    logical(C_BOOL) :: val
  end function
    integer function auxpar_s_f(arg,typ,val,tag) bind(c,name="auxpar_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    character(C_CHAR),dimension(*) :: val
  end function
    integer function auxpar_f_f(arg,typ,val,tag) bind(c,name="auxpar_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_FLOAT) :: val
  end function
    integer function auxpar_i_f(arg,typ,val,tag) bind(c,name="auxpar_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    integer(C_INT) :: val
  end function
    integer function auxpar_g_f(arg,typ,val,tag) bind(c,name="auxpar_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_DOUBLE) :: val
  end function
    integer function auxpar_f_f_a(arg,typ,val,tag) bind(c,name="auxpar_f_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_FLOAT),dimension(*) :: val
  end function
    integer function auxpar_g_f_a(arg,typ,val,tag) bind(c,name="auxpar_g_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_DOUBLE),dimension(*) :: val
  end function
    integer function auxpar_i_f_a(arg,typ,val,tag) bind(c,name="auxpar_i_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    integer(C_INT),dimension(*) :: val
  end function
  
  
  
    !AUXPUTCH
  integer function auxputch_i_f_a(arg,typ,val,tag) bind(c,name="auxputch_i_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    integer(C_INT),dimension(*) :: val
  end function
  integer function auxputch_f_f_a(arg,typ,val,tag) bind(c,name="auxputch_f_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_FLOAT),dimension(*) :: val
  end function
  integer function auxputch_g_f_a(arg,typ,val,tag) bind(c,name="auxputch_g_f_a")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_DOUBLE),dimension(*) :: val
  end function


    integer function auxputch_i_f(arg,typ,val,tag) bind(c,name="auxputch_i_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    integer(C_INT) :: val
  end function
  
    integer function auxputch_s_f(arg,typ,val,tag) bind(c,name="auxputch_s_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    character(C_CHAR),dimension(*) :: val
  end function
    integer function auxputch_f_f(arg,typ,val,tag) bind(c,name="auxputch_f_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_FLOAT) :: val
  end function
    integer function auxputch_l_f(arg,typ,val,tag) bind(c,name="auxputch_l_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    logical(C_BOOL) :: val
  end function
    integer function auxputch_g_f(arg,typ,val,tag) bind(c,name="auxputch_g_f")
    import
    character(C_CHAR),dimension(*) ::arg,typ,tag
    real(C_DOUBLE) :: val
  end function
  
  subroutine init_args(prog_name) bind(c,name="initpar_f")
    import
    character(C_CHAR) , dimension(*) :: prog_name
  end subroutine

  subroutine getch_add_string2(string) bind(c,name="getch_add_string")
    import
    character(C_CHAR), dimension(*) :: string
  end subroutine

  subroutine hclose() bind(c,name="hclose")
    import
  end subroutine
  integer function sep_get_number_data_axes(tag,nd) bind(c,name="sep_get_number_data_axes")
  import
    character(C_CHAR),dimension(*),intent(in) :: tag
    integer :: nd
  end function
  
     integer function sep_get_data_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_get_data_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(out) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(out) :: n
      real(C_FLOAT),intent(out) :: o,d
  end  function
  
  
   integer function sep_put_data_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_put_data_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(in) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(in) :: n
      real(C_FLOAT),intent(in) :: o,d
  end  function
 
end interface

interface tetch
  module procedure tet_i_f,tet_s_f,tet_f_f,tet_l_f,tet_g_f
end interface
interface fetch
  module procedure fet_i_f,fet_s_f,fet_f_f,fet_l_f,fet_g_f
end interface
interface getch
  module procedure get_i_f,get_s_f,get_f_f,get_l_f,get_g_f,get_f_f_a,get_i_f_a,get_g_f_a
end interface
interface putch
  module procedure put_i_f,put_s_f,put_f_f,put_l_f,put_g_f
end interface
interface hetch
  module procedure het_i_f,het_s_f,het_f_f,het_l_f,het_g_f,het_i_f_a,het_f_f_a,het_l_f_a,het_g_f_a
end interface
interface auxpar
  module procedure auxp_i_f,auxp_s_f,auxp_f_f,auxp_l_f,auxp_g_f,auxp_i_f_a,auxp_f_f_a,auxp_g_f_a
end interface

interface auxputch
  module procedure auxpu_i_f,auxpu_s_f,auxpu_f_f,auxpu_l_f,auxpu_g_f,&
    auxpu_i_f_a,auxpu_f_f_a,auxpu_g_f_a
end interface

interface doc
 module procedure doc_0,doc_1
end interface
contains


subroutine doc_0()
  call doc1("")
end subroutine

subroutine doc_1(str)
  character(len=*) :: str
  call doc1(str)
end subroutine

subroutine getch_add_string(str)
  character(len=*) :: str
  call getch_add_string2(str//C_NULL_CHAR)
end subroutine

integer function get_s_f(arg,typ,val)
  character(len=*) arg,typ,val
  get_s_f= getch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
  call c2forstr(val)
end function
integer function get_i_f(arg,typ,val)
  character(len=*) arg,typ
  integer :: val
  get_i_f= getch_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_f_f(arg,typ,val)
  character(len=*) arg,typ
  real :: val
  get_f_f= getch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_g_f(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)) :: val
  get_g_f= getch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_i_f_a(arg,typ,val)
  character(len=*) arg,typ
  integer,dimension(:) :: val
  get_i_f_a= getch_i_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_f_f_a(arg,typ,val)
  character(len=*) arg,typ
  real,dimension(:) :: val
  get_f_f_a= getch_f_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_g_f_a(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)),dimension(:) :: val
  get_g_f_a= getch_g_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function get_l_f(arg,typ,val)
  character(len=*) arg,typ
  logical :: val
  integer :: v
  get_l_f= getch_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,v)
  val=.false.
  if(v==1)  val=.true.;
end function
!TETCH

integer function tet_i_f(arg,typ,val)
  character(len=*) arg,typ
  integer :: val
  tet_i_f= tetch_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

integer function tet_l_f(arg,typ,val)
  character(len=*) arg,typ
  logical(C_BOOL) :: val
  tet_l_f= tetch_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

integer function tet_s_f(arg,typ,val)
  character(len=*) arg,typ,val
  tet_s_f= tetch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
  call c2forstr(val)
end function

integer function tet_f_f(arg,typ,val)
  character(len=*) arg,typ
  real :: val
  tet_f_f= tetch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

integer function tet_g_f(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)) :: val
  tet_g_f= tetch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

!FETCH

integer function fet_i_f(arg,typ,val)
  character(len=*) arg,typ
  integer :: val
  fet_i_f= fetch_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

integer function fet_s_f(arg,typ,val)
  character(len=*) arg,typ,val
  fet_s_f= fetch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
  call c2forstr(val)
end function

integer function fet_f_f(arg,typ,val)
  character(len=*) arg,typ
  real :: val
  fet_f_f= fetch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function fet_g_f(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)) :: val
  fet_g_f= fetch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function fet_l_f(arg,typ,val)
  character(len=*) arg,typ
  logical :: val
  logical*1::v2
  v2=val
  fet_l_f= fetch_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,v2)
end function

integer function put_i_f(arg,typ,val)
  character(len=*) arg,typ
  integer :: val
  put_i_f= putch_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function

integer function put_s_f(arg,typ,val)
  character(len=*) arg,typ,val
  put_s_f= putch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,trim(val)//C_NULL_CHAR)
end function

integer function put_f_f(arg,typ,val)
  character(len=*) arg,typ
  real :: val
  put_f_f= putch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function put_g_f(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)) :: val
  put_g_f= putch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function put_l_f(arg,typ,val)
  character(len=*) arg,typ
  logical :: val
  logical*1 :: v2
  v2=val
  put_l_f= putch_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,v2)
end function



!HETCH
integer function het_s_f(arg,typ,val)
  character(len=*) arg,typ,val
  het_s_f= hetch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
  call c2forstr(val)
end function

integer function het_f_f(arg,typ,val)
  character(len=*) arg,typ
  real :: val
  het_f_f= hetch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function het_g_f(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)) :: val
  het_g_f= hetch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function het_i_f(arg,typ,val)
  character(len=*) arg,typ
  integer :: val
  het_i_f= hetch_i_f(trim(arg)//C_NULL_CHAR,typ,val)
end function
integer function het_l_f(arg,typ,val)
  character(len=*) arg,typ
  logical(C_BOOL) :: val
  het_l_f= hetch_l_f(trim(arg)//C_NULL_CHAR,typ,val)
end function
integer function het_f_f_a(arg,typ,val)
  character(len=*) arg,typ
  real,dimension(:) :: val
  het_f_f_a= hetch_f_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function het_g_f_a(arg,typ,val)
  character(len=*) arg,typ
  real (KIND=KIND(1.0D0)),dimension(:) :: val
  het_g_f_a= hetch_g_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function het_i_f_a(arg,typ,val)
  character(len=*) arg,typ
  integer,dimension(:) :: val
  het_i_f_a= hetch_i_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function
integer function het_l_f_a(arg,typ,val)
  character(len=*) arg,typ
  logical(C_BOOL),dimension(:) :: val
  het_l_f_a= hetch_l_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val)
end function


!AUXPAR

integer function auxp_l_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  logical :: val
  logical*1 :: v2
  v2=val
  auxp_l_f= auxpar_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,v2,trim(tag)//C_NULL_CHAR)
end function

integer function auxp_s_f(arg,typ,val,tag)
  character(len=*) arg,typ,val,tag
  auxp_s_f= auxpar_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
  call c2forstr(val)
end function
integer function auxp_i_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  integer :: val
  auxp_i_f= auxpar_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function
integer function auxp_f_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real :: val
  auxp_f_f= auxpar_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function
integer function auxp_g_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real (KIND=KIND(1.0D0)) :: val
  auxp_g_f= auxpar_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

  
integer function auxp_i_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  integer,dimension(:) :: val
  auxp_i_f_a= auxpar_i_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function
integer function auxp_f_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real,dimension(:) :: val
  auxp_f_f_a= auxpar_f_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function
integer function auxp_g_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real (KIND=KIND(1.0D0)),dimension(:) :: val
  auxp_g_f_a= auxpar_g_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

!AUXPUTCH
integer function auxpu_l_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  logical :: val
  logical*1:: v2
  v2=val
  auxpu_l_f= auxputch_l_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,v2,trim(tag)//C_NULL_CHAR)
end function
integer function auxpu_i_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  integer :: val
  auxpu_i_f= auxputch_i_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

integer function auxpu_s_f(arg,typ,val,tag)
  character(len=*) arg,typ,val,tag
  auxpu_s_f= auxputch_s_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,trim(val)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
end function

integer function auxpu_f_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real :: val
  auxpu_f_f= auxputch_f_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

integer function auxpu_g_f(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real (KIND=KIND(1.0D0)) :: val
  auxpu_g_f= auxputch_g_f(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

integer function auxpu_i_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  integer :: val(:)
  auxpu_i_f_a= auxputch_i_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function


integer function auxpu_f_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real :: val(:)
  auxpu_f_f_a= auxputch_f_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function

integer function auxpu_g_f_a(arg,typ,val,tag)
  character(len=*) arg,typ,tag
  real (KIND=KIND(1.0D0)) :: val(:)
  auxpu_g_f_a= auxputch_g_f_a(trim(arg)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,val,trim(tag)//C_NULL_CHAR)
end function


subroutine initpar()
  integer i
  character(len=9999) :: nm
  call get_command_argument(0,nm)
  call init_args(trim(nm)//C_NULL_CHAR)
  nm=trim(nm)//C_NULL_CHAR
  DO I=1,COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(I,nm)
    call getch_add_string(trim(nm)//C_NULL_CHAR)
  END DO 
end subroutine

  integer function sep_get_data_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(out) :: n
    real,intent(out) :: o,d
    character(len=*) :: label
     sep_get_data_axis_par= sep_get_data_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,label)
     call c2forstr(label)
  end function
  
     integer function sep_put_data_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(in) :: n
    real,intent(in) :: o,d
    character(len=*) ,intent(in):: label
     sep_put_data_axis_par= sep_put_data_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,trim(label)//C_NULL_CHAR)
  end function
  



end module
