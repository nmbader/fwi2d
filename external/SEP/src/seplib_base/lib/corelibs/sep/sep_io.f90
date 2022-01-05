module sep_io_mod
use, intrinsic :: iso_c_binding
use sep_func_mod
implicit none



 
interface 
      integer function sreed_i(tag,ar,sz) bind(c,name="sreed_i")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    integer(C_INT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
      integer function sreed_f(tag,ar,sz) bind(c,name="sreed_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    real(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
        integer function sreed_c(tag,ar,sz) bind(c,name="sreed_c")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    complex(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in) ,value:: sz
  end function
        integer function srite_f(tag,ar,sz) bind(c,name="srite_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    real(C_FLOAT),dimension(*),intent(in) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
       integer function srite_c(tag,ar,sz) bind(c,name="srite_c")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    complex(C_FLOAT),dimension(*),intent(in) :: ar
    integer(C_INT),intent(in),value :: sz
  end function

      integer(kind=8) function sreedll_i(tag,ar,sz) bind(c,name="sreedll_i")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    integer(C_INT),dimension(*),intent(out) :: ar
    integer(C_LONG_LONG),intent(in),value :: sz
  end function
      integer(kind=8) function sreedll_f(tag,ar,sz) bind(c,name="sreedll_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    real(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_LONG_LONG),intent(in),value :: sz
  end function
        integer(kind=8) function sreedll_c(tag,ar,sz) bind(c,name="sreedll_x")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    complex(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_LONG_LONG),intent(in) ,value:: sz
  end function
        integer(kind=8) function sritell_f(tag,ar,sz) bind(c,name="sritell_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    real(C_FLOAT),dimension(*),intent(in) :: ar
    integer(C_LONG_LONG),intent(in),value :: sz
  end function
       integer(kind=8) function sritell_c(tag,ar,sz) bind(c,name="sritell_x")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    complex(C_FLOAT),dimension(*),intent(in) :: ar
    integer(C_LONG_LONG),intent(in),value :: sz
  end function
   integer(kind=8) function sritell_i(tag,ar,sz) bind(c,name="sritell_c")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    integer(C_INT),dimension(*),intent(in) :: ar
    integer(C_LONG_LONG),intent(in) ,value:: sz
  end function

!   integer function input() bind(c,name="sfdinput")
!    import
!  end function
!   integer function output() bind(c,name="sfdoutput")
!    import
!  end function

   integer function srite_i(tag,ar,sz) bind(c,name="srite")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag
    integer(C_INT),dimension(*),intent(in) :: ar
    integer(C_INT),intent(in) ,value:: sz
  end function

  integer function sreed_window_f(tag,nd,ng,n,f,j,sz,buf) bind(c,name="sreed_window")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd


    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    real(C_FLOAT), dimension(*),intent(out) :: buf
  end function
  integer function sreed_window_i(tag,nd,ng,n,f,j,sz,buf) bind(c,name="sreed_window_i")
  import
  character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd

    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    integer(C_INT), dimension(*),intent(out) :: buf
  end function
  integer function sreed_window_c(tag,nd,ng,n,f,j,sz,buf) bind(c,name="sreed_window_c")
  import
  character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd

    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    complex(C_FLOAT), dimension(*),intent(out) :: buf
  end function

    integer function srite_window_f(tag,nd,ng,n,f,j,sz,buf) bind(c,name="srite_window_f")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
     integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd

    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    real(C_FLOAT), dimension(*),intent(in) :: buf
  end function
    integer function srite_window_i(tag,nd,ng,n,f,j,sz,buf)bind(c,name="srite_window_i")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd

    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    integer(C_INT), dimension(*),intent(in) :: buf
  end function
    integer function srite_window_c(tag,nd,ng,n,f,j,sz,buf)bind(c,name="srite_window_c")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: sz
        integer(C_INT), intent(in):: nd

    integer(C_INT),dimension(*),intent(in) :: n,f,j,ng
    complex(C_FLOAT), dimension(*),intent(in) :: buf
  end function
  
  subroutine auxtmpf(tag) bind(c,name="auxtmp")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end subroutine
  integer function auxinf(tag) bind(c,name="fauxin")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end function
  subroutine auxclosef(tag) bind(c,name="auxclose")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end subroutine
  subroutine auxoutf(tag) bind(c,name="auxout")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end subroutine
  subroutine auxinoutf(tag) bind(c,name="auxinout")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end subroutine

    subroutine aux_unlinkf(tag) bind(c,name="aux_unlink")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
  end subroutine
  
    integer function sseek(tag,pos,beg) bind(c,name="sseek")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: pos,beg
  end function
  
      integer function sseek_block(tag,pos,bs,beg) bind(c,name="sseek_block")
    import
    character(C_CHAR),dimension(*), intent(in) :: tag
    integer(C_INT), intent(in),value :: pos,beg,bs
  end function
  
  
        integer function sreed2_f(tag,ar,sz,typ) bind(c,name="sreed2_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag,typ
    real(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
       integer function sreed2_i(tag,ar,sz,typ) bind(c,name="sreed2_i")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag,typ
    integer(C_INT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
 
 
       integer function srite2_f(tag,ar,sz,typ) bind(c,name="srite2_f")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag,typ
    real(C_FLOAT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
       integer function srite2_i(tag,ar,sz,typ) bind(c,name="srite2_i")
    import
    character(C_CHAR),dimension(*),intent(in) ::tag,typ
    integer(C_INT),dimension(*),intent(out) :: ar
    integer(C_INT),intent(in),value :: sz
  end function
  
end interface


interface sreed
  module procedure sre_f_1,sre_f_2,sre_f_3,sre_f_4,sre_f_5,sre_f_0,sre_i_2
  module procedure sre_c_1,sre_c_2,sre_c_3,sre_c_4,sre_c_5,sre_i_0,sre_i_1
    module procedure sre_f_1_ll,sre_f_2_ll,sre_f_3_ll,sre_f_4_ll,sre_f_5_ll,sre_f_0_ll,sre_i_2_ll
  module procedure sre_c_1_ll,sre_c_2_ll,sre_c_3_ll,sre_c_4_ll,sre_c_5_ll,sre_i_0_ll,sre_i_1_ll
end interface

interface sreed2
  module procedure sre2_f_1,sre2_i_1,sre2_i_0
end interface

interface srite2
  module procedure sri2_f_1,sri2_i_1,sri2_i_0
end interface

interface srite
  module procedure sri_f_1,sri_f_2,sri_f_3,sri_f_4,sri_f_5
  module procedure sri_i_2,sri_i_0,sri_i_1
  module procedure sri_c_1,sri_c_2,sri_c_3,sri_c_4,sri_c_5
  
    module procedure sri_f_1_ll,sri_f_2_ll,sri_f_3_ll,sri_f_4_ll,sri_f_5_ll
  module procedure sri_i_2_ll,sri_i_0_ll,sri_i_1_ll
  module procedure sri_c_1_ll,sri_c_2_ll,sri_c_3_ll,sri_c_4_ll,sri_c_5_ll
end interface

interface sreed_window
  module procedure srw_f_1,srw_f_2,srw_f_3,srw_f_4,srw_f_5
  module procedure srw_i_1,srw_i_2,srw_i_3,srw_i_4,srw_i_5
  module procedure srw_c_1,srw_c_2,srw_c_3,srw_c_4,srw_c_5
end interface

interface srite_window
  module procedure sww_f_1,sww_f_2,sww_f_3,sww_f_4,sww_f_5
  module procedure sww_i_1,sww_i_2,sww_i_3,sww_i_4,sww_i_5
  module procedure sww_c_1,sww_c_2,sww_c_3,sww_c_4,sww_c_5
end interface

contains


subroutine auxtmp(tag)
  character(len=*) :: tag
  call auxtmpf(trim(tag)//C_NULL_CHAR)
end subroutine

integer function auxin(tag)
  character(len=*) :: tag
  auxin=auxinf(trim(tag)//C_NULL_CHAR)
end function

subroutine auxclose(tag)
  character(len=*) :: tag
  call auxclosef(trim(tag)//C_NULL_CHAR)
end subroutine
subroutine auxout(tag)
  character(len=*) :: tag
  call auxoutf(trim(tag)//C_NULL_CHAR)
end subroutine
subroutine auxinout(tag)
  character(len=*) :: tag
  call auxinoutf(trim(tag)//C_NULL_CHAR)
end subroutine

subroutine aux_unlink(tag)
  character(len=*) :: tag
  call aux_unlinkf(trim(tag)//C_NULL_CHAR)
end subroutine


!SREED2
integer function sri2_f_1(tag,ar,sz,typ)
  character(len=*) tag,typ
  real,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sri2_f_1= srite2_f(trim(tag)//C_NULL_CHAR,ar,sz,typ)
end function
integer function sri2_i_1(tag,ar,sz,typ)
  character(len=*) tag,typ
  integer,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sri2_i_1= srite2_i(trim(tag)//C_NULL_CHAR,ar,sz,typ)
end function
integer function sri2_i_0(tag,ar,sz,typ)
  character(len=*) tag,typ
  integer   :: ar,x(1)
    integer, intent(in) :: sz
  x(1)=ar
  sri2_i_0= srite2_i(trim(tag)//C_NULL_CHAR,x,sz,typ)
  
end function

!SREED2
integer function sre2_f_1(tag,ar,sz,typ)
  character(len=*) tag,typ
  real,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre2_f_1= sreed2_f(trim(tag)//C_NULL_CHAR,ar,sz,typ)
end function
integer function sre2_i_1(tag,ar,sz,typ)
  character(len=*) tag,typ
  integer,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre2_i_1= sreed2_i(trim(tag)//C_NULL_CHAR,ar,sz,typ)
end function
integer function sre2_i_0(tag,ar,sz,typ)
  character(len=*) tag,typ
  integer   :: ar,x(1)
    integer, intent(in) :: sz

  sre2_i_0= sreed2_i(trim(tag)//C_NULL_CHAR,x,sz,typ)
  ar=x(1)
end function

!INTEGER READS
integer function sre_i_2(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_i_2= sreed_i(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_i_1(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_i_1= sreed_i(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_i_0(tag,ar,sz)
  character(len=*) tag
  integer,intent(out)  :: ar
    integer, intent(in) :: sz
   integer :: x(1)
  sre_i_0= sreed_i(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function

!FLOAT SREEDS
integer function sre_f_1(tag,ar,sz)
  character(len=*) tag
  real,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_f_1= sreed_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_f_0(tag,ar,sz)
  character(len=*) tag
  real,intent(out)  :: ar
    integer, intent(in) :: sz
    real :: x(1)

  sre_f_0= sreed_f(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function
integer function sre_f_2(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_f_2= sreed_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_f_3(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_f_3= sreed_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_f_4(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_f_4= sreed_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_f_5(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_f_5= sreed_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!COMPLEX SREEDS
integer function sre_c_1(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_c_1= sreed_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_c_2(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_c_2= sreed_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_c_3(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_c_3= sreed_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_c_4(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_c_4= sreed_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sre_c_5(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:,:,:),intent(out)  :: ar
    integer, intent(in) :: sz

  sre_c_5= sreed_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!SRITE FLOATS
integer function sri_f_1(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  real,dimension(:),intent(in)  :: ar
  integer, intent(in) :: sz
  sri_f_1= srite_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_f_0(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  real,intent(out)  :: ar
  real :: x(1)
  integer, intent(in) :: sz
  sri_f_0= srite_f(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function
integer function sri_f_2(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_f_2= srite_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_f_3(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_f_3= srite_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_f_4(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_f_4= srite_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_f_5(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_f_5= srite_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!SREED COMPLEX
integer function sri_c_1(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  complex,dimension(:),intent(in)  :: ar
  integer, intent(in) :: sz
  sri_c_1= srite_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_c_2(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_c_2= srite_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_c_3(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_c_3= srite_c(tag,ar,sz)
end function
integer function sri_c_4(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_c_4= srite_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer function sri_c_5(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:,:,:),intent(in)  :: ar
    integer, intent(in) :: sz

  sri_c_5= srite_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function

!SRITE INTEGERS
integer function sri_i_2(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:,:) ,intent(in) :: ar
    integer, intent(in) :: sz

  sri_i_2= srite_i(trim(tag)//C_NULL_CHAR,ar,sz)

end function
integer function sri_i_0(tag,ar,sz)
  character(len=*) tag
  integer,intent(in) :: ar
    integer, intent(in) :: sz
   integer :: x(1)

  x(1)=ar
  sri_i_0= srite_i(trim(tag)//C_NULL_CHAR,x,sz)

end function
integer function sri_i_1(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:) ,intent(in) :: ar
    integer, intent(in) :: sz

  sri_i_1= srite_i(trim(tag)//C_NULL_CHAR,ar,sz)

end function



!!!!!!!!!!!!!!!AAAAAAAAA
integer(kind=8) function sre_i_2_ll(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_i_2_ll= sreedll_i(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_i_1_ll(tag,ar,sz)
  character(len=*) tag
  integer,dimension(:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_i_1_ll= sreedll_i(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_i_0_ll(tag,ar,sz)
  character(len=*) tag
  integer,intent(out)  :: ar
    integer(kind=8), intent(in) :: sz
   integer :: x(1)
  sre_i_0_ll= sreedll_i(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function

!FLOAT SREEDS
integer(kind=8)  function sre_f_1_ll(tag,ar,sz)
  character(len=*) tag
  real,dimension(:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_f_1_ll= sreedll_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_f_0_ll(tag,ar,sz)
  character(len=*) tag
  real,intent(out)  :: ar
    integer(kind=8), intent(in) :: sz
    real :: x(1)

  sre_f_0_ll= sreedll_f(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function
integer(kind=8)  function sre_f_2_ll(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_f_2_ll= sreedll_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_f_3_ll(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_f_3_ll= sreedll_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_f_4_ll(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_f_4_ll= sreedll_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_f_5_ll(tag,ar,sz)
  character(len=*) tag
  real,dimension(:,:,:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_f_5_ll= sreedll_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!COMPLEX SREEDS
integer(kind=8)  function sre_c_1_ll(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_c_1_ll= sreedll_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_c_2_ll(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_c_2_ll= sreedll_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_c_3_ll(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_c_3_ll= sreedll_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_c_4_ll(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_c_4_ll= sreedll_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sre_c_5_ll(tag,ar,sz)
  character(len=*) tag
  complex,dimension(:,:,:,:,:),intent(out)  :: ar
    integer(kind=8), intent(in) :: sz

  sre_c_5_ll= sreedll_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!SRITE FLOATS
integer(kind=8)  function sri_f_1_ll(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  real,dimension(:),intent(in)  :: ar
  integer(kind=8), intent(in) :: sz
  sri_f_1_ll= sritell_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_f_0_ll(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  real,intent(out)  :: ar
  real :: x(1)
  integer(kind=8), intent(in) :: sz
  sri_f_0_ll= sritell_f(trim(tag)//C_NULL_CHAR,x,sz)
  ar=x(1)
end function
integer(kind=8)  function sri_f_2_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_f_2_ll= sritell_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_f_3_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_f_3_ll= sritell_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_f_4_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_f_4_ll= sritell_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_f_5_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  real,dimension(:,:,:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_f_5_ll= sritell_f(trim(tag)//C_NULL_CHAR,ar,sz)
end function
!SREED COMPLEX
integer(kind=8)  function sri_c_1_ll(tag,ar,sz)
  character(len=*),intent(in)  ::tag
  complex,dimension(:),intent(in)  :: ar
  integer(kind=8), intent(in) :: sz
  sri_c_1_ll= sritell_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_c_2_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_c_2_ll= sritell_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_c_3_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_c_3_ll= sritell_c(tag,ar,sz)
end function
integer(kind=8)  function sri_c_4_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_c_4_ll= sritell_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function
integer(kind=8)  function sri_c_5_ll(tag,ar,sz)
  character(len=*),intent(in) :: tag
  complex,dimension(:,:,:,:,:),intent(in)  :: ar
    integer(kind=8), intent(in) :: sz

  sri_c_5_ll= sritell_c(trim(tag)//C_NULL_CHAR,ar,sz)
end function

!SRITE INTEGERS
integer(kind=8)  function sri_i_2_ll(tag,ar,sz)
  character(len=*),intent(in) ::tag
  integer,dimension(:,:) ,intent(in) :: ar
    integer(kind=8), intent(in) :: sz

  sri_i_2_ll= sritell_i(trim(tag)//C_NULL_CHAR,ar,sz)

end function
integer(kind=8)  function sri_i_0_ll(tag,ar,sz)
  character(len=*),intent(in) ::tag
  integer,intent(in) :: ar
    integer(kind=8), intent(in) :: sz
   integer :: x(1)

  x(1)=ar
  sri_i_0_ll= sritell_i(trim(tag)//C_NULL_CHAR,x,sz)

end function
integer(kind=8)  function sri_i_1_ll(tag,ar,sz)
  character(len=*),intent(in) ::tag
  integer,dimension(:) ,intent(in) :: ar
    integer(kind=8), intent(in) :: sz

  sri_i_1_ll= sritell_i(trim(tag)//C_NULL_CHAR,ar,sz)

end function


!!!!!!!!!!!!!!!!!!!!!AAAAAAA





!sreed window
integer function srw_i_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(out), dimension(:) :: buf
  srw_i_1=sreed_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_i_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(out), dimension(:,:) :: buf
  srw_i_2=sreed_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_i_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(out), dimension(:,:,:) :: buf
  srw_i_3=sreed_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function

integer function srw_i_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(out), dimension(:,:,:,:) :: buf
  srw_i_4=sreed_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function

integer function srw_i_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(out), dimension(:,:,:,:,:) :: buf
  srw_i_5=sreed_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function

integer function srw_f_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(out), dimension(:) :: buf
  srw_f_1=sreed_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_f_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(out), dimension(:,:) :: buf
  srw_f_2=sreed_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_f_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(out), dimension(:,:,:) :: buf
  srw_f_3=sreed_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function

integer function srw_f_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(out), dimension(:,:,:,:,:) :: buf
  srw_f_5=sreed_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_f_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(out), dimension(:,:,:,:) :: buf
  srw_f_4=sreed_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function

integer function srw_c_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(out), dimension(:) :: buf
  srw_c_1=sreed_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_c_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(out), dimension(:,:) :: buf
  srw_c_2=sreed_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_c_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(out), dimension(:,:,:) :: buf
  srw_c_3=sreed_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_c_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(out), dimension(:,:,:,:,:) :: buf
  srw_c_5=sreed_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function srw_c_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(out), dimension(:,:,:,:) :: buf
  srw_c_4=sreed_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function



!srite window
integer function sww_i_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(in), dimension(:) :: buf
  sww_i_1=srite_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_i_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(in), dimension(:,:) :: buf
  sww_i_2=srite_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_i_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(in), dimension(:,:,:) :: buf
  sww_i_3=srite_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_i_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(in), dimension(:,:,:,:) :: buf
  sww_i_4=srite_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_i_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  integer,intent(in), dimension(:,:,:,:,:) :: buf
  sww_i_5=srite_window_i(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_f_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(in), dimension(:) :: buf
  sww_f_1=srite_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_f_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(in), dimension(:,:) :: buf
  sww_f_2=srite_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_f_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(in), dimension(:,:,:) :: buf
  sww_f_3=srite_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_f_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(in), dimension(:,:,:,:) :: buf
  sww_f_4=srite_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_f_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  real,intent(in), dimension(:,:,:,:,:) :: buf
  sww_f_5=srite_window_f(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_c_1(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(in), dimension(:) :: buf
  sww_c_1=srite_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_c_2(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(in), dimension(:,:) :: buf
  sww_c_2=srite_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_c_3(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(in), dimension(:,:,:) :: buf
  sww_c_3=srite_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_c_4(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(in), dimension(:,:,:,:) :: buf
  sww_c_4=srite_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function
integer function sww_c_5(tag,nd,ng,n,f,j,sz,buf)
  character(len=*),intent(in) :: tag
  integer,intent(in) :: nd,sz
  integer,intent(in),dimension(:) :: n,f,j,ng
  complex,intent(in), dimension(:,:,:,:,:) :: buf
  sww_c_5=srite_window_c(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,sz,buf)
end function







end module
