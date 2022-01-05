module sep_3d_mod
use, intrinsic :: iso_c_binding
use sep_func_mod
implicit none



interface 
  integer function sep_put_keyf(tag,nm,typ,val,index) bind(c,name="sep_put_key")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR),dimension(*),intent(in) :: nm,typ,val
    integer(C_INT),intent(in) :: index
  end function

  integer function sep_get_grid_windowf(tag,nd,ng,n,f,j,ar) bind(c,name="sep_get_grid_window")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
    integer(C_INT),intent(in) :: ng(*),nd,n(*),f(*),j(*)
    integer(C_INT),intent(out) :: ar(*)
  end function
  integer function sep_put_grid_windowf(tag,nd,ng,n,f,j,ar) bind(c,name="sep_put_grid_window")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
    integer(C_INT),intent(in) :: nd,ng(*),n(*),f(*),j(*)
    integer(C_INT),intent(in) :: ar(*)
  end function
  integer function sep_get_key_typef(tag,ikey,nm) bind(c,name="sep_get_key_type")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
        character(C_CHAR), dimension(*),intent(out) :: nm

    integer(C_INT),intent(in) :: ikey
  end function
  integer function sep_get_key_indexf(tag,nm,ikey) bind(c,name="sep_get_key_index")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag,nm
    integer(C_INT),intent(out) :: ikey
  end function
    integer function sep_get_number_keysf(tag,ikey) bind(c,name="sep_get_number_keys")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
    integer(C_INT),intent(out) :: ikey
  end function
    integer function sep_put_number_keysf(tag,ikey) bind(c,name="sep_put_number_keys")
    import
    character(C_CHAR), dimension(*),intent(in) :: tag
    integer(C_INT),intent(in) :: ikey
  end function

  integer function sep_put_val_headers_if(tag,n,f,ar) bind(c,name="sep_put_val_headers_if")
    import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT),intent(in) :: n,f
    integer(C_INT),dimension(*),intent(in) :: ar
  end function
    integer function sep_put_val_headers_ff(tag,n,f,ar) bind(c,name="sep_put_val_headers_ff")
    import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT),intent(in) :: n,f
    real(C_FLOAT),dimension(*),intent(in) :: ar
  end function
    integer function sep_get_val_headers_if(tag,n,f,ar) bind(c,name="sep_get_val_headers_if")
    import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT),intent(in) :: n,f
    integer(C_INT),dimension(*),intent(out) :: ar
  end function
      integer function sep_get_val_headers_ff(tag,n,f,ar) bind(c,name="sep_get_val_headers_ff")
          import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT),intent(in) :: n,f
    real(C_FLOAT),dimension(*),intent(out) :: ar
  end function
  
  


  integer function sep_tag_is_pipef(tag1) bind(c,name="sep_tag_is_pipe")
        import
character(C_CHAR), dimension(*) :: tag1
  end function

  integer function sep_copy_data_pointerf(tag1,tag2) bind(c,name="sep_copy_data_pointer")
        import
character(C_CHAR), dimension(*) :: tag1,tag2
  end function


  integer function sep_copy_header_keysf(tag1,tag2) bind(c,name="sep_copy_header_keys")
       import
 character(C_CHAR), dimension(*) :: tag1,tag2
  end function
  

  integer function sep_get_number_grid_axesf(tag,nd) bind(c,name="sep_get_number_grid_axes")
    import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT) :: nsz
  end function

  integer function sep_get_number_header_axesf(tag,nd) bind(c,name="sep_get_number_header_axes")
    import
    character(C_CHAR), dimension(*) :: tag
    integer(C_INT) :: nsz
  end function

  subroutine sep_copy_gfff(tag1,tag2) bind(c,name="sep_copy_gff")
   import
    character(C_CHAR), dimension(*) :: tag1,tag2
  end subroutine

  subroutine sep_set_no_headersf(tag) bind(c,name="sep_set_no_headers")
       import
 character(C_CHAR), dimension(*) :: tag
  end subroutine
  
  subroutine sep_set_regular_gridf(tag) bind(c,name="sep_set_regular_grid")
       import
 character(C_CHAR), dimension(*) :: tag
  end subroutine
  
  
    subroutine sep_set_no_gridf(tag) bind(c,name="sep_set_no_grid")
       import
 character(C_CHAR), dimension(*) :: tag
  end subroutine
  
integer function sep_get_header_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_get_header_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(out) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(out) :: n
      real(C_FLOAT),intent(out) :: o,d
  end  function
   
   integer function sep_get_grid_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_get_grid_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(out) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(out) :: n
      real(C_FLOAT),intent(out) :: o,d
  end  function
   
   
   integer function sep_put_header_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_put_header_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(in) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(in) :: n
      real(C_FLOAT),intent(in) :: o,d
  end  function
  integer function sep_reorder_dataf(t1,t2,n1,n2,ord) bind(c,name="sep_reorder_data")
    import
    character(C_CHAR),dimension(*),intent(in) :: t1,t2
    integer(C_INT),intent(in),value :: n1,n2
    integer(C_INT),intent(in) ::ord(*)
  end function
   integer function sep_put_grid_axis_parf(tag,iax,n,o,d,label) bind(c,name="sep_put_grid_axis_par")
        import
  character(C_CHAR), dimension(*),intent(in) :: tag
    character(C_CHAR), dimension(*),intent(in) :: label
    integer(C_INT),intent(in) :: iax
    integer(C_INT),intent(in) :: n
      real(C_FLOAT),intent(in) :: o,d
  end  function
   
  
  
  
    subroutine init_3d() bind(c,name="init_3d")
    import
  end subroutine

     subroutine sep_3d_close() bind(c,name="sep_3d_close")
    import
  end subroutine

end interface
  interface sep_put_val_headers
    module procedure sep_put_val_headers_i1,sep_put_val_headers_i2,putvhf
  end interface
  interface sep_get_val_headers
    module procedure getheadi1,getheadi2,getheadf1,getheadf2
  end interface
   interface sep_put_grid_window
     module procedure sputgridi1,sputgridi2,sputgridi3,sputgridi4,sputgridi5,sputgridi6,sputgrid7
  end interface
  
  interface sep_get_grid_window
    module procedure getwind1,getwind2,getwind3,getwind4
  end interface
  


contains
  integer function sep_copy_header_keys(tag1,tag2)
    character(len=*) ::tag1,tag2
    sep_copy_header_keys= sep_copy_header_keysf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR)
  end function
  
    integer function sep_copy_data_pointer(tag1,tag2)
    character(len=*) ::tag1,tag2
    sep_copy_data_pointer= sep_copy_data_pointerf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR)
  end function
  
  subroutine sep_copy_gff(tag1,tag2)
    character(len=*) ::tag1,tag2
    call sep_copy_gfff(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR)
  end subroutine
  subroutine sep_set_no_headers(tag)
    character(len=*) ::tag
    call sep_set_no_headersf(trim(tag)//C_NULL_CHAR)
  end subroutine
  subroutine sep_set_regular_grid(tag)
    character(len=*) ::tag
    call sep_set_regular_gridf(trim(tag)//C_NULL_CHAR)
  end subroutine
  subroutine sep_set_no_grid(tag)
    character(len=*) ::tag
    call sep_set_no_gridf(trim(tag)//C_NULL_CHAR)
  end subroutine

  integer function sep_tag_is_pipe(tag)
    character(len=*) :: tag
    sep_tag_is_pipe=sep_tag_is_pipef(trim(tag)//C_NULL_CHAR)
  end function
  integer function sep_get_number_grid_axes(tag,nd)
    character(len=*) :: tag
    integer nd
    sep_get_number_grid_axes=sep_get_number_grid_axesf(tag,nd)
  end function

  integer function sep_get_number_header_axes(tag,nd)
    character(len=*) :: tag
    integer  :: nd
    sep_get_number_header_axes=sep_get_number_header_axesf(trim(tag)//C_NULL_CHAR,nd)
  end function
  integer function sep_put_val_headers_i1(tag,nr,nh,vals)
    character(len=*) :: tag
    integer :: nr,nh,vals(:)
    sep_put_val_headers_i1=sep_put_val_headers_if(trim(tag)//C_NULL_CHAR,nr,nh,vals)
  end function
  integer function putvhf(tag,nr,nh,vals)
    character(len=*) :: tag
    integer :: nr,nh
    real :: vals(:)
    putvhf=sep_put_val_headers_ff(trim(tag)//C_NULL_CHAR,nr,nh,vals)
  end function
  integer function sep_put_val_headers_i2(tag,nr,nh,vals)
    character(len=*) :: tag
    integer :: nr,nh,vals(:,:)
    sep_put_val_headers_i2=sep_put_val_headers_if(trim(tag)//C_NULL_CHAR,nr,nh,vals)
  end function
  
  integer function sputgridi1(tag ,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:)
    sputgridi1=sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgridi2(tag ,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:)
    sputgridi2=sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgridi3(tag ,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:,:)
    sputgridi3=sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgridi4(tag,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:,:,:)
    sputgridi4=sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgridi5(tag,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:,:,:,:)
    sputgridi5= sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgridi6(tag,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:,:,:,:,:)
     sputgridi6= sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sputgrid7(tag,nd,ng,n,f,j,vals)
    character(len=*) :: tag
    integer :: nd,ng(:),n(:),f(:),j(:),vals(:,:,:,:,:,:,:)
     sputgrid7= sep_put_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  integer function sep_get_grid_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(out) :: n
    real,intent(out) :: o,d
    character(len=*) :: label
     sep_get_grid_axis_par= sep_get_grid_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,label)
     call c2forstr(label)
  end function

    integer function sep_get_header_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(out) :: n
    real,intent(out) :: o,d
    character(len=*) :: label
     sep_get_header_axis_par= sep_get_header_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,label)
     call c2forstr(label)
  end function

  integer function sep_put_grid_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(in) :: n
    real,intent(in) :: o,d
    character(len=*),intent(in) :: label
     sep_put_grid_axis_par= sep_put_grid_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,trim(label)//C_NULL_CHAR)
  end function

    integer function sep_put_header_axis_par(tag,iax,n,o,d,label)
    character(len=*),intent(in) :: tag
    integer,intent(in) ::iax
    integer,intent(in) :: n
    real,intent(in) :: o,d
    character(len=*),intent(in):: label
     sep_put_header_axis_par= sep_put_header_axis_parf(trim(tag)//C_NULL_CHAR,iax,n,o,d,trim(label)//C_NULL_CHAR)
  end function
  
    integer function sep_get_key_index(tag,nm,ikey) 
    character(len=*),intent(in) :: tag,nm
    integer,intent(out) :: ikey
      sep_get_key_index=sep_get_key_indexf(trim(tag)//C_NULL_CHAR,trim(nm)//C_NULL_CHAR,ikey)
  end function
    integer function sep_put_number_keys(tag,ikey) 
    character(len=*),intent(in) :: tag
    integer,intent(in) :: ikey
      sep_put_number_keys=sep_put_number_keysf(trim(tag)//C_NULL_CHAR,ikey)
  end function
      integer function sep_get_number_keys(tag,ikey) 
    character(len=*),intent(in) :: tag
    integer,intent(out) :: ikey
      sep_get_number_keys=sep_get_number_keysf(trim(tag)//C_NULL_CHAR,ikey)
  end function
  integer function getheadi1(tag,f,n,head)
    character(len=*) :: tag
    integer,intent(in) :: f,n
    integer,dimension(:),intent(out) :: head
    getheadi1=sep_get_val_headers_if(trim(tag)//C_NULL_CHAR,f,n,head)
  end function
    integer function getheadi2(tag,f,n,head)
    character(len=*) :: tag
    integer,intent(in) :: f,n
    integer,dimension(:,:),intent(out) :: head
    getheadi2=sep_get_val_headers_if(trim(tag)//C_NULL_CHAR,f,n,head)
  end function
    integer function getheadf1(tag,f,n,head)
    character(len=*) :: tag
    integer,intent(in) :: f,n
    real,dimension(:),intent(out) :: head
    getheadf1=sep_get_val_headers_ff(trim(tag)//C_NULL_CHAR,f,n,head)
  end function
    integer function getheadf2(tag,f,n,head)
    character(len=*) :: tag
    integer,intent(in) :: f,n
    real,dimension(:,:),intent(out) :: head
    getheadf2=sep_get_val_headers_ff(trim(tag)//C_NULL_CHAR,f,n,head)
  end function
  integer function sep_get_key_type(tag,index,typ)
    character(len=*),intent(in) :: tag
    integer,intent(in) :: index
    character(len=*),intent(out) :: typ
    sep_get_key_type=sep_get_key_typef(trim(tag)//C_NULL_CHAR,index,typ)
    call c2forstr(typ)
  end function
  integer function sep_reorder_data(t1,t2,n1,n2,ord)
    character(len=*),intent(in) :: t1,t2
    integer,intent(in) :: n1,n2,ord(:)
    sep_reorder_data=sep_reorder_dataf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n1,n2,ord)
  end function
  

    integer function getwind1(tag,nd,ng,n,f,j,vals)
    character(len=*), intent(in) ::tag
    integer,intent(in) :: nd,ng(:),n(:),f(:),j(:)
    integer,intent(out)  :: vals(:)
    getwind1=sep_get_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
    integer function getwind2(tag,nd,ng,n,f,j,vals)
    character(len=*), intent(in) ::tag
    integer,intent(in) :: nd,ng(:),n(:),f(:),j(:)
    integer,intent(out)  :: vals(:,:)
    getwind2=sep_get_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
    integer function getwind3(tag,nd,ng,n,f,j,vals)
    character(len=*), intent(in) :: tag
    integer,intent(in) :: nd,ng(:),n(:),f(:),j(:)
    integer,intent(out)  :: vals(:,:,:)
    getwind3=sep_get_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
    integer function getwind4(tag,nd,ng,n,f,j,vals)
    character(len=*), intent(in) :: tag
    integer,intent(in) :: nd,ng(:),n(:),f(:),j(:)
    integer,intent(out)  :: vals(:,:,:,:)
    getwind4=sep_get_grid_windowf(trim(tag)//C_NULL_CHAR,nd,ng,n,f,j,vals)
  end function
  
  
  
        integer function sep_put_key(tag,nm,typ,val,ind)
    character(len=*), intent(in) :: tag
    integer,intent(in) :: ind
    character(len=*),intent(in)  :: nm,typ,val
    sep_put_key=sep_put_keyf(trim(tag)//C_NULL_CHAR,trim(nm)//C_NULL_CHAR,&
      trim(typ)//C_NULL_CHAR,trim(val)//C_NULL_CHAR,ind)
  end function
end module

