!
!sep3d_struct_mod
module sep3d_struct_mod 
  use sep_3d_mod
  use sep_mod
  use, intrinsic :: iso_c_binding
  implicit none
  integer, private, save :: counter
  
  
  interface
  
       integer function sep3d_set_header_vals_sf(t1,t2,v) bind(c,name="sep3d_set_header_vals_sf")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     real(C_FLOAT),dimension(*),intent(in) :: v
    end function
           integer function sep3d_close_tagf(t1,t2) bind(c,name="sep3d_close_tag")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
    end function
               integer function sep3d_deletef(t1) bind(c,name="sep3d_delete")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1
    end function
         integer function sep3d_set_header_vals_si(t1,t2,v) bind(c,name="sep3d_set_header_vals_si")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(in) :: v
    end function
  
       integer function sep3d_set_header_vals_if(t1,t2,v) bind(c,name="sep3d_set_header_vals_if")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1
          integer(C_INT) ,value :: t2
     real(C_FLOAT),dimension(*),intent(in) :: v
    end function
    
    
         integer function sep3d_set_header_vals_ii(t1,t2,v) bind(c,name="sep3d_set_header_vals_ii")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1
     integer(C_INT) ,value :: t2
     integer(C_INT),dimension(*),intent(out) :: v
    end function
  
     integer function sep3d_grab_header_vals_sf(t1,t2,v) bind(c,name="sep3d_grab_header_vals_sf")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     real(C_FLOAT),dimension(*),intent(out) :: v
    end function
    
    
         integer function sep3d_grab_header_vals_si(t1,t2,v) bind(c,name="sep3d_grab_header_vals_si")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(out) :: v
    end function
  
       integer function sep3d_grab_header_vals_if(t1,t2,v) bind(c,name="sep3d_grab_header_vals_if")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1
          integer(C_INT) ,value :: t2
     real(C_FLOAT),dimension(*),intent(out) :: v
    end function
    
    
         integer function sep3d_grab_header_vals_ii(t1,t2,v) bind(c,name="sep3d_grab_header_vals_ii")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1
     integer(C_INT) ,value :: t2
     integer(C_INT),dimension(*),intent(in) :: v
    end function
  
   integer function sep3d_riteff(t1,t2,n,f,j,v,i1,i2,i3,i4) bind(c,name="sep3d_riteff")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(in) :: n,f,j
     real(C_FLOAT),dimension(*),intent(in) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3,i4
    end function
       integer function sep3d_riteif(t1,t2,n,f,j,v,i1,i2,i3,i4) bind(c,name="sep3d_riteif")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(in) :: n,f,j
     integer(C_INT),dimension(*),intent(in) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3,i4
    end function
       integer function sep3d_ritecf(t1,t2,n,f,j,v,i1,i2,i3,i4) bind(c,name="sep3d_ritecf")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(in) :: n,f,j
     complex(C_FLOAT),dimension(*),intent(in) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3,i4
    end function
           integer function sep3d_reedcf(t1,t2,i1,i2,i3,v) bind(c,name="sep3d_read_datacf")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     complex(C_FLOAT),dimension(*),intent(out) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3
    end function
    
               integer function sep3d_reedif(t1,t2,i1,i2,i3,v) bind(c,name="sep3d_read_dataif")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     integer(C_INT),dimension(*),intent(out) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3
    end function
               integer function sep3d_reedff(t1,t2,i1,i2,i3,v) bind(c,name="sep3d_read_dataff")
     import
     character(C_CHAR),dimension(*) ,intent(in) :: t1,t2
     real(C_FLOAT),dimension(*),intent(out) :: v
     integer(C_INT),intent(in),value :: i1,i2,i3
    end function
     
             integer function sep3d_grab_key_indexf(in,lab,val) bind(c,name="sep3d_grab_key_index")
    import
    character(C_CHAR),dimension(*),intent(in) :: in,lab
    integer(C_INT),intent(out) :: val
 end function
              integer function sep3d_rite_statusf(in,val,val2) bind(c,name="sep3d_set_rite_status")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in) ,value :: val,val2
 end function
           integer function sep3d_grab_axis_indexf(in,lab,val) bind(c,name="sep3d_grab_axis_index")
    import
    character(C_CHAR),dimension(*),intent(in) :: in,lab
    integer(C_INT),intent(out) :: val
 end function
         integer function sep3d_grab_axisf(in,index,n,o,d,label,unit) bind(c,name="sep3d_grab_axis")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: n
    integer(C_INT),intent(in),value :: index
    real(C_FLOAT),intent(out) :: o,d
    character,dimension(*),intent(out) :: label,unit
 end function
 
          integer function sep3d_set_axisf(in,index,n,o,d,label,unit) bind(c,name="sep3d_set_axis")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: n
    integer(C_INT),intent(in),value :: index
    real(C_FLOAT),intent(in),value :: o,d
    character,dimension(*),intent(in) :: label,unit
 end function
 
          integer function sep3d_change_dimsf(in,ndim,n) bind(c,name="sep3d_change_dims")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in)  :: n(*)
    integer(C_INT),intent(in),value :: ndim
 end function
           integer function sep3d_set_ngf(in,ndim) bind(c,name="sep3d_set_ng")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ndim
 end function
            integer function sep3d_grab_ngf(in,ndim) bind(c,name="sep3d_grab_ng")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ndim
 end function
        integer function sep3d_grab_coord_valsf(in,ind) bind(c,name="sep3d_grab_coord_vals")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out),dimension(*) :: ind
 end function
         integer function sep3d_inorderf(in) bind(c,name="sep3d_set_inorder")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
 end function
          integer function sep3d_print_infof(in) bind(c,name="sep3d_print_info")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
 end function
         integer function sep3d_set_coord_valsf(in,ind) bind(c,name="sep3d_set_coord_vals")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),dimension(*) :: ind
 end function
       integer function sep3d_grab_windf(in,n,f,j) bind(c,name="sep3d_grab_wind")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out),dimension(*) :: n,f,j
 end function
        integer function sep3d_set_windf(in,n,f,j) bind(c,name="sep3d_set_wind")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out),dimension(*) :: n,f,j
 end function
 
         integer function sep3d_read_headersf(t1,t2,n,f,j,v) bind(c,name="sep3d_read_headers")
    import
    character(C_CHAR),dimension(*),intent(in) :: t1,t2
    integer(C_INT),intent(in),dimension(*) :: n,f,j
    integer(C_INT),intent(out) :: v
 end function
 
            integer function sep3d_grab_keyf(in,index,nm,typ,fmt) bind(c,name="sep3d_grab_key")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(out) :: nm,typ,fmt
    integer(C_INT),intent(in),value :: index
 end function
             integer function sep3d_set_keyf(in,index,nm,typ,fmt) bind(c,name="sep3d_set_key")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: nm,typ,fmt
    integer(C_INT),intent(in),value :: index
 end function
          integer function sep3d_grab_usagef(in,typ) bind(c,name="sep3d_grab_usage")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(out) :: typ
 end function
          integer function sep3d_rite_formatf(in,typ) bind(c,name="sep3d_rite_format")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: typ
 end function
 
           integer function sep3d_rite_ntracesf(in,fl) bind(c,name="sep3d_rite_ntraces")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
             integer function sep3d_ge_spacef(in,fl) bind(c,name="sep3d_ge_space")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
              integer function sep3d_copy_gridf(in,fl) bind(c,name="sep3d_copy_grid")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
               integer function sep3d_copy_headersf(in,fl) bind(c,name="sep3d_copy_headers")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
                integer function sep3d_grab_grid_blockf(in,val) bind(c,name="sep3d_grab_grid_block")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),dimension(*),intent(out) :: val
 end function
                 integer function sep3d_set_grid_valsf(in,val) bind(c,name="sep3d_set_grid_vals")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),dimension(*),intent(in) :: val
 end function
 
 
               integer function sep3d_copy_coordsf(in,fl) bind(c,name="sep3d_copy_coords")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
            integer function sep3d_conformf(in,fl) bind(c,name="sep3d_conform")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
            integer function sep3d_copy_structf(in,fl) bind(c,name="sep3d_copy_struct")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: fl
 end function
 
         integer function sep_num_thread() bind(c,name="sep_num_thread")
    import
 end function
          integer function sep_thread_num() bind(c,name="sep_thread_num")
    import
 end function
         integer function sep3d_count_ntracesf(in) bind(c,name="sep3d_count_ntraces")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
 end function
        integer function sep3d_grab_file_typef(in,typ) bind(c,name="sep3d_grab_file_type")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(out) :: typ
 end function
         integer function sep3d_set_file_typef(in,typ) bind(c,name="sep3d_set_file_type")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: typ
 end function
          integer function sep3d_set_data_typef(in,typ) bind(c,name="sep3d_set_data_type")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(in) :: typ
 end function
      integer function sep3d_grab_data_typef(in,typ) bind(c,name="sep3d_grab_data_type")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    character(C_CHAR),dimension(*),intent(out) :: typ
 end function
      integer function sep3d_grab_ndimsf(in,ntr) bind(c,name="sep3d_grab_ndims")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ntr
 end function
     integer function sep3d_grab_nkeysf(in,ntr) bind(c,name="sep3d_grab_nkeys")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ntr
    end function
         integer function sep3d_set_nkeysf(in,ntr) bind(c,name="sep3d_set_nkeys")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
    end function
         integer function sep3d_alloc_coordf(in,ntr) bind(c,name="sep3d_alloc_coord")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
 end function
    integer function sep3d_grab_nhf(in,ntr) bind(c,name="sep3d_grab_nh")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ntr
 end function
     integer function sep3d_grab_ncoordf(in,ntr) bind(c,name="sep3d_grab_ncoord")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ntr
 end function
     integer function sep3d_set_ndimsf(in,ntr) bind(c,name="sep3d_set_ndims")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
 end function
      integer function sep3d_set_nhf(in,ntr) bind(c,name="sep3d_set_nh")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
 end function
      integer function sep3d_set_ntracesf(in,ntr) bind(c,name="sep3d_set_ntraces")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
 end function
       integer function sep3d_grab_drnf(in,ntr) bind(c,name="sep3d_grab_drn")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out),dimension(*) :: ntr
 end function
       integer function sep3d_set_drnf(in,ntr) bind(c,name="sep3d_set_drn")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),dimension(*) :: ntr
 end function
      integer function sep3d_fileexistf(in,ntr) bind(c,name="sep3d_fileexist")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(in),value :: ntr
 end function
    integer function sep3d_grab_ntracesf(in,ntr) bind(c,name="sep3d_grab_ntraces")
    import
    character(C_CHAR),dimension(*),intent(in) :: in
    integer(C_INT),intent(out) :: ntr
 end function
  integer function sep3d_par_initf(in,out) bind(c,name="sep3d_par_init")
    import
    character(C_CHAR),dimension(*),intent(in) :: in,out
 end function
    integer function sep3d_tag_initf(tag1,tag2,usage) bind(c,name="sep3d_tag_init")
    import
    character(C_CHAR),dimension(*),intent(in) :: tag1,tag2,usage
 end function
     integer function sep3d_tag_init_threadf(tag1,tag2,usage,ind) bind(c,name="sep3d_tag_init_thread")
    import
    character(C_CHAR),dimension(*),intent(in) :: tag1,tag2,usage
    integer(C_INT), value :: ind
 end function
     integer function sep3d_struct_initf(tag1,tag2,typ) bind(c,name="sep3d_struct_init")
    import
    character(C_CHAR),dimension(*),intent(in):: tag1,tag2,typ
 end function
  end interface
  

  integer, external, private  :: sep3d_read_list

!!$=head1 NAME
!!$
!!$sep3df - sep3d f90 structure (superset)
!!$
!!$=head1 SYNOPSIS
!!$
!!$type(sep3d)
!!$
!!$=head1 DESCRIPTION
!!$
!!!$type sep3d{
!!$
!!$        character(len=128) :: tag              #tag to sep3df, MUST  BE SPECIFIED
!!$
!!$        integer, dimension(:),pointer  :: n #length of regular axes
!!$
!!$        real, dimension(:) ,pointer    :: o,d #initial and sampling  of regular axes
!!$
!!$        character(len=128),dimension(:),pointer::label,unit #label and unit for axes
!!$
!!$        integer :: nkeys,ndims           #number of keys and regular dimensions
!!$
!!$        integer :: drn #data record number (-1 if in same order)
!!$
!!$        integer :: ntraces # number of traces in the dataset
!!$
!!$
!!$        character(len=128),pointer,dimension(:) :: keyname,keytype,keyfmt #keyname,type and format
!!$
!!$        character(len=128) :: usage #usage (INPUT,OUTPUT,SCRATCH) for dataset
!!$
!!$        character(len=128) :: data_format #data format (FLOAT,BYTE,COMPLEX,INTEGER)
!!$
!!$        character(len=128) :: file_format #file format (REGULAR,HEADER,GRID)
!!$
!!$
!!$
!!$
!!$
!!$
!!$=head1 COMMENTS
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
  type sep3d 
    character(len=128) :: tag              
    !tag to sep3df, MUST  BE SPECIFIED
    integer, dimension(:),pointer  :: n !length of regular axes
    integer, dimension(:),pointer  :: nwind,fwind,jwind 
    !window parameters
    real, dimension(:) ,pointer    :: o,d 
    !initial and sampling  of regular axes
    character(len=128),dimension(:),pointer::label,unit 
    !label and unit for axes
    integer :: nkeys,ndims           !number of keys and regular dimensions
    integer :: drn !data record number (-1 if in same order)
    integer :: ntraces ! number of traces in the dataset
    character(len=128),pointer,dimension(:) :: keyname,keytype,keyfmt
!keyname,type and format
    character(len=128) :: usage 
    !usage (INPUT,OUTPUT,SCRATCH) for dataset
    character(len=128) :: data_format 
    !data format (FLOAT,BYTE,COMPLEX,INTEGER)
    character(len=128) :: file_format 
    !file format (REGULAR,HEADER,GRID)
  end type
  interface init_sep3d
  module procedure tag_init,struct_init,par_init
  end interface 
!!$
!!$=head1 NAME
!!$
!!$ init_sep3d - initialize a SEP3d type
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call init_sep3d(instruct, outstruct, usage)>
!!$
!!$C<call init_sep3d(intag,outstruct,usage)>
!!$
!!$C<call init_sep3d(outstruct,usage,data_format,file_format,n,o,d,label,unit,ntraces,keyname,keytype,keyfmt,nh)>
!!$
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item IMPORTANT FOR FORTRAN USERS
!!$    include C<ctag="out"> as an argument when initializing output.
!!$    (have you noticed a f90 file mysteriously in your directory?)
!!$
!!$=item instruct  -  sep3d
!!$
!!$      Sep3d to structure to copy from
!!$
!!$=item usage     -  char*
!!$
!!$      Usage of the output tag
!!$
!!$=item tag       -  char*
!!$
!!$      Tag to initilize from
!!$
!!$=item data_format- char*
!!$
!!$      Data format (FLOAT,INTEGER,BYTE,COMPLEX)
!!$
!!$=item file_format- char*
!!$
!!$      File format (GRID,HEADER,REGULAR)
!!$
!!$=item n-  int*
!!$
!!$      (optional) axis dimensions
!!$
!!$=item o-  float*
!!$
!!$     (optional)  first sample axis
!!$
!!$=item d-  float*
!!$
!!$      (optional)  sampling of axis
!!$
!!$=item label-  char*
!!$
!!$      (optional)  label for axis
!!$
!!$=item unit-  char*
!!$
!!$      (optional)  unit for axis
!!$
!!$=item ntraces   -  int
!!$
!!$      (optional)  number of trac3es
!!$
!!$=item keyname  - char**
!!$
!!$      (optional)  keyname for dataset
!!$
!!$=item keytype  -  char**
!!$
!!!$      (optional)  keytype for dataset
!!$
!!$=item keyfmt -  char**
!!!$
!!$      (optional)  keyfmt for dataset
!!$
!!$=item nh-  int
!!$
!!!$      (optional)  number of headers to store
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!!$
!!$=item outstruct -  sep3d
!!$
!!$      sep3d structure to create
!!$
!!$=back
!!$
!!!$
!!$=head1 DESCRIPTION
!!$
!!$Initialize structure
!!$
!!$
!!$=head1 LIBRARY
!!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
!!$=head1 NAME
!!$
!!$sep3d_grab_key_vals - grab header values from C structure
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_grab_key_vals(struct, locate, values)>
!!$   must have already made a call to C<sep3d_grab_headers('in',input,nh)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!!$=item struct   -    sep3d
!!$
!!$      structure to grab from
!!$
!!$=item locate   -    (char*/int)
!!$
!!$      locate through keyname or keyindex
!!$
!!$=back
!!$
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!!$=over 4
!!$
!!$=item values   -    (float*/int*)
!!$
!!$              array of header values
!!!$
!!$=back
!!$
!!!$=head1 DESCRIPTION
!!$
!!$Grab header values
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_set_key_vals>
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
  interface sep3d_grab_key_vals
  module procedure getkey_i_c,getkey_i_i,getkey_f_c,getkey_f_i
  end interface
!!$=head1 NAME
!!$
!!$sep3d_set_key_vals  - set header values in C structure
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_set_key_vals(struct, locate, values)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct   -    sep3d
!!$
!!$      structure to grab from
!!$
!!!$=item locate   -    (char*/int)
!!$
!!$      locate through keyname or keyindex
!!$
!!!$=item values   -    (float*/int*)
!!$
!!$      array of header values
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Set header values
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_grab_key_vals>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  interface sep3d_set_key_vals
  module procedure putkey_i_c,putkey_i_i,putkey_f_c,putkey_f_i
  end interface
!!$
!!$=head1 NAME
!!$
!!$sep3d_read_data -read in data
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic=sep3d_read_data(tag,struct,data,fwind,jwind,nwind)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag      -  char*
!!$
!!$      name of the tag to read from
!!$
!!$=item struct   -  sep3d
!!$
!!$      structure to read from
!!$
!!$=item fwind    -  int*
!!$
!!$      (optional) begining of window to read
!!$
!!$=item jwind    -  int*
!!$
!!$      (optional) sampling to read
!!$
!!$=item nwind    -  int*
!!$
!!$      (optional) number of elem to read
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item data     -  float/int/complex
!!$
!!$       data to read
!!$
!!$=back
!!$
!!$=head1 RETURN VALUES
!!$
!!$=over 4
!!$
!!$=item true    -
!!$
!!$      if success
!!$
!!$=item false   -
!!$
!!$      if it fails
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Read data
!!$
!!$
!!$=head1 COMMENTS
!!$
!!$Must read in the headers first
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_write_data>,L<sep3d_grab_headers>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  interface sep3d_read_data
  module procedure  get_data_2i,get_data_2c,get_data_2f
  module procedure  get_data_3i,get_data_3c,get_data_3f
  module procedure  get_data_4i,get_data_4c,get_data_4f
  end interface
!!$
!!$=head1 NAME
!!$
!!$sep3d_write_data - write out data
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic=sep3d_write_data(tag,struct,data,fwind,jwind,nwind,write_headers,write_grid)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag      -  char*
!!$
!!$      name of the tag to read from
!!$
!!$=item struct   -  sep3d
!!$
!!$      structure to read from
!!$
!!$=item fwind    -  int*
!!$
!!$      (optional) begining of window to read
!!$
!!$=item jwind    -  int*
!!$
!!$      (optional) sampling to read
!!$
!!$=item nwind    -  int*
!!$
!!$      (optional) number of elem to read
!!$
!!$=item data     -  float/int/complex
!!$
!!$      data tor read
!!$
!!$=item write_headers - logical
!!$
!!$      wheter or not write headers
!!$
!!$=item write_grid - logical
!!$
!!$      wheter or not write grid
!!$
!!$=back
!!$
!!$=head1 RETURN VALUES
!!$
!!$=over 4
!!$
!!$=item true    -
!!$
!!$      if success
!!$
!!$=item false   -
!!$
!!$      if fails
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Write data
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_read_data>,L<sep3d_grab_headers>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  interface sep3d_reed_data
    module procedure read_f_1,read_f_2,read_f_3,read_f_4,read_f_5
    module procedure read_c_1,read_c_2,read_c_3,read_c_4,read_c_5
        module procedure read_i_1,read_i_2,read_i_3,read_i_4,read_i_5

end interface
interface sep3d_rite
  module procedure rite_f_1,rite_f_2,rite_f_3,rite_f_4,rite_f_5
  module procedure rite_c_1,rite_c_2,rite_c_3,rite_c_4,rite_c_5
    module procedure rite_i_1,rite_i_2,rite_i_3,rite_i_4,rite_i_5

end interface

  

  interface sep3d_write_data
  module procedure  put_data_2i,put_data_2c,put_data_2f,put_data_n
  module procedure  put_data_3i,put_data_3c,put_data_4f
  module procedure  put_data_4i,put_data_4c,put_data_3f
  end interface
!
!
!
  contains
!!$
!!$=head1 NAME
!!$
!!$ sep3d_grab_sep3d - synchronize f90 structure with C structure
!!$
!!$=head1 SYNOPSIS
!!$
!!$ call sep3d_grab_sep3d(sep3dc,sep3df)
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item sep3dc  -  char*
!!$
!!$      pointer to C sep3d structure
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item sep3df  -  sep3d
!!$
!!$      fortran sep3d structure to copy to
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Copies sep3d C structure to its fortran equivilant
!!$
!!$=head1 SEE ALSO
!!$
!!$ sep3d_set_sep3d
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_grab_sep3d(sep3dc,sep3df)
    type(sep3d) :: sep3df
    character(len=*),optional    :: sep3dc
    character(len=128) :: tag
    integer     :: ierr,i1
    if (valid_structure(sep3df) .and. (sep3df%ndims.ne.0 .or. sep3df%nkeys.ne.0&
      &)) then
      call sep3d_clean(sep3df)
      if (.not. present(sep3dc)) then
        tag=sep3df%tag
      else
        tag =sep3dc
      end if
    else if (.not. present(sep3dc)) then
      call seperr("must specify a tag when first reading a structure")
    else
      tag=sep3dc(1:MIN(LEN(sep3dc),LEN(tag)))
    end if 
    ierr= sep3d_grab_ndims(tag,sep3df%ndims)
    if (ierr.ne.0) then
      call seperr("trouble grabing ndim for tag ")
    end if
    if (sep3df%ndims>0) then
      allocate(sep3df%n(sep3df%ndims), sep3df%o(sep3df%ndims), sep3df%d&
        &(sep3df%ndims))
      allocate(sep3df%nwind(sep3df%ndims), sep3df%fwind(sep3df%ndims),&
        & sep3df%jwind(sep3df%ndims))
      allocate(sep3df%label(sep3df%ndims), sep3df%unit(sep3df%ndims)&
        &,stat=ierr)
      if (ierr.ne.0) then
        call seperr("trouble allocating n,o,d in structure")
      end if
      do i1=1,sep3df%ndims 
         ierr=sep3d_grab_axis(tag,i1,sep3df%n(i1),sep3df%o(i1),sep3df%d&
          &(i1),sep3df%label(i1),sep3df%unit(i1))
      end do 
      if (0.ne. sep3d_grab_wind(sep3dc,sep3df%nwind,sep3df%fwind&
        &,sep3df%jwind)) then
        call seperr("trouble grabbing window parameters")
      end if
    end if
    ierr=sep3d_grab_nkeys(tag,sep3df%nkeys)
    if (ierr.ne.0) then
      call seperr("trouble grabbing nkeys ")
    end if
    if (sep3df%nkeys>0) then
      allocate(sep3df%keyname(sep3df%nkeys),sep3df%keytype(sep3df%nkeys&
        &)) 
      allocate(sep3df%keyfmt(sep3df%nkeys),stat=ierr)
      if (ierr.ne.0) then
        call seperr("trouble allocating keys in sep3df")
      end if
      do i1=1,sep3df%nkeys 
        if (0.ne.sep3d_grab_key(tag,i1,sep3df%keyname(i1),sep3df%keytype&
          &(i1),sep3df%keyfmt(i1))) then
          sep3df%keyname(i1)="Unspecified"
          sep3df%keytype(i1)="Unspecified"
          sep3df%keyfmt(i1)="Unspecified"
        end if
      end do
    end if
    ierr=sep3d_grab_usage(tag,sep3df%usage)
    ierr=sep3d_grab_file_type(tag,sep3df%file_format)
    ierr=sep3d_grab_data_type(tag,sep3df%data_format)
    sep3df%tag=tag
    if (0.ne. sep3d_grab_ntraces(tag,sep3df%ntraces)) then
      sep3df%ntraces=0
    end if
  end subroutine 
!!$
!!$=head1 NAME
!!$
!!$ sep3d_set_sep3d - synchronize C structure with f90 structure
!!$
!!$=head1 SYNOPSIS
!!$
!!$<call sep3d_set_sep3d(sep3df)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item sep3df  -  sep3d
!!$
!!$      fortran sep3d structure to copy from
!!$
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Copies sep3d F90 structure to its C equivilant
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_grab_sep3d>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_set_sep3d(sep3df)
    type(sep3d) :: sep3df
    character(len=128) :: sep3dc
    integer :: i1,ierr
    if (.not. valid_structure(sep3df)) then
      call seperr("must specify usage before copy structure")
    end if
    sep3dc=sep3df%tag
    if (0.ne. sep3d_par_init(sep3dc,sep3df%usage)) then
      call seperr("trouble initializing tag")
    end if
    if (sep3df%ndims>0) then
      ierr=sep3d_set_ndims(sep3dc,sep3df%ndims)
      do i1=1,sep3df%ndims 
         ierr=sep3d_set_axis(sep3dc,i1,sep3df%n(i1),sep3df%o(i1),sep3df%d&
          &(i1),sep3df%label(i1),sep3df%unit(i1))
      end do 
      if (0.ne. sep3d_set_wind(sep3dc,sep3df%nwind,sep3df%fwind,sep3df%jwind&
        &)) then
        call seperr("trouble setting window parameters")
      end if
    end if
    if (sep3df%nkeys>0) then
      ierr= sep3d_set_nkeys(sep3dc,sep3df%nkeys)
      do i1=1,sep3df%nkeys 
        if (sep3df%keyname(i1)(1:18).ne."data_record_number") then
          ierr=sep3d_set_key(sep3dc,i1,sep3df%keyname(i1),sep3df%keytype&
            &(i1),sep3df%keyfmt(i1))
        end if
      end do
    end if
    ierr= sep3d_set_file_type(sep3dc,sep3df%file_format)
    ierr=sep3d_set_data_type(sep3dc,sep3df%data_format)
    if (0.ne.sep3df%ntraces) then
      if (0.ne. sep3d_set_ntraces(sep3dc,sep3df%ntraces)) then
        call seperr("trouble putting into C the number of traces")
      end if
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$valid_structure - check if sep3d structure is valid
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logical =valid_structure(struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct - sep3d
!!$
!!$       Structure to check the validity of
!!$
!!$=back
!!$
!!$=head1 RETURN VALUE
!!$
!!$=over 4
!!$
!!$=item .true.
!!$
!!$      if valid
!!$
!!$=item .false.
!!$
!!$      if not valid
!!$
!!$=back
!!$
!!$
!!$=head1 DESCRIPTION
!!$
!!$Check to see if structure is valid
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3df>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function valid_structure(sep3df)
    type(sep3d) :: sep3df
    valid_structure = .FALSE.    
     if(LLE(sep3df%usage,"INPUT")) valid_structure = .TRUE.
    if(LLE(sep3df%usage,"OUTPUT")) valid_structure = .TRUE.
    if(LLE(sep3df%usage,"SCRATCH")) valid_structure = .TRUE.
  end function 
  subroutine sep3df_delete(sep3df,tag1,tag2,tag3,tag4,tag5)
    type(sep3d) :: sep3df
    character(len=*),optional :: tag1,tag2,tag3,tag4,tag5
    if (present(tag1)) then
      if (0.ne.sep3d_close_tag(tag1,sep3df%tag)) then
        write(0,*) "tag=",trim(tag1)
        call seperr("trouble closing tag")
      end if
    end if
    if (present(tag2)) then
      if (0.ne.sep3d_close_tag(tag2,sep3df%tag)) then
        write(0,*) "tag=",trim(tag2)
        call seperr("trouble closing tag")
      end if
    end if
    if (present(tag3)) then
      if (0.ne.sep3d_close_tag(tag3,sep3df%tag)) then
        write(0,*) "tag=",trim(tag3)
        call seperr("trouble closing tag")
      end if
    end if
    if (present(tag4)) then
      if (0.ne.sep3d_close_tag(tag4,sep3df%tag)) then
        write(0,*) "tag=",trim(tag4)
        call seperr("trouble closing tag")
      end if
    end if
    if (present(tag5)) then
      if (0.ne.sep3d_close_tag(tag5,sep3df%tag)) then
        write(0,*) "tag=",trim(tag5)
        call seperr("trouble closing tag")
      end if
    end if
    if (0.ne.sep3d_delete(trim(sep3df%tag))) then
      call seperr("trouble deleting tag")
    end if
    call sep3d_clean(sep3df)
  end subroutine 
  subroutine sep3d_clean(sep3df)
    type(sep3d) :: sep3df
!        call sep3d_clear_headers(sep3df%tag)
!        call sep3d_clear_grid(sep3df%tag)
    if (sep3df%ndims.ne.0) then
      deallocate(sep3df%n,sep3df%o,sep3df%d,sep3df%label,sep3df%unit&
        &,sep3df%nwind,sep3df%fwind,sep3df%jwind)
    end if
    if (sep3df%nkeys.ne.0) then
      deallocate(sep3df%keyname,sep3df%keytype,sep3df%keyfmt)
    end if
    sep3df%usage="burned"
  end subroutine 
  subroutine struct_init(instruct,outstruct,usage,ctag)
    type(sep3d) :: instruct,outstruct
    character(len=*) :: usage
    character(len=128) :: temp,temp2
    character(len=*),optional :: ctag
    if (present(ctag)) then
      temp=ctag
    else
      write (temp,"(a,i1)") "sepf",counter
    end if 
    counter=counter+1
! call mkrandom_string(temp,temp2)
    if (0.ne.sep3d_struct_init(instruct%tag,temp,usage)) then
      call seperr("trouble initializing tag")
    end if
    call sep3d_grab_sep3d(temp,outstruct)
  end subroutine 
  subroutine tag_init(intag,outstruct,usage)
    type(sep3d) :: outstruct
    character(len=*) ::  intag,usage
    character(len=7) :: temp
    logical :: myloc
    integer :: ithread,nsz,i,nthread
    real, allocatable :: buf(:)
!write (temp,"(a,i1)") "sepf",counter;counter+=1
    ithread= sep_thread_num()
    nthread= sep_num_thread()
    if (0.ne. sep3d_tag_init_thread(intag,intag,usage,0)) then
      call seperr("trouble initializing tag")
    end if
    call sep3d_grab_sep3d(intag,outstruct)
  end subroutine 
  subroutine par_init(struct,usage,data_format,file_format,n,o,d,label&
    &,unit,ntraces,keyname,keytype,keyfmt,nh,ctag)
    type(sep3d) :: struct
    character(len=*) :: usage,data_format,file_format
    character(len=128) :: temp
    integer,optional :: n(:),ntraces,nh
    character(len=*), optional :: label(:),unit(:),keyname(:),keytype(:&
      &),keyfmt(:)
    real, optional :: d(:),o(:)
    integer :: ndim
    character(len=*),optional :: ctag
    logical :: do
    do=.true.
!if(len_trim(tag) > 16)
!  if(tag(1:17)=="/scr1/bob/stemp.H") do=.false.
    if (do) then
      if (present(ctag)) then
        temp=ctag
      else
        write (temp,"(a,i1)") "sepf",counter
        counter=counter+1
      end if 
      if (0.ne.sep3d_par_init(temp,usage)) then
        call seperr("trouble initializing structure")
      end if
      call sep3d_grab_sep3d(temp,struct)
      struct%file_format=file_format
      struct%data_format=data_format
      ndim=0
      if (present(n)) then
        ndim=size(n)
      end if
      if (present(o)) then
        if (ndim.ne.0) then
          ndim=size(o)
        else if (size(o).ne.ndim) then
          call seperr("size of axis parts must be consistent")
        end if
      end if
      if (present(d)) then
        if (ndim.ne.0) then
          ndim=size(d)
        else if (size(d).ne.ndim) then
          call seperr("size of axis parts must be consistent")
        end if
      end if
      if (present(label)) then
        if (ndim.ne.0) then
          ndim=size(label)
        else if (size(label).ne.ndim) then
          call seperr("size of axis parts must be consistent")
        end if
      end if
      if (present(unit)) then
        if (ndim.ne.0) then
          ndim=size(unit)
        else if (size(unit).ne.ndim) then
          call seperr("size of axis parts must be consistent")
        end if
      end if
      if (ndim.ne.0) then
        call axis_allocate(struct,ndim)
        if (present(n)) then
          struct%n=n
        else
          struct%n=1
        end if 
        if (present(o)) then
          struct%o=o
        else
          struct%o=0.
        end if 
        if (present(d)) then
          struct%d=d
        else
          struct%d=1.
        end if 
        if (present(label)) then
          struct%label=label
        else
          struct%label=C_NULL_CHAR
        end if 
        if (present(unit)) then
          struct%unit=unit
        else
          struct%unit="Undefined"//C_NULL_CHAR
        end if 
        struct%nwind=n
        struct%fwind=0
        struct%jwind=1
      end if
      if (present(keyname)) then
        if (.not. present(keytype) .or. .not. present(keyfmt)) then
          call seperr("must supply keyname, keytype, and keyfmt")
        else if (size(keyname) .ne. size(keyfmt) .or. size(keyname&
          &).ne.size(keyfmt)) then
          call seperr("keyname,keytype, and keyfmt must be the same size"&
            &)
        end if
        call key_allocate(struct,size(keyname))
        struct%keyname=keyname
        struct%keytype=keytype
        struct%keyfmt=keyfmt
      end if
      if (present(ntraces)) then
        struct%ntraces=ntraces
      end if
      if (present(nh)) then
        if (0.ne. sep3d_set_nh(struct%tag,nh)) then
          call seperr("trouble setting nh ")
        end if
      end if
      call sep3d_set_sep3d(struct)
    end if
  end subroutine 
  subroutine key_allocate(struct,nkeys)
    type(sep3d) :: struct
    integer :: nkeys
    if (.not.valid_structure(struct)) then
      call seperr("key_allocate: Structure has not been initialied")
    end if
    if (struct%nkeys.ne.0) then
      deallocate(struct%keyname,struct%keytype,struct%keyfmt)
    end if
    allocate(struct%keyname(nkeys),struct%keyfmt(nkeys),struct%keytype&
      &(nkeys))
    struct%keytype="scalar_float"
    struct%keyfmt="xdr_float"
    struct%nkeys=nkeys
  end subroutine 
  subroutine axis_allocate(struct,ndim)
    type(sep3d) :: struct
    integer :: ndim
    if (.not.valid_structure(struct)) then
      call seperr("axis_allocate: Structure has not been initialied")
    end if
    if (struct%ndims.ne.0) then
      deallocate(struct%n,struct%o,struct%d,struct%label,struct%unit&
        &,struct%nwind,struct%fwind,struct%jwind)
    end if
    allocate(struct%n(ndim),struct%o(ndim),struct%d(ndim),struct%label&
      &(ndim),struct%unit(ndim),struct%nwind(ndim),struct%fwind(ndim&
      &),struct%jwind(ndim))
    struct%ndims=ndim
    struct%n=1
    struct%o=0.
    struct%d=1.
    struct%label=" "
    struct%unit=" "
    struct%label(1)=C_NULL_CHAR
    struct%unit(1)=C_NULL_CHAR
  end subroutine 
  subroutine putkey_i_c(struct,keyname,header)
    type(sep3d) :: struct
    character(len=*) :: keyname
    integer,dimension(:) :: header
    integer :: nh
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_nh(struct%tag,nh)) then
      if (0.ne. sep3d_set_nh(struct%tag,size(header))) then
        call seperr("trouble setting number of headers in structure\n"&
          &)
      end if
    else if (nh> size(header)) then
      write(0,*) "size wrong ",nh,size(header)
      call seperr("header array not the same size as stored nh")
    end if
    if (0.ne.sep3d_set_header_vals_si(trim(struct%tag)//C_NULL_CHAR,trim(keyname)//C_NULL_CHAR,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine putkey_f_c(struct,keyname,header)
    type(sep3d) :: struct
    character(len=*) :: keyname
    real,dimension(:) :: header
    integer :: nh
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_nh(struct%tag,nh)) then
      if (0.ne. sep3d_set_nh(struct%tag,size(header))) then
        call seperr("trouble setting number of headers in structure\n"&
          &)
      end if
    else if (nh>size(header)) then
      write(0,*) "check",nh,size(header)
      call seperr("header array not the same size as stored nh")
    end if
    if (0.ne.sep3d_set_header_vals_sf(trim(struct%tag)//C_NULL_CHAR,trim(keyname)//C_NULL_CHAR,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine putkey_i_i(struct,keyindex,header)
    type(sep3d) :: struct
    integer :: keyindex,nh
    integer,dimension(:) :: header
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_nh(struct%tag,nh)) then
      if (0.ne. sep3d_set_nh(struct%tag,size(header))) then
        call seperr("trouble setting number of headers in structure\n"&
          &)
      end if
    else if (nh>size(header)) then
      call seperr("header array not the same size as stored nh")
    end if
    if (0.ne.sep3d_set_header_vals_ii(trim(struct%tag)//C_NULL_CHAR,keyindex,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine putkey_f_i(struct,keyindex,header)
    type(sep3d) :: struct
    integer :: keyindex,nh
    real,dimension(:) :: header
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_nh(struct%tag,nh)) then
      if (0.ne. sep3d_set_nh(struct%tag,size(header))) then
        call seperr("trouble setting number of headers in structure\n"&
          &)
      end if
    else if (nh>size(header)) then
      call seperr("header array not the same size as stored nh")
    end if
    if (0.ne.sep3d_set_header_vals_if(trim(struct%tag)//C_NULL_CHAR,keyindex,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine getkey_i_c(struct,keyname,header)
    type(sep3d) :: struct
    character(len=*) :: keyname
    integer,dimension(:) :: header
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_header_vals_si(trim(struct%tag)//C_NULL_CHAR,trim(keyname)//C_NULL_CHAR,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine getkey_f_c(struct,keyname,header)
    type(sep3d) :: struct
    character(len=*) :: keyname
    real,dimension(:) :: header
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_header_vals_sf(trim(struct%tag)//C_NULL_CHAR,trim(keyname)//C_NULL_CHAR,header)) then
      call seperr("trouble reading key")
    end if
  end subroutine 
  subroutine getkey_i_i(struct,keyindex,header)
    type(sep3d) :: struct
    integer :: keyindex
    integer,dimension(:) :: header
    real,dimension(:),allocatable :: headr
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%keytype(keyindex).eq."scalar_float") then
      allocate(headr(size(header)))
      if (0.ne.sep3d_grab_header_vals_if(trim(struct%tag)//C_NULL_CHAR,keyindex,headr))&
        & then
        call seperr("trouble reading key")
      end if
      header=nint(headr)
      deallocate(headr)
    else
      if (0.ne.sep3d_grab_header_vals_ii(trim(struct%tag)//C_NULL_CHAR,keyindex,header))&
        & then
          call seperr("trouble reading key")
      end if
    end if
  end subroutine 
  subroutine getkey_f_i(struct,keyindex,header)
    type(sep3d) :: struct
    integer :: keyindex
    real,dimension(:) :: header
    integer,dimension(:),allocatable :: headi
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%keytype(keyindex).eq."scalar_int") then
      allocate(headi(size(header)))
      if (0.ne.sep3d_grab_header_vals_ii(trim(struct%tag)//C_NULL_CHAR,keyindex,headi))&
        & then
        call seperr("trouble reading key")
      end if
      header=real(headi)
      deallocate(headi)
    else
      if (0.ne.sep3d_grab_header_vals_if(trim(struct%tag)//C_NULL_CHAR,keyindex,header))&
        & then
        call seperr("trouble reading key")
      end if
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_grab_headers - grab headers
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_grab_headers(tag,struct,nh,nwind,fwind,jwind)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag   -  char*
!!$
!!$      tag to read from
!!$
!!$=item struct-  sep3d
!!$
!!$      struct to read into
!!$
!!$=item nwind -  int*
!!$
!!$      number of elements
!!$
!!$=item fwind -  int*
!!$
!!$      first element
!!$
!!$=item jwind -  int*
!!$
!!$      sampling of elements
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item nh    -  int
!!$
!!$      number of headers read
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Read in headers to structure
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_read_data>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
!!$
  subroutine sep3d_grab_headers(tag,struct,nh,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    integer, optional, dimension(:) :: nwind,fwind,jwind
    integer, allocatable, dimension(:) :: ntemp,ftemp,jtemp,ng
    real, allocatable,dimension(:) :: og,dg
    character(len=128), allocatable,dimension(:) :: label,unit
    integer :: ndim,i,nh,nthread
    nthread= sep_num_thread()
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (0.ne.sep3d_grab_ndims(struct%tag,ndim)) then
      call seperr("trouble getting grid axis")
    end if
    ndim=ndim-1
    if(ndim>0) then
      allocate(ntemp(ndim),ftemp(ndim),jtemp(ndim))
      ntemp=1
      ftemp=0
      jtemp=1
      if (present(nwind)) then
!          if(size(nwind)<ndim) call seperr("nwind not the right size")
        ntemp(1:size(nwind))=nwind(1:size(nwind))
      else
        ntemp=-1
      end if 
      if (present(fwind)) then
!        if(size(fwind)<ndim .and. any(fwind(ndim:1)!=0) ) call seperr("fwind not the right size")
        ftemp(1:size(fwind))=fwind(1:size(fwind))
      else
        ftemp=-1
      end if 
      if (present(jwind)) then
        jtemp(1:size(jwind))=jwind(1:size(jwind))
      else
        jtemp=-1
      end if 
      call fix_vector(struct%n(2:),ftemp,jtemp,ntemp)
      struct%nwind(2:)=ntemp
      struct%fwind(2:)=ftemp
      struct%jwind(2:)=jtemp
    else
      allocate(ntemp(1),ftemp(1),jtemp(1))
      ntemp=1;ftemp=0;jtemp=1
    end if 
      if (0.ne.sep3d_read_headers(tag,struct%tag,ntemp,ftemp,jtemp,nh))&
        & then
        call seperr("trouble reading headers")
      end if
      deallocate(ftemp,jtemp,ntemp)
    end subroutine 
  logical function  get_data_3f(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    real,dimension(:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh
    integer :: nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "trouble fixing axis 1"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "FLOAT") then
      call seperr("data type is not float")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_3f=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_3c(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    complex,dimension(:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not.  fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "COMPLEX") then
      call seperr("data type is not complex")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_3c=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_3i(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    integer,dimension(:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    get_data_3i=.false.
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "INTEGER") then
      call seperr("data type is integer")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_3i=.true.
  end function 
  logical function  get_data_4f(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    real,dimension(:,:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh
    integer :: nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "trouble fixing axis 1"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "FLOAT") then
      call seperr("data type is not float")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_4f=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_4c(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    complex,dimension(:,:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not.  fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "COMPLEX") then
      call seperr("data type is not complex")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_4c=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_4i(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    integer,dimension(:,:,:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    get_data_4i=.false.
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "INTEGER") then
      call seperr("data type is integer")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_4i=.true.
  end function 
  logical function  get_data_2f(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    real,dimension(:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh
    integer :: nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "trouble fixing axis 1"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "FLOAT") then
      call seperr("data type is not float")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_2f=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_2c(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    complex,dimension(:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not.  fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "COMPLEX") then
      call seperr("data type is not complex")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_2c=.true.
!deallocate(trnum)
  end function 
  logical function  get_data_2i(tag,struct,data,nwind,fwind,jwind)
    character(len=*) :: tag
    type(sep3d) :: struct
    integer,dimension(:,:) :: data
    integer ::drn,ierr
    integer, dimension(:), allocatable ::trnum
    integer, optional :: fwind,jwind,nwind
    integer :: ft,nt,jt,nh,nthread
    nthread= sep_num_thread()
    get_data_2i=.false.
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (.not. fix1(struct%n(1),ft,jt,nt)) then
      write(0,*) "axis 1 problem"
      return
    end if
    struct%nwind(1)=nt
    struct%fwind(1)=ft
    struct%jwind(1)=jt
    if (.not. valid_structure(struct)) then
      call seperr("structure has not been initialized")
    end if
    if (struct%data_format .ne. "INTEGER") then
      call seperr("data type is integer")
    end if
    if (size(data) < product(struct%nwind) .and. struct%file_format(1:1&
      &).ne."G") then
      write(0,*) sep_thread_num(),"=thread Window size greater than&
        & buffer size",product(struct%nwind),size(data)
      call seperr("read request beyond buffer size")
    end if
    if (0.ne.sep3d_reed_data(tag,struct%tag,nt,ft,jt,data)) then
      call seperr("trouble reading data")
    end if
    get_data_2i=.true.
  end function 
! fix_window
! Usage
! fix_window(n_data,fwind,jwind,nwind)
! Input Parameter# n_data - integer(:) array describing the extent of the dataset
! fwind - integer(:) array of minimum values# jwind - integer(:) array of skip values
! nwind - integer(:) array of number of values tor read#
! Description
!
! Routine that checks if window params are valid.  By passing -1 it
! will also do the right thing in calculating reasonable values for# these parameters
!
! Category
! Base library - utility
!
  logical function fix1(n_data,f,j,n)
    integer,intent(in) :: n_data
    integer,intent(inout) ::  f,j,n
    fix1=.false.
    if (n_data <= f) then
      write(0,*) n_data,"=ndata f=",f
      write(0,*) "utils:invalid f parameter  "
      return
    end if
    if (n_data < n) then
      write(0,*) n_data,"=ndata n=",n
      write(0,*) "utils:invalid n parameter"
      return
    end if
    if (j .eq. -1) then
      j=1
    end if
    if (n .eq. -1) then
      if (f .eq. -1) then
        f=0
      end if
      n= (n_data-1-f)/j +1
!else if (f==-1) f=(n-1)/j +1
    else if (f.eq.-1) then
      f=0
    end if
    if (n_data < (1 + f + j * (n -1))) then
      write(0,*) "utils:invalid window parameters"
      return
    end if
    fix1=.true.
  end function 
  subroutine fix_vector(n_data,f,j,n)
    integer,dimension(:),intent(in)  ::  n_data
    integer,dimension(:),intent(inout)  ::  f,j,n
    integer  :: i1
    if (size(n_data) .ne.size(f) .or.  size(f) .ne. size(j) .or.  size&
      &(j) .ne. size(n)) then
      call seperr("utils:incompatible fix window vector lengths")
    end if
    do i1=1,size(n)
      if (.not. fix1(n_data(i1),f(i1),j(i1),n(i1))) then
        write(0,*) "trouble with axis ",i1," parameters"
        call seperr("")
      end if
    end do
  end subroutine 
  

  logical function put_data_4i(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    integer, dimension(:,:,:,:),target :: data
    integer, pointer      :: tempp(:,:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data)/size(data,1)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if

    if (0.ne.rite_i_4(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_4i=.true.
  end function 
  logical function put_data_4f(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    real, dimension(:,:,:,:),target:: data
    real, dimension(1) :: tempd
    real, pointer      :: tempp(:,:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer, target   :: gt(1)
    integer  :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data)/size(data,1)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_4f=.true.
  end function 
  logical function put_data_4c(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    complex, dimension(:,:,:,:),target :: data
    complex, dimension(1) :: tempd
    complex, pointer      :: tempp(:,:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer           :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    gp=size(data)/size(data,1)
    tempp=>data
    esize=8
    wd=1
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_4c=.true.
  end function 
  logical function put_data_3i(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    integer, dimension(:,:,:),target :: data
    integer, pointer      :: tempp(:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data)/size(data,1)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_3i=.true.
  end function 
  logical function put_data_3f(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    real, dimension(:,:,:),target:: data
    real, dimension(1) :: tempd
    real, pointer      :: tempp(:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer, target   :: gt(1)
    integer  :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data)/size(data,1)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_3f=.true.
  end function 
  logical function put_data_3c(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    complex, dimension(:,:,:),target :: data
    complex, dimension(1) :: tempd
    complex, pointer      :: tempp(:,:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer           :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    gp=size(data)/size(data,1)
    tempp=>data
    esize=8
    wd=1
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_3c=.true.
  end function 
  logical function put_data_2i(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    integer, dimension(:,:),target :: data
    integer, pointer      :: tempp(:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data,2)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_2i=.true.
  end function 
  logical function put_data_2f(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    real, dimension(:,:),target:: data
    real, dimension(1) :: tempd
    real, pointer      :: tempp(:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer, target   :: gt(1)
    integer  :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    esize=4
    tempp=>data
    wd=1
    gp=size(data,2)
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)
    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind(1:n)
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind(1:n)
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind(1:n)
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_2f=.true.
  end function 
  logical function put_data_2c(tag,struct,data,nwind,fwind,jwind&
    &,write_headers,write_grid)
    complex, dimension(:,:),target :: data
    complex, dimension(1) :: tempd
    complex, pointer      :: tempp(:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer           :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    gp=size(data,2)
    tempp=>data
    esize=8
    wd=1
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    put_data_2c=.true.
  end function 
  logical function put_data_n(tag,struct,nwind,fwind,jwind,write_headers&
    &,write_grid)
    complex, dimension(1,1),target :: tempd
    complex, pointer      :: tempp(:,:)
    character(len=*) :: tag
    type(sep3d)      :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer,allocatable,dimension(:) :: nt(:),ft(:),jt(:)
    logical, optional :: write_headers,write_grid
    integer           :: n
    integer,target    :: gt(1)
    integer, pointer  :: gp
    integer           :: wd,wg,wh,ntm,nh,ierr,esize,ng
    allocate(gp)
    gp=0
    esize=8
    tempp=>tempd
    wd=0
    if (.not. valid_structure(struct)) then
      call seperr("not a valid sep structure\n")
    end if
!if(struct%usage=="INPUT") call seperr("can not write out input struct \n");
!    n=struct%ndims
    n=sep3d_ndims(struct)

    allocate(nt(n),ft(n),jt(n))
    if (present(nwind)) then
      nt=nwind
    else
      nt=-1
    end if 
    if (present(fwind)) then
      ft=fwind
    else
      ft=-1
    end if 
    if (present(jwind)) then
      jt=jwind
    else
      jt=-1
    end if 
    if (present(write_grid)) then
      if (struct%file_format.ne."GRID" .and. write_grid) then
        call seperr("can not write grid when file_format not grid \n")
      end if
      if (struct%file_format.eq."GRID" .and. .not. write_grid) then
        call seperr("must write grid when file_format grid \n")
      end if
      if (write_grid) then
        wg=1
      else
        wg=0
      end if
    else
      if (struct%file_format.ne."GRID") then
        wg=0
      else
        wg=1
      end if
    end if 
    ierr=sep3d_grab_nh(struct%tag,nh)
    if (present(write_headers)) then
      if (struct%file_format.eq."REGULAR" .and. write_headers) then
        call seperr("file_format set to regular can not write out headers&
          & \n")
      end if
      if (struct%file_format.ne."REGULAR" .and. .not. write_headers)&
        & then
        call seperr("file_format set not regular must write out headers&
          & \n")
      end if
      if (write_headers) then
        wh=1
      else
        wh=0
      end if
    else
      if (struct%file_format.ne."REGULAR") then
        wh=1
      else
        wh=0
      end if
    end if 
    call fix_vector(struct%n(1:sep3d_ndims(struct)),ft,jt,nt)
    if (wd.eq.1 .and. struct%file_format(1:1).ne."G") then
      if (size(tempp)<product(nt)) then
        write(0,*) sep_thread_num(),"=thread window size=",product(nt)&
          &," buffer size=",size(tempp)
        call seperr("window size > buffer size")
      end if
    end if
    if (0.ne.sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))&
      & then
      call seperr("trouble writing out the data \n")
    end if
    deallocate(nt,ft,jt)
    deallocate(gp)
    put_data_n=.true.
  end function 
!/*<
!init_sepf90
!
!USAGE
!call init_sep3d()
!
!
!DESCRIPTION
!Initialize sep3d
!
!>*/
  subroutine init_sepf90()
    counter=0
  end subroutine 
!!$
!!$=head1 NAME
!!$
!!$print_sep3d - print info about the sep3d struct
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call print_sep3d(struct,level)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure to get info about
!!$
!!$=item level   -  int
!!$
!!$      level of verbosity
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Get info about strucut
!!$
!!$head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
  subroutine print_sep3d(struct,level)
    type(sep3d) :: struct 
    integer, optional :: level
    integer :: lv,i,ierr
    if (present(level)) then
      lv=level
    else
      lv=0
    end if 
    write(0,*) "----------FORTRAN 90 STATUS----------"
    write(0,*) "tag=",struct%tag(1:15)
    write(0,*) "data_format=",struct%data_format(1:15)
    write(0,*) "file_format=",struct%file_format(1:15)
    write(0,*) "usage=",struct%usage(1:len_trim(struct%usage))
    if (struct%ndims>0) then
      write(0,*)  "axis      n       o      d     label unit"
      do i=1,struct%ndims
        write(0,*) i,struct%n(i),struct%o(i),struct%d(i),struct%label(i&
          &)(1:15)," ",struct%unit(i)(1:15)
      end do
    else
      write(0,*) "NO AXES ALLOCATED"
    end if 
    if (struct%nkeys>0) then
      write(0,*)  "key number   name    type format"
      do i=1,struct%nkeys
        write(0,*) i,struct%keyname(i)(1:15)," ",struct%keytype(i)(1:15&
          &)," ",struct%keyfmt(i)(1:15)
      end do
    else
      write(0,*) "NO KEYS ALLOCATED"
    end if 
    write(0,*) "ntraces",struct%ntraces
    write(0,*) "----------------- C  STATUS----------"
    ierr=sep3d_print_info(struct%tag)
  end subroutine 
!!$
!!$=head1 NAME
!!$
!!$sep3d_write_description - write out format file info
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<sep3d_write_description(tag,struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag      -  char*
!!$
!!$      tag to write out
!!$
!!$=item struct   -  sep3d
!!$
!!$      structure to write out
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!!$
!!$Write out description to tag
!!$
!!$=head1 SEE ALSO
!!$
!!!$B<sep3d_write>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!!$=cut
!!$
  subroutine sep3d_write_description(tag,structure)
    type(sep3d) :: structure
    character(len=*) :: tag
    if (.not. valid_structure(structure)) then
      call seperr("you are attempting to write out an invalid strucutre"&
        &)
    end if
    call sep3d_set_sep3d(structure)
    if (0.ne.sep3d_set_ntraces(structure%tag,structure%ntraces)) then
      call seperr("trouble writing out structure to disk1")
    end if
    if (0.ne.sep3d_rite_format(tag,structure%tag)) then
      call seperr("trouble writing out structure to disk2")
    end if
        if (0.ne.sep3d_rite_ntraces(tag,structure%tag)) then
      call seperr("trouble writing out structure to disk3")
    end if
  end subroutine 
  subroutine sep3d_set_write_status(structure,data,header)
    type(sep3d) :: structure
    logical :: header,data
    integer :: h,d
    if (.not. valid_structure(structure)) then
      call seperr("you are attempting to write out an invalid strucutre"&
        &)
    end if
    call sep3d_set_sep3d(structure)
    if (header) then
      h=1
    else
      h=0
    end if 
    if (data) then
      d=1
    else
      d=0
    end if 
    if (0.ne.sep3d_rite_status(structure%tag,d,h)) then
      call seperr("trouble writing out structure to disk ")
    end if
  end subroutine 
!!$
!!$=head1 NAME
!!$
!!$sep3d_add_drn - add data record number
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_add_drn(struct,drn)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure to add drn too
!!$
!!$=item drn     -  int*
!!$
!!$      drn to add
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Add data record number to sep3d structure
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
!!$
!!$
  subroutine sep3d_add_drn(structure,drn)
    type(sep3d) :: structure
    integer :: drn(:)
    if (.not. valid_structure(structure)) then
      call seperr("you are attempting to write out an invalid strucutre"&
        &)
    end if
    if (0.ne.sep3d_set_drn (structure%tag,drn)) then
      call seperr("trouble setting drn")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_key_index - try to find key in structure
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logical sep3d_key_index(struct,keyname,keyindex)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure to find key in
!!$
!!$=item keyname - char*
!!$
!!$       name of the key to find
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item keyindex- integer
!!$
!!$      index of the key
!!$
!!$=back
!!$
!!$=head1 RETURN VALUES
!!$
!!$=over 4
!!$
!!$=item .true.  -
!!$
!!$      found key
!!$
!!$=item .false. -
!!$
!!$       didn't find key
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Attempte to find the given key
!!$
!!$=head1 SEE ALSO
!!$
!!!$L<sep3d_axis_index>
!!$
!!$=head1 LIBRARY
!!$
!!!$B<supersetf90>
!!$
!!$=cut
!>
  logical  function sep3d_key_index(struct,keyname,keyindex)
    type(sep3d) :: struct
    integer :: keyindex
    character(len=*) ::  keyname
    integer :: ierr
    if (.not. valid_structure(struct)) then
      call seperr("you are attempting to write out an invalid strucutre"&
        &)
    end if
    call sep3d_set_sep3d(struct)
    ierr=sep3d_grab_key_index(struct%tag,keyname,keyindex)
    if (ierr.eq.0) then
      sep3d_key_index=.true.
    else
      sep3d_key_index=.false.
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_axis_index - try to axis key in structure
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logical sep3d_axis_index(struct,axisname,axisindex)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure to find key in
!!$
!!$=item axisname - char*
!!$
!!$      name of the axis to find
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item axisindex- integer
!!$
!!$      index of the axis
!!$
!!$=back
!!$
!!$=head1 RETURN VALUES
!!$
!!$=over 4
!!$
!!$=item .true.  -
!!$
!!$      found axis
!!$
!!$=item .false. -
!!$
!!$      didn't find axis
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Attempte to find the given axis
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_key_index>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical  function sep3d_axis_index(struct,axisname,axisindex)
    type(sep3d) :: struct
    integer :: axisindex
    character(len=*) ::  axisname
    integer :: ierr
    if (.not. valid_structure(struct)) then
      call seperr("you are attempting to write out an invalid strucutre"&
        &)
    end if
    call sep3d_set_sep3d(struct)
    ierr=sep3d_grab_axis_index(struct%tag,axisname,axisindex)
    if (ierr.eq.0) then
      sep3d_axis_index=.true.
    else
      sep3d_axis_index=.false.
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_set_number_coords - set the number of coords to store
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_set_number_coords(struct,nh)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=item nh      -  integer
!!!$
!!$      number of coords
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!!$Set number of coords to hold in structure
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_store_grid_values>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!!$=cut
  subroutine sep3d_set_number_coords(struct,num)
    type(sep3d) :: struct
    integer :: num
    if (.not. valid_structure(struct)) then
      call seperr("invalid structure ")
    end if
    if (0.ne.sep3d_alloc_coord(struct%tag,num)) then
      call seperr("trouble setting number of coords ")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_set_number_headers - set the number of headers to store
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!!$C<call sep3d_set_number_headers(struct,nh)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!!$      structure
!!$
!!$=item nh      -  integer
!!$
!!$      number of headers
!!!$
!!$=back
!!$
!!!$=head1 DESCRIPTION
!!$
!!$Set number of headers to hold in structure
!!$
!!$=head1 SEE ALSO
!!!$
!!$L<sep3d_store_grid_values>
!!$
!!$=head1 LIBRARY
!!$
!!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_set_number_headers(struct,num)
    type(sep3d) :: struct
    integer :: num
    if (.not. valid_structure(struct)) then
      call seperr("invalid structure ")
    end if
    if (0.ne.sep3d_set_nh(struct%tag,num)) then
      call seperr("trouble setting number of headers ")
    end if
    call sep3d_set_number_coords(struct,num)
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_grab_grid_values - store  a grid
!!$
!!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_grab_grid_values(struct,grid)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item grid - C<int(:)>
!!$
!!$      grid
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab grid values in C strucutre
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_set_number_headers>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_grab_grid_values(struct,vals)
    type(sep3d) :: struct
    integer     :: vals(:)
    if (.not. valid_structure(struct)) then
      call seperr("invalid structure ")
    end if
    if (0.ne.sep3d_grab_grid_block(struct%tag,vals)) then
      call seperr("trouble storing grid values")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_store_grid_values - store  a grid
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_store_grid_values(struct,grid)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=item grid - C<int(:)>
!!$
!!$      grid
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Store grid values in C strucutre
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep3d_set_number_headers>
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_store_grid_values(struct,vals)
    type(sep3d) :: struct
    integer     :: vals(:)
    if (.not. valid_structure(struct)) then
      call seperr("invalid structure ")
    end if
    if (0.ne.sep3d_set_ng(struct%tag,size(vals))) then
      call seperr("trouble setting number of grid values")
    end if
    if (0.ne.sep3d_set_grid_vals(struct%tag,vals)) then
      call seperr("trouble storing grid values")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_with_drn - sets
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<ndim= sep3d_with_drn(without,with)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item without  -  sep3d
!!$
!!$      structure  without drn
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item with  -  sep3d
!!$
!!$      structure  without drn included
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Returns the sep3d structure with the drn record number included in the
!!$  correct key position
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_with_drn(real_struct,drn_struct)
    type(sep3d) :: real_struct,drn_struct
    integer :: ierr,i,index,ito
    call init_sep3d(real_struct,drn_struct,real_struct%usage)
    ierr=sep_get_key_index(real_struct%tag,"data_record_number",index)
    if (ierr.eq.0) then
!we have a data_record_number
      call key_allocate(drn_struct,real_struct%nkeys+1)
      do i=1,real_struct%nkeys 
        if (i.eq. index) then
!drn key location
          drn_struct%keyname(i)="data_record_number"
          drn_struct%keytype(i)="scalar_int"
          drn_struct%keyfmt(i)="xdr_int"
        else
          if (i<index) then
            ito=i
          else
            ito=i+1
          end if 
          drn_struct%keyname(ito)=real_struct%keyname(i)
          drn_struct%keytype(ito)=real_struct%keytype(i)
          drn_struct%keyfmt(ito)=real_struct%keyfmt(i)
        end if
      end do
    end if
    call sep3d_set_sep3d(drn_struct)
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_ndims - number of dimensions in sep3d datasets
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<ndim= sep3d_ndims(struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$The last dimension with axis length greater than 1
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  integer function sep3d_ndims(struct)
    type(sep3d) :: struct
    integer :: i
    if (.not. valid_structure(struct)) then
      sep3d_ndims=-1
    else
      sep3d_ndims=struct%ndims
      do while (sep3d_ndims.ne.1 .and.  struct%n(sep3d_ndims).eq.1)
        sep3d_ndims= sep3d_ndims-1
      end do
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_ndims - number of dimensions in sep3d datasets
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_inorder(struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Establish that the data and the headers are synched
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_inorder(struct)
    type(sep3d) :: struct
    if (0.ne.sep3d_set_inorder(struct%tag)) then
      write(0,*) "trouble setting in order ",trim(struct%tag)
      call seperr("")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_rite_num_traces - write number of traces in a dataset
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<sep3d_rite_num_traces(tag,struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag  -   sepfile
!!$
!!$      tag to write to
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Write  the number of traces in a dataset to a tag
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_rite_num_traces(tag,struct)
    type(sep3d) :: struct
    character(len=*) :: tag
    if (0.ne.sep3d_rite_ntraces(tag,struct%tag)) then
      write(0,*) "trouble writing number of traces for ",trim(tag)
      call seperr("")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_reshape -  reshape a dataset
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call  sep3d_reshape(struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=item n  -  int*
!!$
!!$      axis remapping
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Change the dimensions of a dataset
!!$
!!$=head1 EXAMPLES
!!$
!!$ Given data%n=(/10,20,10/)  n=(/1,2,3,3/) --> data%n(/10,20,10,1/)
!!$
!!$ Given data%n=(/10,20,10/)  n=(/1,3/) --> data%n(/10,200/)
!!$
!!$ Given data%n=(/10,20,10/)  n=(/1,2,2,3/) --> data%n(/10,1,20,10/)
!!$
!!$
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine  sep3d_reshape(struct,n)
    type(sep3d) :: struct
    integer :: n(:)
    call sep3d_set_sep3d(struct)
    if (0.ne.sep3d_change_dims(struct%tag,size(n),n)) then
      write(0,*) "trouble changing dimensions for tag ",trim(struct%tag&
        &)
      call seperr("")
    end if
    call sep3d_grab_sep3d(struct%tag,struct)
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_header_copy - Copy the headers from one dataset to another
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_header_copy(to,from)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item to  -  sep3d
!!$
!!$      structure
!!$
!!$=item from  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Copy  the headers from one dataset to another
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_coord_copy(struct1,struct2)
    type(sep3d) :: struct1,struct2
    if (0.ne.sep3d_copy_coords(struct1%tag,struct2%tag)) then
      write(0,*) "trouble copying coords from ",trim(struct1%tag), " to&
        & ",trim(struct2%tag)
      call seperr("")
    end if
  end subroutine 
  subroutine sep3d_header_copy(struct1,struct2)
    type(sep3d) :: struct1,struct2
    if (0.ne.sep3d_copy_headers(struct1%tag,struct2%tag)) then
      write(0,*) "trouble copying header from ",trim(struct1%tag), " to&
        & ",trim(struct2%tag)
      call seperr("")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_grid_copy - Copy the grid
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_grid_copy(in,out)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item in  -  sep3d
!!$
!!$      structure
!!$
!!$=item out  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Copy the grid from in to out
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_grid_copy(struct1,struct2)
    type(sep3d) :: struct1,struct2
    if (0.ne.sep3d_copy_grid(struct1%tag,struct2%tag)) then
      write(0,*) "trouble copying grid from ",trim(struct1%tag), " to "&
        &,trim(struct2%tag)
      call seperr("")
    end if
  end subroutine 
!!$=head1 NAME
!!$
!!$sep3d_ge_space - see if space contains another space
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_ge_space(struct1,struct2)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct1  -  sep3d
!!$
!!$      structure
!!$
!!$=item struct2  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Check to see if struct2 fits within the space occupied by struct1
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_ge_space(sep3d1,sep3d2)
    type(sep3d) :: sep3d1,sep3d2
    sep3d_ge_space=.true.
    if (0.ne.sep3dge_space(sep3d1%tag,sep3d2%tag)) then
      sep3d_ge_space=.false.
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_conform - see if two spaces are the same size
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_conform(struct1,struct2)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct1  -  sep3d
!!$
!!$      structure
!!$
!!$=item struct2  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Check to see if two spaces have the same n,o, and d
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_conform(sep3d1,sep3d2)
    type(sep3d) :: sep3d1,sep3d2
    sep3d_conform=.true.
    if (0.ne.sep3dconform(sep3d1%tag,sep3d2%tag)) then
      sep3d_conform=.false.
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_copy - a copy a structure and its contents
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_copy(struct1,struct2)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct1  -  sep3d
!!$
!!$      structure
!!$
!!$=item struct2  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Copy from struct1 to struct2
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_copy(struct_in,struct_out)
    type(sep3d) ::  struct_in,struct_out
    integer :: ierr
    sep3d_copy=.false.
    ierr=sep3d_copy_struct(struct_in%tag,struct_out%tag)
    call sep3d_grab_sep3d(struct_out%tag,struct_out)
    if (ierr.eq.0) then
      sep3d_copy=.true.
    end if
  end function 
!!$=head1 NAME
!!$
!!$sep3d_grab_coords -  Grab the coordinates
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_grab_coords(struct,coords)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=item coords  -  int (:,:)
!!$
!!$      coordinates
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab the coordinates [axis(2:)] of the dataset currently in memory
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_grab_coords(struct,coords)
    type(sep3d) :: struct
    integer :: coords(:,:)
    sep3d_grab_coords=.false.
    if (0.ne.sep3d_grab_coord_vals(struct%tag,coords)) then
      write(0,*) "trouble grabbing coord vals"
      return
    end if
    sep3d_grab_coords=.true.
  end function 
!!$=head1 NAME
!!$
!!$sep3d_update_ntraces - Update the traces of a dataset
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_update_ntraces(struct)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ All write statements count the number of traces written out. This routine will
!!$ transfer that count to the ntraces of the dataset
!!$
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_update_ntraces(struct)
    type(sep3d) :: struct
    sep3d_update_ntraces=.false.
    if (0.ne.sep3d_count_ntraces(struct%tag)) then
      write(0,*) "trouble grabbing coord vals"
      return
    end if
    call sep3d_grab_sep3d(struct%tag,struct)
    sep3d_update_ntraces=.true.
  end function 
!!$=head1 NAME
!!$
!!$sep3d_set_coords - Specify the coordinates of a dataset in memory
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<logic= sep3d_set_coords(struct,coords)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct  -  sep3d
!!$
!!$      structure
!!$
!!$=item coords  -  int(:,:)
!!$
!!$      coordinates
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Set the coordinates of the dataset in memory
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  logical function sep3d_set_coords(struct,coords)
    type(sep3d) :: struct
    integer :: coords(:,:)
    sep3d_set_coords=.false.
!call sep3d_set_number_headers(struct,size(coords,2))
    if (0.ne.sep3d_alloc_coord(struct%tag,size(coords,2))) then
      write(0,*) "trouble setting size of coord"
      return
    end if
    if (0.ne.sep3d_set_coord_vals(struct%tag,coords)) then
      write(0,*) "trouble setting coord vals"
      return
    end if
    sep3d_set_coords=.true.
  end function 
!!$=head1 NAME
!!$
!!$sep3d_set_window - specify a subsection of data for a tag
!!$
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep3d_set_window(struct_in,nwind,fwind,jwin)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item struct_in  -  sep3d
!!$
!!$      structure
!!$
!!$=item nwind  -  int (:)
!!$
!!$      number of elements in a window (optional)
!!$
!!$=item fwind  -  int (:)
!!$
!!$      first  elements in a window (optional)
!!$
!!$=item jwind  -  int (:)
!!$
!!$      sampling in a window (optional)
!!$
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$  Specify a subsection for a dataset
!!$
!!$=head1 SEE ALSO
!!$
!!$
!!$=head1 LIBRARY
!!$
!!$B<supersetf90>
!!$
!!$=cut
  subroutine sep3d_set_window(struct,nwind,fwind,jwind)
    type(sep3d) :: struct
    integer,optional :: nwind(:),fwind(:),jwind(:)
    integer :: ndim,i
    ndim=sep3d_ndims(struct)
    struct%nwind=-1
    struct%fwind=-1
    struct%jwind=-1
    if (present(nwind)) then
      struct%nwind=nwind
    end if
    if (present(fwind)) then
      struct%fwind=fwind
    end if
    if (present(jwind)) then
      struct%jwind=jwind
    end if
    do i=1,ndim 
      if (.not. fix1(struct%n(i),struct%fwind(i),struct%jwind(i)&
        &,struct%nwind(i))) then
        call seperr("illegal windowing parameters")
      end if
    end do 
    if (0.ne. sep3d_set_wind(struct%tag,struct%nwind,struct%fwind&
      &,struct%jwind)) then
      call seperr("trouble setting window parameters")
    end if
  end subroutine 
  subroutine sep3d_add_key(struct,keyname,keytype,keyfmt)
    type(sep3d) :: struct
    character(len=*) :: keyname,keytype,keyfmt
    character(len=128),pointer :: keynold(:),keytold(:),keyfold(:)
    integer ::i
    if (struct%nkeys >0) then
      allocate(keynold(struct%nkeys))
      keynold=struct%keyname
      allocate(keytold(struct%nkeys))
      keytold=struct%keytype
      allocate(keyfold(struct%nkeys))
      keyfold=struct%keyfmt
      struct%nkeys=struct%nkeys+1 
      call key_allocate(struct,struct%nkeys)
      struct%keyname(1:struct%nkeys-1)=keynold
      struct%keytype(1:struct%nkeys-1)=keytold
      struct%keyfmt(1:struct%nkeys-1)=keyfold
      deallocate(keynold,keytold,keyfold)
    else
      struct%nkeys=1
      call key_allocate(struct,1)
    end if 
    struct%keyname(struct%nkeys)=keyname
    struct%keytype(struct%nkeys)=keytype
    struct%keyfmt(struct%nkeys)=keyfmt
  end subroutine 
  logical function sep3d_file_exist(tag,local)
    integer :: ilocal,i
    character(len=*) :: tag
    logical, optional :: local
    sep3d_file_exist=.false.
    ilocal=0
    if (present(local)) then
      if (local) then
        ilocal=1
      end if
    end if
    i=sep3d_fileexist(tag,ilocal)
    if (i.eq.1) then
      sep3d_file_exist=.true.
    end if
  end function 
  
  integer function sep3d_par_init(tag1,tag2)
    character(len=*) ,intent(in) :: tag1,tag2
    integer lentag1, lentag2
    lentag1 = LEN_TRIM(tag1)
    lentag2 = LEN_TRIM(tag2)
    sep3d_par_init=sep3d_par_initf(tag1(1:lentag1)//C_NULL_CHAR,tag2(1:lentag2)//C_NULL_CHAR)
 end function
 
    integer function sep3d_struct_init(tag1,tag2,typ)
    character(len=*) ,intent(in) :: tag1,tag2,typ
    sep3d_struct_init=sep3d_struct_initf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR,&
      trim(typ)//C_NULL_CHAR)
 end function
    integer function sep3d_tag_init_thread(tag1,tag2,typ,ind)
    character(len=*) ,intent(in) :: tag1,tag2,typ
    integer,intent(in) :: ind
    sep3d_tag_init_thread=sep3d_tag_init_threadf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR,&
     trim(typ)//C_NULL_CHAR,ind)
 end function
   integer function sep3d_tag_init(tag1,tag2,typ)
    character(len=*) ,intent(in) :: tag1,tag2,typ
    sep3d_tag_init=sep3d_tag_initf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR)
 end function
     integer function sep3d_grab_nkeys(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) :: ntr
    sep3d_grab_nkeys=sep3d_grab_nkeysf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
      integer function sep3d_set_nkeys(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_set_nkeys=sep3d_set_nkeysf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
      integer function sep3d_grab_nh(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) :: ntr
    sep3d_grab_nh=sep3d_grab_nhf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
       integer function sep3d_grab_ncoord(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) :: ntr
    sep3d_grab_ncoord=sep3d_grab_ncoordf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
       integer function sep3d_set_nh(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_set_nh=sep3d_set_nhf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
        integer function sep3d_set_ndims(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_set_ndims=sep3d_set_ndimsf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
        integer function sep3d_set_ntraces(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_set_ntraces=sep3d_set_ntracesf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
         integer function sep3d_grab_drn(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) ::ntr(:)
    sep3d_grab_drn=sep3d_grab_drnf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
         integer function sep3d_set_drn(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) ::ntr(:)
    sep3d_set_drn=sep3d_set_drnf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
        integer function sep3d_fileexist(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_fileexist=sep3d_fileexistf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
      integer function sep3d_grab_ndims(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) :: ntr
    sep3d_grab_ndims=sep3d_grab_ndimsf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
    integer function sep3d_grab_ntraces(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(out) :: ntr
    sep3d_grab_ntraces=sep3d_grab_ntracesf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
     integer function sep3d_alloc_coord(tag1,ntr)
    character(len=*) ,intent(in) :: tag1
    integer,intent(in) :: ntr
    sep3d_alloc_coord=sep3d_alloc_coordf(trim(tag1)//C_NULL_CHAR,ntr)
 end function
     integer function sep3d_grab_data_type(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(out) :: typ
    sep3d_grab_data_type=sep3d_grab_data_typef(trim(tag)//C_NULL_CHAR,typ)
    call c2forstr(typ)
 end function
      integer function sep3d_set_data_type(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_set_data_type=sep3d_set_data_typef(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR)
 end function
       integer function sep3d_rite_status(tag,c1,c2)
    character(len=*) ,intent(in) :: tag
   integer(C_INT),intent(in):: c1,c2
    sep3d_rite_status=sep3d_rite_statusf(trim(tag)//C_NULL_CHAR,c1,c2)
 end function
       integer function sep3d_set_file_type(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_set_file_type=sep3d_set_file_typef(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR)
 end function
      integer function sep3d_grab_file_type(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(out) :: typ
    sep3d_grab_file_type=sep3d_grab_file_typef(trim(tag)//C_NULL_CHAR,typ)
call c2forstr(typ)
 end function
       integer function sep3d_count_NTRACES(tag)
    character(len=*) ,intent(in) :: tag
    sep3d_count_ntraces=sep3d_count_ntracesf(trim(tag)//C_NULL_CHAR)
 end function
       integer function sep3d_grab_usage(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(out) :: typ
    sep3d_grab_usage=sep3d_grab_usagef(trim(tag)//C_NULL_CHAR,typ)
    call c2forstr(typ)
 end function
        integer function sep3d_close_tag(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_close_tag=sep3d_close_tagf(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR)
 end function
         integer function sep3d_delete(tag)
    character(len=*) ,intent(in) :: tag
    sep3d_delete=sep3d_deletef(trim(tag)//C_NULL_CHAR)
 end function
        integer function sep3d_copy_struct(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_copy_struct=sep3d_copy_structf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
         integer function sep3dconform(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3dconform=sep3d_conformf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
         integer function sep3dge_space(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3dge_space=sep3d_ge_spacef(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
           integer function sep3d_copy_headers(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_copy_headers=sep3d_copy_headersf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
            integer function sep3d_copy_coords(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_copy_coords=sep3d_copy_coordsf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
          integer function sep3d_copy_grid(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_copy_grid=sep3d_copy_gridf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
           integer function sep3d_rite_format(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_rite_format=sep3d_rite_formatf(trim(tag)//C_NULL_CHAR,trim(tag)//C_NULL_CHAR)
 end function
        integer function sep3d_rite_ntraces(tag,typ)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
    sep3d_rite_ntraces=sep3d_rite_ntracesf(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR)
 end function
         integer function sep3d_grab_axis_index(tag,typ,ind)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
   integer,intent(out) :: ind
    sep3d_grab_axis_index=sep3d_grab_axis_indexf(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,ind)
 end function
 
          integer function sep3d_grab_key_index(tag,typ,ind)
    character(len=*) ,intent(in) :: tag
   character(len=*),intent(in) :: typ
   integer,intent(out) :: ind
    sep3d_grab_key_index=sep3d_grab_key_indexf(trim(tag)//C_NULL_CHAR,trim(typ)//C_NULL_CHAR,ind)
 end function
        integer function sep3d_grab_key(tag,index,nm,typ,fmt)
    character(len=*) ,intent(in) :: tag
    integer,intent(in) :: index
   character(len=*),intent(out) :: nm,typ,fmt
    sep3d_grab_key=sep3d_grab_keyf(trim(tag)//C_NULL_CHAR,index,nm,typ,fmt)
    call c2forstr(nm)
        call c2forstr(typ)
    call c2forstr(fmt)

 end function
         integer function sep3d_set_key(tag,index,nm,typ,fmt)
    character(len=*) ,intent(in) :: tag
    integer,intent(in) :: index
   character(len=*),intent(in) :: nm,typ,fmt
    sep3d_set_key=sep3d_set_keyf(trim(tag)//C_NULL_CHAR,index,trim(nm)//C_NULL_CHAR, &
      trim(typ)//C_NULL_CHAR,trim(fmt)//C_NULL_CHAR)
 end function
     integer function sep3d_grab_wind(tag1,n,f,j)
    character(len=*) ,intent(in) :: tag1
    integer,dimension(:) :: n,j,f
    sep3d_grab_wind=sep3d_grab_windf(trim(tag1)//C_NULL_CHAR,n,f,j)
 end function
      integer function sep3d_set_wind(tag1,n,f,j)
    character(len=*) ,intent(in) :: tag1
    integer,dimension(:) :: n,j,f
    sep3d_set_wind=sep3d_set_windf(trim(tag1)//C_NULL_CHAR,n,f,j)
 end function
       integer function sep3d_read_headers(tag1,tag2,n,f,j,v)
    character(len=*) ,intent(in) :: tag1,tag2
    integer,dimension(:),intent(in) :: n,j,f
    integer,intent(out):: v
    sep3d_read_headers=sep3d_read_headersf(trim(tag1)//C_NULL_CHAR,trim(tag2)//C_NULL_CHAR,n,f,j,v)
 end function
 
       integer function sep3d_set_coord_vals(tag,vals)
    character(len=*) ,intent(in) :: tag
    integer,dimension(:,:),intent(in) :: vals
    sep3d_set_coord_vals=sep3d_set_coord_valsf(trim(tag)//C_NULL_CHAR,vals)
 end function
        integer function sep3d_set_inorder(tag)
    character(len=*) ,intent(in) :: tag
    sep3d_set_inorder=sep3d_inorderf(trim(tag)//C_NULL_CHAR)
 end function
         integer function sep3d_print_info(tag)
    character(len=*) ,intent(in) :: tag
    sep3d_print_info=sep3d_print_infof(trim(tag)//C_NULL_CHAR)
 end function
        integer function sep3d_grab_coord_vals(tag,vals)
    character(len=*) ,intent(in) :: tag
    integer,dimension(:,:),intent(out) :: vals
    sep3d_grab_coord_vals=sep3d_grab_coord_valsf(trim(tag)//C_NULL_CHAR,vals)
 end function
      integer function sep3d_grab_axis(tag,ind,n,o,d,lab,unit)
    character(len=*) ,intent(in) :: tag
    integer,intent(in) :: ind
    integer,intent(out) :: n
    real,intent(out) :: o,d
    character(len=*),intent(out) :: lab,unit
    sep3d_grab_axis=sep3d_grab_axisf(trim(tag)//C_NULL_CHAR,ind,n,o,d,lab,unit)
    call c2forstr(lab)
    call c2forstr(unit)
 end function
 
       integer function sep3d_set_axis(tag,ind,n,o,d,lab,unit)
    character(len=*) ,intent(in) :: tag
    integer,intent(in) :: ind
    integer,intent(in) :: n
    real,intent(in) :: o,d
    character(len=*),intent(in) :: lab,unit
    sep3d_set_axis=sep3d_set_axisf(trim(tag)//C_NULL_CHAR,ind,n,o,d,trim(lab)//C_NULL_CHAR,&
      trim(unit)//C_NULL_CHAR)
 end function
 integer function sep3d_change_dims(tag,ndim,n)
   character(len=*) ,intent(in) :: tag
    integer,intent(in) :: ndim,n(:)
    sep3d_change_dims=sep3d_change_dimsf(trim(tag)//C_NULL_CHAR,ndim,n)
end function
 integer function sep3d_set_ng(tag,ndim)
   character(len=*) ,intent(in) :: tag
    integer,intent(in) :: ndim
    sep3d_set_ng=sep3d_set_ngf(trim(tag)//C_NULL_CHAR,ndim)
end function
 integer function sep3d_grab_ng(tag,ndim)
   character(len=*) ,intent(in) :: tag
    integer,intent(out) :: ndim
    sep3d_grab_ng=sep3d_grab_ngf(trim(tag)//C_NULL_CHAR,ndim)
end function
 integer function sep3d_grab_grid_block(tag,n)
   character(len=*) ,intent(in) :: tag
    integer,intent(out) :: n(:)
    sep3d_grab_grid_block=sep3d_grab_grid_blockf(trim(tag)//C_NULL_CHAR,n)
end function
 integer function sep3d_set_grid_vals(tag,n)
   character(len=*) ,intent(in) :: tag
    integer,intent(in) :: n(:)
    sep3d_set_grid_vals=sep3d_set_grid_valsf(trim(tag)//C_NULL_CHAR,n)
end function
integer function read_f_1(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 real,intent(out) :: v(:)
   read_f_1=sep3d_reedff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
end function
  integer function read_f_2(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 real,intent(out) :: v(:,:)
   read_f_2=sep3d_reedff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
   end function
  integer function read_f_3(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 real,intent(out) :: v(:,:,:)
   read_f_3=sep3d_reedff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_f_4(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 real,intent(out) :: v(:,:,:,:)
   read_f_4=sep3d_reedff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_f_5(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 real,intent(out) :: v(:,:,:,:,:)
   read_f_5=sep3d_reedff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function
      
      
      
integer function read_i_1(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  integer,intent(out) :: v(:)
   read_i_1=sep3d_reedif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
end function
  integer function read_i_2(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 integer,intent(out) :: v(:,:)
   read_i_2=sep3d_reedif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
   end function
  integer function read_i_3(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  integer,intent(out) :: v(:,:,:)
   read_i_3=sep3d_reedif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_i_4(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  integer,intent(out) :: v(:,:,:,:)
   read_i_4=sep3d_reedif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_i_5(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  integer,intent(out) :: v(:,:,:,:,:)
   read_i_5=sep3d_reedif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function
      
      integer function read_c_1(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  complex,intent(out) :: v(:)
   read_c_1=sep3d_reedcf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
end function
  integer function read_c_2(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
 complex,intent(out) :: v(:,:)
   read_c_2=sep3d_reedcf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
   end function
  integer function read_c_3(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  complex,intent(out) :: v(:,:,:)
   read_c_3=sep3d_reedcf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_c_4(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  complex,intent(out) :: v(:,:,:,:)
   read_c_4=sep3d_reedcf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function

  integer function read_c_5(t1,t2,i1,i2,i3,v)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3
  complex,intent(out) :: v(:,:,:,:,:)
   read_c_5=sep3d_reedcf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,i1,i2,i3,v)
      end function
      

  integer function rite_f_1(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 real,intent(in) :: v(:)
   rite_f_1=sep3d_riteff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_f_2(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 real,intent(in) :: v(:,:)
   rite_f_2=sep3d_riteff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_f_3(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 real,intent(in) :: v(:,:,:)
   rite_f_3=sep3d_riteff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_f_4(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 real,intent(in) :: v(:,:,:,:)
   rite_f_4=sep3d_riteff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
  integer function rite_f_5(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 real,intent(in) :: v(:,:,:,:,:)
   rite_f_5=sep3d_riteff(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function


  integer function rite_c_1(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 complex,intent(in) :: v(:)
   rite_c_1=sep3d_ritecf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_c_2(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 complex,intent(in) :: v(:,:)
   rite_c_2=sep3d_ritecf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_c_3(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 complex,intent(in) :: v(:,:,:)
   rite_c_3=sep3d_ritecf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_c_4(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 complex,intent(in) :: v(:,:,:,:)
   rite_c_4=sep3d_ritecf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
  integer function rite_c_5(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 complex,intent(in) :: v(:,:,:,:,:)
   rite_c_5=sep3d_ritecf(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function


  integer function rite_i_1(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 integer,intent(in) :: v(:)
   rite_i_1=sep3d_riteif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_i_2(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 integer,intent(in) :: v(:,:)
   rite_i_2=sep3d_riteif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_i_3(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 integer,intent(in) :: v(:,:,:)
   rite_i_3=sep3d_riteif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
        integer function rite_i_4(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 integer,intent(in) :: v(:,:,:,:)
   rite_i_4=sep3d_riteif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function
  integer function rite_i_5(t1,t2,n,f,j,v,i1,i2,i3,i4)
 character(len=*),intent(in) :: t1,t2
 integer,intent(in) :: i1,i2,i3,n(:),f(:),j(:),i4
 integer,intent(in) :: v(:,:,:,:,:)
   rite_i_5=sep3d_riteif(trim(t1)//C_NULL_CHAR,trim(t2)//C_NULL_CHAR,n,f,j,v,i1,i2,i3,i4)
      end function

end module 
