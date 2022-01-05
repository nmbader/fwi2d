! interface routines  
module sep
  use sep_mod
  use, intrinsic :: iso_c_binding

  implicit none
  integer,private :: ierr
  character(len=10),private :: putch_no,alt_putch
  character(len=256),private :: altfile
  
  

!!$
!!$=head1 NAME
!!$
!!$from_history - grab parameters from history file
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call from_history(par,value,default)>
!!$
!!$C<call from_history(n1,n2,n3,n4,esize,compress)>
!!$
!!$C<call from_history(n,esize)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item par     - char* 
!!$
!!$      Parameter to grab
!!$
!!$=item default - void  
!!$
!!$      (optional) Default value if not found in history file
!!$
!!$=item compress - logical 
!!$
!!$      (optional) [.false.] Compress n2,n3,n4 ->n1
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item value void     
!!$
!!$      Value of parameter in history file
!!$
!!$=item n     C<int(:)>   
!!$
!!$      Dimension of dataset
!!$
!!$=item n1,n2,n3,n4 integer
!!$
!!$      Dimension of dataset
!!$
!!$=item esize integer
!!$
!!$      Esize of dataset (# of bytes per element)
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab parameters from history file
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be float,C<float(:)>,int, or C<int(:)>. Must be same
!!$ for value and default.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<hetch>, L<from_aux>, L<from_param>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface from_history
     module procedure from_history_int
     module procedure from_history_int_array
     module procedure from_history_real
     module procedure from_history_real_array
     module procedure from_history_char
     module procedure from_history_dim
     module procedure from_history_dim_array
     module procedure from_conj_int
     module procedure from_conj_real
     module procedure from_history_double
     module procedure from_history_double_array
  end interface

!!$
!!$=head1 NAME
!!$
!!$from_aux - grab parameters from auilary file
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call from_aux(tag,par,value,default)>
!!$
!!$C<call from_aux(tag,n1,n2,n3,n4,esize,compress)>
!!$
!!$C<call from_aux(tag,n,esize)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item tag     - char*  
!!$
!!$      Tag of file 
!!$
!!$=item par     - char* 
!!$
!!$      Parameter to grab
!!$
!!$=item default - void  
!!$
!!$      (optional) Default value if not found in history file
!!$
!!$=item compress - logical 
!!$
!!$      (optional) [.false.] Compress n2,n3,n4 ->n1
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item value void     
!!$
!!$      Value of parameter in history file
!!$
!!$=item n     C<int(:)>   
!!$
!!$      Dimension of dataset
!!$
!!$=item n1,n2,n3,n4 integer 
!!$
!!$      Dimension of dataset
!!$
!!$=item esize integer       
!!$
!!$      Esize of dataset (# of bytes per element)
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab parameters from history file
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be float,C<float(:)>,int, or C<int(:)>. Must be same
!!$ for value and default.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<auxpar>, L<from_history>, L<from_param>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface from_aux
     module procedure from_aux_bool
     module procedure from_aux_int
     module procedure from_aux_int_array
     module procedure from_aux_real
     module procedure from_aux_real_array
     module procedure from_aux_double
     module procedure from_aux_double_array
     module procedure from_aux_char
     module procedure from_aux_dim
     module procedure from_aux_dim_array
  end interface

  interface from_old
     module procedure from_old_int
     module procedure from_old_real
  end interface

!!$
!!$=head1 NAME
!!$
!!$from_param - grab parameters from comand line parameter
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call from_param(par,value,default)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item par     - char* 
!!$
!!$      Parameter to grab
!!$
!!$=item default - void  
!!$
!!$      (optional) Default value if not found in history file
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item value void     
!!$
!!$      Value of parameter in history file
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab parameters from history file
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be logical, char*, float,C<float(:)>,int, or C<int(:)>. Must be same
!!$ for value and default.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<getch>, L<from_history>, L<from_aux>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface from_param
!     module procedure from_param_long
     module procedure from_param_int
     module procedure from_param_int_array
     module procedure from_param_real
     module procedure from_param_real_array
     module procedure from_param_bool
     module procedure from_param_char
  end interface

  interface from_par
     module procedure from_par_int
     module procedure from_par_int_array
     module procedure from_par_real
     module procedure from_par_real_array
     module procedure from_par_bool
     module procedure from_par_char
  end interface

!!$
!!$=head1 NAME
!!$
!!$from_either - grab parameters from comand line or history parameter
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call from_either(par,value,default)>
!!$
!!$C<call from_either(n1,n2,n3,esize)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item par     - char* 
!!$
!!$      Parameter to grab
!!$
!!$=item default - void  
!!$
!!$      (optional) Default value if not found in history file
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item value void     
!!$
!!$      Value of parameter in history file
!!$
!!$=item n1,n2,n3 int   
!!$
!!$      Datasets dimensions
!!$
!!$=item esize    integer   
!!$
!!$      Datasets esize
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Grab parameters from history file
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be float or int. Must be same
!!$ for value and default.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<fetch>, L<from_history>, L<from_param>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface from_either
     module procedure from_either_int
     module procedure from_either_real
     module procedure from_either_dim
     module procedure from_either_bool
     module procedure from_either_char
  end interface

!!$
!!$=head1 NAME
!!$
!!$to_history - store parameters in history file
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call to_history(par,value,file)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item par     - char* 
!!$
!!$      Parameter to store
!!$
!!$=item value     - void 
!!$
!!$      Value to store
!!$
!!$=item file     - char* 
!!$
!!$      (optional) Tag to store in (defaults to "out")
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Store parameters in sep tag
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be float,C<float(:)>,int, char, or C<int(:)>.
!!$File parameter only availible for int, char and float options.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<putch>, L<to_history>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface to_history
     module procedure to_history_int
     module procedure to_history_int_array
     module procedure to_history_real
     module procedure to_history_real_array
     module procedure to_history_double
     module procedure to_history_double_array
     module procedure to_history_bool
     module procedure to_history_char
  end interface

!!$
!!$=head1 NAME
!!$
!!$sep_read - read a seplib dataset
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep_read(array, file,dim1,esize)>
!!$
!!$C<call sep_read(array, file,dim1,dim2,esize)>
!!$
!!$C<call sep_read(array, file,dim1,dim2,dim3,dim4,esize)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item array     -void  
!!$
!!$      Array to read into
!!$
!!$=item file     - char* 
!!$
!!$      (optional) Tag to read data from (defaults to in)
!!$
!!$=item dim1,dim2,dim3,dim4-int  
!!$
!!$      (optional)  Dimension of datasets (defaults to size of array.)
!!$
!!$=item esize    - int    
!!$
!!$      (optional)   number of bytes per element defaults to size(array)
!!$
!!$=back
!!$
!!$=head1 OUTPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item array     -void  
!!$
!!$      Array to read into
!!$
!!$=back
!!$
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Read in SEPlib dataset
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be C<float(:)>,C<float(:,:)>,C<float(:,:,:)>,C<float(:,:,:.:)>
!!$ C<complex(:,:,:,:)>, C<complex(:,:,:)>, C<complex(:,:)>, and C<complex(:)>
!!$ Number of dims correspond to dimension of the dataset
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep_write>, L<sreed>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface sep_read
     module procedure sep_read_1d
     module procedure sep_read_2d
     module procedure sep_read_3d
     module procedure sep_read_4d
     module procedure sep_read_5d
     module procedure sep_read_complex_1d
     module procedure sep_read_complex_2d 
     module procedure sep_read_complex_3d 
     module procedure sep_read_complex_4d 
     module procedure sep_read_complex_5d
  end interface

!!$
!!$=head1 NAME
!!$
!!$sep_write - write a seplib dataset
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep_write(array, file,dim1,esize)>
!!$
!!$C<call sep_write(array, file,dim1,dim2,esize)>
!!$
!!$C<call sep_write(array, file,dim1,dim2,dim3,dim4,esize)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item array     -void  
!!$
!!$      Array to write from
!!$
!!$=item file     - char* (optional) 
!!$
!!$      Tag to write data to (defaults to out)
!!$
!!$=item dim1,dim2,dim3,dim4-int  (optional)  
!!$
!!$      Dimension of datasets (defaults to size of array)
!!$
!!$=item esize    - int    (optional)   
!!$
!!$      number of bytes per element defaults to size(array)
!!$
!!$=item array     -void  
!!$
!!$      Array to write out
!!$
!!$=back
!!$
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Write out a SEPlib dataset
!!$
!!$=head1 COMMENTS
!!$
!!$Void can be C<float(:)>,C<float(:,:)>,C<float(:,:,:)>,C<float(:,:,:.:)>
!!$ C<complex(:,:,:,:)>, C<complex(:,:,:)>, C<complex(:,:)>, and C<complex(:)>
!!$ Number of dims correspond to dimension of the dataset
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep_read>, L<srite>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  interface sep_write
     module procedure sep_write_1d
     module procedure sep_write_2d
     module procedure sep_write_3d
     module procedure sep_write_4d
     module procedure sep_write_5d
     module procedure sep_write_complex_1d
     module procedure sep_write_complex_2d 
     module procedure sep_write_complex_3d 
     module procedure sep_write_complex_4d 
     module procedure sep_write_complex_5d 
  end interface

  private :: strip_blanks
  private :: num_digits

contains

!!$
!!$=head1 NAME
!!$
!!$sep_init - Initialize a seplib dataset
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep_init(SOURCE)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item source     -char*  
!!$
!!$      Location of source
!!$
!!$=back
!!$
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Initializees SEPlib I/O.
!!$ Calls doc, initpar, noiee, etc
!!$
!!$=head1 SEE ALSO
!!$
!!$L<initpar>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut
  subroutine sep_init (source)
    character (len=*), intent (in), optional :: source
    ! call initpar ()
    character (len=9998) :: tmp
    integer :: nargs, i
    call get_command_argument(0,tmp)
    call init_args(trim(tmp))
    nargs = command_argument_count()
    do i=1,nargs
       call get_command_argument(i,tmp)
       call getch_add_string(trim(tmp)//C_NULL_CHAR)
    end do
    if (present (source)) call doc (source)
  end subroutine sep_init

!!$
!!$=head1 NAME
!!$
!!$sep_close - close a seplib format files
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call sep_close()>
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Calls hclose.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<hclose>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut
  subroutine sep_close ()
    call hclose ()
  end subroutine sep_close
!!$
!!$=head1 NAME
!!$
!!$sep_dimension - returns number of dimensions in dataset
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<ndim=sep_dimension(tag)>
!!$
!!$=head1 INPUT PARAMETERS
!!$
!!$=over 4
!!$
!!$=item file - char* (optional) 
!!$
!!$      tag, defaults to "in"
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$ Gets number of dimension in dataset.
!!$
!!$=head1 SEE ALSO
!!$
!!$L<sep_get_number_data_axes>
!!$
!!$=head1 LIBRARY
!!$
!!$B<sep2df90>
!!$
!!$=cut

  function sep_dimension (file) result (nd)
    integer                                  :: nd
    character (len=*), intent (in), optional :: file

    integer                                  :: stat

    if (present (file)) then
       stat = sep_get_number_data_axes(file,nd)
    else
       stat = sep_get_number_data_axes("in"//C_NULL_CHAR,nd)
    end if
    if (stat /= 0) call erexit("Trouble obtaining number of data axes")
  end function sep_dimension

  function axisname (i,var) result (name)
    character (len=3) :: name
    integer           :: i
    character (len=*) :: var
    optional          :: var
    if (present (var)) then
       if (i < 10) then
          write (name,"(a,i1)") var,i
       else
          write (name,"(a,i2)") var,i
       end if
    else
       if (i < 10) then
          write (name,"(a,i1)") "n",i
       else
          write (name,"(a,i2)") "n",i
       end if
    end if
  end function axisname

  subroutine from_history_dim (n1, n2, n3, n4, esize, compress)
    integer, intent (out)   :: n1, n2, n3, n4
    integer, intent (in)    :: esize
    logical, intent (in)    :: compress
    optional                :: n2, n3, n4, esize, compress

    integer                 :: esz, i1, ni
    integer, parameter      :: maxn1 = 2000000
    logical                 :: axis1

    call from_history ("n1",n1)
    if (present (n2)) call from_history ("n2",n2,1)
    if (present (n3)) call from_history ("n3",n3,1)
    if (present (n4)) call from_history ("n4",n4,1)
    call from_history ("esize", esz)
    if (present (esize)) then
       if (esize == 48 .and. esz /= 4) then
          if (esz /= 8) then
             call erexit('esize /= 4 or 8')
          else
             n1 = n1*2
          end if
       else if (esz /= esize) then
          call erexit("wrong esize")
       end if
    else
       if (esz /= 4) call erexit('esize /= 4')
    end if

    if (.not. present (compress)) return
    if (.not. compress) return
    if (present (n3)) n3 = 1
    if (present (n2)) n2 = 1
    do i1=2, sep_dimension ()
       call from_history (axisname (i1), ni)
       if(.not. present (n2) .or. (n1 * ni < maxn1 .and. axis1)) then
          n1 = n1*ni;
       else
          axis1=.false. 
          n2 = n2*ni;
       end if
    end do
  end subroutine from_history_dim

  subroutine from_history_dim_array (n, esize)
    integer, dimension (:), intent (out) :: n
    integer, intent (in)                 :: esize
    optional                             :: esize

    integer                              :: i1, esz

    do i1=1, size (n)
       call from_history (axisname (i1), n (i1), 1)
    end do

    call from_history ("esize", esz)
    if (present (esize)) then
       if (esize == 48 .and. esz /= 4) then
          if (esz /= 8) then
             call erexit('esize /= 4 or 8')
          else
             n(1) = n(1)*2
          end if
       else if (esz /= esize) then
          call erexit("wrong esize")
       end if
    else
       if (esz /= 4) call erexit('esize /= 4')
    end if
  end subroutine from_history_dim_array

  subroutine from_history_int (name, value, default)
    character (len = *), intent (in) :: name
    integer, intent (out)            :: value
    integer, intent (in), optional   :: default



    if(hetch(name,'i',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing history value: " // name)
       end if
    end if
  end  subroutine from_history_int

  subroutine from_aux_char (file,name, value, default)
    character (len = *), intent (in)  :: name, default,file
    character (len = *), intent (out) :: value
    optional                          :: default

    if(auxpar(name,'s',value,file) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing history value: " // name)
       end if
    end if
  end  subroutine from_aux_char

  subroutine from_history_char (name, value, default)
    character (len = *), intent (in)  :: name, default
    character (len = *), intent (out) :: value
    optional                          :: default


    if(hetch(name,'s',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing history value: " // name)
       end if
    end if
  end  subroutine from_history_char


  subroutine from_history_int_array (name, value, default)
    character (len = *), intent (in)               :: name
    integer, intent (out), dimension (:)           :: value
    integer, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt


    if(hetch(name,'i',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_history_int (varname, value (i), default (i))
          else
             call from_history_int (varname, value (i))
          end if
       end do
    end if
  end  subroutine from_history_int_array

  subroutine from_history_real (name, value, default)
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default


    if(hetch(name,'r',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing history value: " // name)
       end if
    end if
  end  subroutine from_history_real

  subroutine from_history_real_array (name, value, default)
    character (len = *), intent (in)               :: name
    real, intent (out), dimension (:)           :: value
    real, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(hetch(name,'r',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_history_real (varname, value (i), default (i))
          else
             call from_history_real (varname, value (i))
          end if
       end do
    end if
  end  subroutine from_history_real_array

  subroutine from_history_double (name, value, default)
    character (len = *), intent (in)  :: name
    double precision, intent (out)                :: value
    double precision, intent (in), optional       :: default


    if(hetch(name,'g',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing history value: " // name)
       end if
    end if
  end  subroutine from_history_double

  subroutine from_history_double_array (name, value, default)
    character (len = *), intent (in)               :: name
    double precision, intent (out), dimension (:)           :: value
    double precision, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt


    if(hetch(name,'g',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_history_double (varname, value (i), default (i))
          else
             call from_history_double (varname, value (i))
          end if
       end do
    end if
  end  subroutine from_history_double_array

  subroutine from_conj_int (flag, name, value1, def1, value2, def2)
    logical, intent (in)              :: flag
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value1, value2
    integer, intent (in)              :: def1, def2

    if (flag) then
       call from_history (name, value1)
       call from_old (name, value2, def2)
       call to_history (name, value2)
    else
       call from_history (name, value2)
       call from_old (name, value1, def1)
       call to_history (name, value1)
    end if
  end  subroutine from_conj_int

  subroutine from_conj_real (flag, name, value1, def1, value2, def2)
    logical, intent (in)              :: flag
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value1, value2
    real, intent (in)                 :: def1, def2

    if (flag) then
       call from_history (name, value1)
       call from_old (name, value2, def2)
       call to_history (name, value2)
    else
       call from_history (name, value2)
       call from_old (name, value1, def1)
       call to_history (name, value1)
    end if
  end  subroutine from_conj_real

  subroutine from_old_int (name, value, default)
    character (len = *), intent (in) :: name
    integer, intent (out)            :: value
    integer, intent (in), optional   :: default


    if(tetch(name,'i',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing old value: " // name)
       end if
    end if
  end  subroutine from_old_int

  subroutine from_old_real (name, value, default)
    character (len = *), intent (in) :: name
    real, intent (out)               :: value
    real, intent (in), optional      :: default


    if(tetch(name,'r',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing old value: " // name)
       end if
    end if
  end  subroutine from_old_real

  subroutine from_either_dim (n1, n2, n3, esize)
    integer, intent (out)   :: n1, n2, n3
    integer, intent (in)    :: esize
    optional                :: n2, n3, esize

    integer                 :: esz

    call from_either ("n1",n1)
    if (present (n2)) call from_either ("n2",n2)
    if (present (n3)) call from_either ("n3",n3)
    call from_either ("esize", esz)
    if (present (esize)) then
       if (esz /= esize) call erexit("wrong esize")
    else
       if (esz /= 4    ) call erexit('esize /= 4')
    end if
  end subroutine from_either_dim

  subroutine from_either_int (name, value, default)
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value
    integer, intent (in), optional    :: default



    if(fetch(name,'i',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing either value: " // name)
       end if
    end if
    if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From either: " // name,'i',value,altfile)
      else
       ierr= putch("From either: " // name,'i',value)
      end if
    end if
  end  subroutine from_either_int

  subroutine from_either_char (name, value, default)
    character (len = *), intent (in)  :: name, default
    character (len = *), intent (out) :: value
    optional                          :: default



    if(fetch(name,'s',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing either value: " // name)
       end if
    end if
    if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From either: " // name,'s',value,altfile)
      else
       ierr= putch("From either: " // name,'s',value)
      end if
    end if
  end  subroutine from_either_char

  subroutine from_either_real (name, value, default)
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default



    if(fetch(name,'r',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing either value: " // name)
       end if
    end if
       if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From either: " // name,'r',value,altfile)
      else
       ierr= putch("From either: " // name,'r',value)
      end if
    end if
  end  subroutine from_either_real

  subroutine from_either_bool (name, value, default)
    character (len = *), intent (in)  :: name
    logical, intent (out)             :: value
    logical, intent (in), optional    :: default



    if(fetch(name,'l',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing either value: " // name)
       end if
    end if
    if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From either: " // name,'l',value,altfile)
      else
       ierr= putch("From either: " // name,'l',value)
      end if
    end if
  end  subroutine from_either_bool

  subroutine from_param_int (name, value, default)
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value
    integer, intent (in), optional    :: default


    if(getch(name,'i',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
    if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'i',value,altfile)
      else
       ierr= putch("From par: " // name,'i',value)
      end if
    end if

  end  subroutine from_param_int

  subroutine from_param_bool (name, value, default)
    character (len = *), intent (in)  :: name
    logical, intent (out)             :: value
    logical, intent (in), optional    :: default
    logical                           :: val


    if(getch(name,'l',val) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    else
       value = val
    end if
        if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'l',value,altfile)
      else
       ierr= putch("From par: " // name,'l',value)
      end if
    end if

  end  subroutine from_param_bool

  subroutine from_param_int_array (name, value, default)
    character (len = *),    intent (in)           :: name
    integer, dimension (:), intent (out)          :: value
    integer, dimension (:), intent (in), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

   
    if(getch(name,'i',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_param_int (varname, value (i), default (i))
          else
             call from_param_int (varname, value (i))
          end if
       end do
    else
       call to_history_int_array ("From par: " // name,value)
    end if
  end  subroutine from_param_int_array

!  subroutine from_param_long (name, value, default)
!    character (len = *), intent (in)  :: name
!    integer*16, intent (out)                :: value
!    integer*16, intent (in), optional       :: default
!

!
!    if(getch(name,'m',value) == 0) then
!       if (present (default)) then
!          value = default
!       else
!          call erexit("missing parameter value: " // name)
!       end if
!    end if
!       if(putch_no .ne. "no putch") then
!      if(alt_putch(1:3)=="YeS") then
!       ierr= auxputch("From par: " // name,'m',value,altfile)
!      else
!       ierr= putch("From par: " // name,'m',value)
!      end if
!    end if
!
!  end  subroutine from_param_long
  subroutine from_param_real (name, value, default)
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default


    if(getch(name,'r',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
       if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'r',value,altfile)
      else
       ierr= putch("From par: " // name,'r',value)
      end if
    end if

  end  subroutine from_param_real

  subroutine from_param_real_array (name, value, default)
    character (len = *), intent (in)           :: name
    real, dimension (:), intent (out)          :: value
    real, dimension (:), intent (in), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(getch(name,'r',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_param_real (varname, value (i), default (i))
          else
             call from_param_real (varname, value (i))
          end if
       end do
    else
       call to_history_real_array ("From par: " // name,value)
    end if
  end  subroutine from_param_real_array

  subroutine from_param_char (name, value, default)
    character (len = *), intent (in)  :: name, default
    character (len = *), intent (out) :: value
    optional                          :: default



    if(getch(name,'s',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
      if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'s',value,altfile)
      else
       ierr= putch("From par: " // name,'s',value)
      end if
    end if

  end  subroutine from_param_char




!!!!!!!!!!!!OBSOLETE, USE FROM_PARAM

  subroutine from_par_int (name, value, default)
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value
    integer, intent (in), optional    :: default

    if(getch(name,'i',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
    if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'i',value,altfile)
      else
       ierr= putch("From par: " // name,'i',value)
      end if
    end if

  end  subroutine from_par_int

  subroutine from_par_bool (name, value, default)
    character (len = *), intent (in)  :: name
    logical, intent (out)             :: value
    logical, intent (in), optional    :: default
    logical                           :: val



    if(getch(name,'l',val) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    else
       value = val
    end if
        if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'l',value,altfile)
      else
       ierr= putch("From par: " // name,'l',value)
      end if
    end if

  end  subroutine from_par_bool

  subroutine from_par_int_array (name, value, default)
    character (len = *),    intent (in)           :: name
    integer, dimension (:), intent (out)          :: value
    integer, dimension (:), intent (in), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt


    if(getch(name,'i',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_par_int (varname, value (i), default (i))
          else
             call from_par_int (varname, value (i))
          end if
       end do
    else
       call to_history_int_array ("From par: " // name,value)
    end if
  end  subroutine from_par_int_array

  subroutine from_par_real (name, value, default)
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default


    if(getch(name,'r',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
         if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'f',value,altfile)
      else
       ierr= putch("From par: " // name,'f',value)
      end if
    end if

  end  subroutine from_par_real

  subroutine from_par_real_array (name, value, default)
    character (len = *), intent (in)           :: name
    real, dimension (:), intent (out)          :: value
    real, dimension (:), intent (in), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(getch(name,'r',value) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1, size (value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_par_real (varname, value (i), default (i))
          else
             call from_par_real (varname, value (i))
          end if
       end do
    else
       call to_history_real_array ("From par: " // name,value)
    end if
  end  subroutine from_par_real_array

  subroutine from_par_char (name, value, default)
    character (len = *), intent (in)  :: name, default
    character (len = *), intent (out) :: value
    optional                          :: default



    if(getch(name,'s',value) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit("missing parameter value: " // name)
       end if
    end if
        if(putch_no .ne. "no putch") then
      if(alt_putch(1:3)=="YeS") then
       ierr= auxputch("From par: " // name,'s',value,altfile)
      else
       ierr= putch("From par: " // name,'s',value)
      end if
    end if
  end  subroutine from_par_char




























  function exist_file (name) result (test)
    logical                           :: test
    character (len = *), intent (in)  :: name

    character (len = 80)              :: file



    test = (getch (name,'s',file) /= 0)
  end function exist_file

  subroutine to_history_int (name, value, file)
    character (len = *), intent (in)  :: file, name
    integer, intent (in)              :: value
    optional                          :: file

    if (present (file)) then
       ierr= auxputch(name,'i',value,file)
    else
    if(putch_no .ne. "no putch") &
       ierr= putch(name,'i',value)
    end if
  end  subroutine to_history_int

  subroutine to_history_int_array (name, value,tag)
    character (len = *),    intent (in) :: name
    integer, dimension (:), intent (in) :: value
    character(len=*) , optional :: tag

    character (len=128)                 :: tag_out
    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(present(tag)) then
      tag_out=tag
    else
      tag_out="out"
    end if

    varfmt = ' '
    varname = ' '
    do i=1, size(value)
       call num_digits(i,j)
       write(varfmt,'(a,i10,a)')'(a,i',j,')'
       write(varname,varfmt) name,i
       call to_history_int (varname, value (i),tag_out)
    end do

   end subroutine to_history_int_array

  subroutine strip_blanks (string)
    character (len = *)            :: string
    character (len = len (string)) :: new
    integer                        :: i, j, k
    j = 1; k = len (string)
    do i = 1, len (string)
       if (string(i:i) /= ' ') then
          new (j:j) = string (i:i)
          j = j + 1
       else
          new (k:k) = ' '
          k = k - 1
       end if
    end do
    string = new
  end subroutine strip_blanks

  subroutine num_digits(num,ndig)
    integer, intent(in) :: num
    integer, intent(out) :: ndig
    integer j

    ndig=1
    j=10
    do while(ndig < 10)
      if(num < j) return
      ndig = ndig+1
      j = j*10
    enddo
  end subroutine num_digits

  subroutine to_history_real (name, value, file)
    character (len = *), intent (in)  :: file, name
    real, intent (in)                 :: value
    optional                          :: file

    if (present (file)) then
       ierr= auxputch(name,'r',value,file)
    else
    if(putch_no .ne. "no putch") &
       ierr= putch (name,'r',value)
    end if
  end  subroutine to_history_real

 subroutine to_history_real_array (name, value,tag)
    character (len = *), intent (in) :: name
    real, dimension (:), intent (in) :: value
    character (len=*), optional      :: tag
    character (len = 1024)           :: tag_out
    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(present(tag)) then
       tag_out=tag
    else
      tag_out="out"
    end if

    varfmt = ' '
    varname = ' '
    do i=1, size(value)
       call num_digits(i,j)
       write(varfmt,'(a,i10,a)')'(a,i',j,')'
       write(varname,varfmt) name,i
       call to_history_real (varname, value (i),tag_out)
    end do
  end  subroutine to_history_real_array

  subroutine to_history_double (name, value, file)
    character (len = *), intent (in)  :: file, name
    double precision, intent (in)     :: value
    optional                          :: file

    if (present (file)) then
       ierr= auxputch(name,'g',value,file)
    else
    if(putch_no .ne. "no putch") &
       ierr= putch (name,'g',value)
    end if
  end  subroutine to_history_double

 subroutine to_history_double_array (name, value,tag)
    character (len = *), intent (in) :: name
    double precision, dimension (:), intent (in) :: value
    character (len=*), optional     :: tag
    character (len = 130)           :: tag_out
    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(present(tag)) then
       tag_out=tag
    else
      tag_out="out"
    end if
    varfmt = ' '
    varname = ' '
    do i=1,size(value)
       call num_digits(i,j)
       write(varfmt,'(a,i10,a)')'(a,i',j,')'
       write(varname,varfmt) name,i
       call to_history_double (varname, value (i),tag_out)
    enddo
  end  subroutine to_history_double_array

  subroutine to_history_bool (name, value, file)
    character (len = *), intent (in)  :: file, name
    logical, intent (in)              :: value
    optional                          :: file

    if (present (file)) then
       ierr= auxputch(name,'l',value,file)
    else
    if(putch_no .ne. "no putch") &
       ierr= putch (name,'l',value)
    end if
  end  subroutine to_history_bool

  subroutine to_history_char (name, value, file)
    character (len = *), intent (in)  :: file, name, value
    optional                          :: file

    if (present (file)) then
       ierr= auxputch(name,'s',value,file)
    else
    if(putch_no .ne. "no putch") &
       ierr= putch (name,'s',value)
    end if
  end  subroutine to_history_char

  subroutine from_aux_dim (file, n1, n2, n3, esize)
    character (len = *), intent (in)   :: file
    integer, intent (out)              :: n1, n2, n3
    integer, intent (in)               :: esize
    optional                           :: n2, n3, esize

    integer                            :: esz

    call from_aux (file, "n1",n1)
    if (present (n2)) call from_aux (file, "n2",n2,1)
    if (present (n3)) call from_aux (file, "n3",n3,1)
    call from_aux (file, "esize", esz)
    if (present (esize)) then
       if (esz /= esize) call erexit("wrong esize")
    else
       if (esz /= 4    ) call erexit('esize /= 4')
    end if
  end subroutine from_aux_dim

  subroutine from_aux_dim_array (file, n, esize)
    character (len=*), intent (in)       :: file
    integer, dimension (:), intent (out) :: n
    integer, intent (in)                 :: esize
    optional                             :: esize

    integer                              :: i1, esz

    do i1=1, size (n)
       call from_aux (file, axisname (i1), n (i1), 1)
    end do

    call from_aux (file, "esize", esz)
    if (present (esize)) then
       if (esize == 48 .and. esz /= 4) then
          if (esz /= 8) then
             call erexit('esize /= 4 or 8')
          else
             n(1) = n(1)*2
          end if
       else if (esz /= esize) then
          call erexit("wrong esize")
       end if
    else
       if (esz /= 4) call erexit('esize /= 4')
    end if
  end subroutine from_aux_dim_array

  subroutine from_aux_bool (file, name, value, default)
    character (len = *), intent (in)   :: file, name
    logical, intent (out)              :: value
    logical, intent (in), optional     :: default

    if(auxpar(name,'1',value,file) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit ("missing history value: " // name // &
          " from file: " // file )
       end if
    end if
  end  subroutine from_aux_bool

  subroutine from_aux_int (file, name, value, default)
    character (len = *), intent (in)   :: file, name
    integer, intent (out)              :: value
    integer, intent (in), optional     :: default



    if(auxpar(name,'i',value,file) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit ("missing history value: " // name // &
          " from file: " // file )
       end if
    end if
    !    ierr= putch ("From aux(" // file // "): " // file // "_" // name, &
    !         'i',value)
  end  subroutine from_aux_int

  subroutine from_aux_int_array (file, name, value, default)
    character (len = *), intent (in)               :: file, name
    integer, intent (out), dimension (:)           :: value
    integer, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt

    if(auxpar(name,'i',value,file) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1,size(value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_aux_int (file, varname, value (i), default (i))
          else
             call from_aux_int (file, varname, value (i))
          end if
       end do
    else
       !       call to_history_int_array (&
       !       "From Aux(" // file // "): " // file // "_" // name,value)
    end if
  end  subroutine from_aux_int_array

  subroutine from_aux_real (file, name, value, default)
    character (len = *), intent (in)  :: file, name
    real, intent (out)                :: value
    real, intent (in), optional       :: default



    if(auxpar(name,'r',value,file) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit ("missing history value: " // name // &
          " from file: " // file )
       end if
    end if
    !    ierr= putch ("From aux(" // file // "): " // file // "_" // name, &
    !         'r',value)
  end  subroutine from_aux_real

  subroutine from_aux_real_array (file, name, value, default)
    character (len = *), intent (in)               :: file, name
    real, intent (out), dimension (:)           :: value
    real, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt


    if(auxpar(name,'r',value,file) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1,size(value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_aux_real (file, varname, value (i), default (i))
          else
             call from_aux_real (file, varname, value (i))
          end if
       end do
    else
       !    	call to_history_real_array (&
       !        "From aux(" // file // "): " // file // "_" // name,value)
    end if
  end  subroutine from_aux_real_array

  subroutine from_aux_double (file, name, value, default)
    character (len = *), intent (in)  :: file, name
    double precision, intent (out)                :: value
    double precision, intent (in), optional       :: default



    if(auxpar(name,'g',value,file) == 0) then
       if (present (default)) then
          value = default
       else
          call erexit ("missing history value: " // name // &
          " from file: " // file )
       end if
    end if
    !    ierr= putch ("From aux(" // file // "): " // file // "_" // name, &
    !         'g',value)
  end  subroutine from_aux_double

  subroutine from_aux_double_array (file, name, value, default)
    character (len = *), intent (in)               :: file, name
    double precision, intent (out), dimension (:)           :: value
    double precision, intent (in),  dimension (:), optional :: default

    integer                                        :: i,j
    character (len = len(name) + 10)               :: varname
    character (len=15)                             :: varfmt


    if(auxpar(name,'g',value,file) == 0) then
       varfmt = ' '
       varname = ' '
       do i=1,size(value)
          call num_digits(i,j)
          write(varfmt,'(a,i10,a)')'(a,i',j,')'
          write(varname,varfmt) name,i
          if (present (default)) then
             call from_aux_double (file, varname, value (i), default (i))
          else
             call from_aux_double (file, varname, value (i))
          end if
       end do
    else
       !    	call to_history_double_array (&
       !        "From aux(" // file // "): " // file // "_" // name,value)
    end if
  end  subroutine from_aux_double_array

  subroutine sep_read_1d (array, file, dim, esize)
    real, dimension (:), intent (inout)         :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim, esize

    integer                                     :: bites

    if (present (dim)) then
       bites = dim
    else
       bites = size (array)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_1d

  subroutine sep_read_2d (array, file, dim1, dim2, esize)
    real, dimension (:,:), intent (inout)       :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, esize

    integer                                     :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_2d

  subroutine sep_read_3d (array, file, dim1, dim2, dim3, esize)
    real, dimension (:,:,:), intent (inout)     :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, dim3, esize

    integer                                     :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_3d

  subroutine sep_read_4d (array, file, dim1, dim2, dim3, dim4, esize)
    real, dimension (:,:,:,:), intent (inout)  :: array 
    character (len = *), intent (in), optional :: file
    integer, intent (in), optional             :: dim1,dim2,dim3,dim4, esize

    integer                                    :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (dim4)) then
       bites = bites * dim4
    else
       bites = bites * size (array,4)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_4d

  subroutine sep_read_5d (array, file, dim1, dim2, dim3, dim4, dim5, esize)
    real, dimension (:,:,:,:,:), intent (inout)  :: array 
    character (len = *), intent (in), optional :: file
    integer, intent (in), optional             :: dim1,dim2,dim3,dim4,dim5, esize

    integer                                    :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (dim4)) then
       bites = bites * dim4
    else
       bites = bites * size (array,4)
    end if

    if (present (dim5)) then
       bites = bites * dim5
    else
       bites = bites * size (array,5)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_5d

  subroutine sep_read_complex_1d (array, file, dim, esize)
    complex, dimension (:), intent (inout)      :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim, esize

    integer                                     :: bites

    if (present (dim)) then
       bites = dim
    else
       bites = size (array)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_complex_1d

  subroutine sep_read_complex_2d (array, file, dim1, dim2, esize)
    complex, dimension (:,:), intent (inout)    :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, esize

    integer     :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_complex_2d

  subroutine sep_read_complex_3d (array, file, dim1, dim2, dim3, esize)
    complex, dimension (:,:,:), intent (inout)  :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, dim3, esize

    integer                                     :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_complex_3d

  subroutine sep_read_complex_4d (array, file, dim1, dim2, dim3, dim4, esize)
    complex, dimension (:,:,:,:), intent (inout)  :: array 
    character (len = *), intent (in), optional :: file
    integer, intent (in), optional             :: dim1,dim2,dim3,dim4, esize

    integer                                    :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (dim4)) then
       bites = bites * dim4
    else
       bites = bites * size (array,4)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_complex_4d

  subroutine sep_read_complex_5d (array, file, dim1, dim2, dim3, dim4, dim5, esize)
    complex, dimension (:,:,:,:,:), intent (inout)  :: array 
    character (len = *), intent (in), optional :: file
    integer, intent (in), optional             :: dim1,dim2,dim3,dim4,dim5, esize

    integer                                    :: bites

    if (present (dim1)) then
       bites = dim1
    else
       bites = size (array,1)
    end if

    if (present (dim2)) then
       bites = bites * dim2
    else
       bites = bites * size (array,2)
    end if

    if (present (dim3)) then
       bites = bites * dim3
    else
       bites = bites * size (array,3)
    end if

    if (present (dim4)) then
       bites = bites * dim4
    else
       bites = bites * size (array,4)
    end if

    if (present (dim5)) then
       bites = bites * dim5
    else
       bites = bites * size (array,5)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr= sreed(file, array, bites)
    else
       ierr= sreed('in', array, bites)
    end if
  end subroutine sep_read_complex_5d

  subroutine sep_write_1d (array, file, dim, esize)
    real, dimension (:), intent (in)            :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim, esize

    integer      :: bites

    if (present (dim)) then
       bites = dim
    else
       bites = size (array)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 4
    end if

    if (present (file)) then
       ierr = srite(file, array, bites)
    else
       ierr = srite("out", array, bites)
    end if
  end subroutine sep_write_1d

  subroutine sep_write_complex_1d (array, file, dim, esize)
    complex, dimension (:)         :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim, esize

    integer      :: bites

    if (present (dim)) then
       bites = dim
    else
       bites = size (array)
    end if

    if (present (esize)) then
       bites = bites * esize
    else
       bites = bites * 8
    end if

    if (present (file)) then
       ierr = srite(file, array, bites)
    else
       ierr = srite("out", array, bites)
    end if
  end subroutine sep_write_complex_1d

  subroutine sep_write_2d (array, file, dim1, dim2, esize)
    real, dimension (:,:), intent (in)          :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, esize

    integer                                     :: i, n1, n2, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 4
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    do i = 1, n2
       call sep_write_1d (array (:,i), f, n1, e)
    end do
  end subroutine sep_write_2d

  subroutine sep_write_complex_2d (array, file, dim1, dim2, esize)
    complex, dimension (:,:), intent (in)       :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, esize

    integer                                     :: i, n1, n2, e
    character (len = 80)                        :: f

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 8
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    do i = 1, n2
       call sep_write_complex_1d (array (:,i), f, n1, e)
    end do
  end subroutine sep_write_complex_2d

  subroutine sep_write_complex_3d (array, file, dim1, dim2, dim3, esize)
    complex, dimension (:,:,:), intent (in)     :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, dim3, esize

    integer                                     :: i3, i2, n1, n2, n3, e
    character (len = 80)                        :: f

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 8
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    do i3 = 1, n3
       do i2 = 1, n2
          call sep_write_complex_1d (array (:,i2,i3), f, n1, e)
       end do
    end do
  end subroutine sep_write_complex_3d

  subroutine sep_write_complex_4d (array, file, dim1, dim2, dim3, dim4, esize)
    complex, dimension (:,:,:,:), intent (in)   :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1,dim2,dim3,dim4, esize

    integer                                     :: i2,i3,i4, n1,n2,n3,n4, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 8
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    if (present (dim4)) then
       n4 = dim4
    else
       n4 = size (array, 4)
    end if

    do i4 = 1, n4
       do i3 = 1, n3
          do i2 = 1, n2
             call sep_write_complex_1d (array (:,i2,i3,i4), f, n1, e)
          end do
       end do
    end do
  end subroutine sep_write_complex_4d

  subroutine sep_write_complex_5d (array, file, dim1, dim2, dim3, dim4, dim5, esize)
    complex, dimension (:,:,:,:,:), intent (in)   :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1,dim2,dim3,dim4,dim5, esize

    integer                                     :: i2,i3,i4,i5, n1,n2,n3,n4,n5, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 8
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    if (present (dim4)) then
       n4 = dim4
    else
       n4 = size (array, 4)
    end if

    if (present (dim5)) then
       n5 = dim5
    else
       n5 = size (array, 5)
    end if

    do i5 = 1, n5
       do i4 = 1, n4
          do i3 = 1, n3
             do i2 = 1, n2
                call sep_write_complex_1d (array (:,i2,i3,i4,i5), f, n1, e)
             end do
          end do
       end do
    end do
  end subroutine sep_write_complex_5d

  subroutine sep_write_3d (array, file, dim1, dim2, dim3, esize)
    real, dimension (:,:,:), intent (in)        :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1, dim2, dim3, esize

    integer                                     :: i2, i3, n1, n2, n3, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 4
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    do i3 = 1, n3
       do i2 = 1, n2
          call sep_write_1d (array (:,i2,i3), f, n1, e)
       end do
    end do
  end subroutine sep_write_3d

  subroutine sep_write_4d (array, file, dim1, dim2, dim3, dim4, esize)
    real, dimension (:,:,:,:), intent (in)      :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1,dim2,dim3,dim4, esize

    integer                                     :: i2,i3,i4, n1,n2,n3,n4, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 4
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    if (present (dim4)) then
       n4 = dim4
    else
       n4 = size (array, 4)
    end if

    do i4 = 1, n4
       do i3 = 1, n3
          do i2 = 1, n2
             call sep_write_1d (array (:,i2,i3,i4), f, n1, e)
          end do
       end do
    end do
  end subroutine sep_write_4d

  subroutine sep_write_5d (array, file, dim1, dim2, dim3, dim4, dim5, esize)
    real, dimension (:,:,:,:,:), intent (in)      :: array 
    character (len = *), intent (in), optional  :: file
    integer, intent (in), optional              :: dim1,dim2,dim3,dim4,dim5, esize

    integer                                     :: i2,i3,i4,i5, n1,n2,n3,n4,n5, e
    character (len = 80)                        :: f 

    if (present (file)) then
       f = file
    else
       f = "out"
    end if

    if (present (esize)) then
       e = esize
    else
       e = 4
    end if

    if (present (dim1)) then
       n1 = dim1
    else
       n1 = size (array, 1)
    end if

    if (present (dim2)) then
       n2 = dim2
    else
       n2 = size (array, 2)
    end if

    if (present (dim3)) then
       n3 = dim3
    else
       n3 = size (array, 3)
    end if

    if (present (dim4)) then
       n4 = dim4
    else
       n4 = size (array, 4)
    end if

    if (present (dim5)) then
       n5 = dim5
    else
       n5 = size (array, 5)
    end if

    do i5 = 1, n5
       do i4 = 1, n4
          do i3 = 1, n3
             do i2 = 1, n2
                call sep_write_1d (array (:,i2,i3,i4,i5), f, n1, e)
             end do
          end do
       end do
    end do
  end subroutine sep_write_5d



 subroutine set_yes_putch()
	putch_no="yes it"
 end subroutine

 subroutine set_no_putch()
	putch_no="no putch"
 end subroutine

 subroutine set_alternate_putch(tag)
  character(len=*)  :: tag
  alt_putch="YeS"
  altfile=tag
 end subroutine



end module sep
