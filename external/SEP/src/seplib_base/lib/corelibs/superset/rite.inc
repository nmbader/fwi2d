

if(.not. valid_structure(struct)) call erexit("not a valid sep structure\n");
!if(struct%usage=="INPUT") call erexit("can not write out input struct \n");

n=struct%ndims
allocate(nt(n),ft(n),jt(n))

if(present(nwind)) nt=nwind
else nt=-1

if(present(fwind)) ft=fwind
else ft=-1

if(present(jwind)) jt=jwind
else jt=-1


if(present(write_grid)){
	if(struct%file_format!="GRID" && write_grid) 
		call erexit("can not write grid when file_format not grid \n");
	if(struct%file_format=="GRID" && .not. write_grid) 
		call erexit("must write grid when file_format grid \n");
	if(write_grid) wg=1;
	else  wg=0;
}
else{
	if(struct%file_format!="GRID") wg=0
  else wg=1
}



ierr=sep3d_grab_nh(struct%tag,nh)
if(present(write_headers)){
	if(struct%file_format=="REGULAR" && write_headers)
	 call erexit("file_format set to regular can not write out headers \n");
	if(struct%file_format!="REGULAR" && .not. write_headers)
	 call erexit("file_format set not regular must write out headers \n");
	if(write_headers) wh=1
	else wh=0
}
else{
	if(struct%file_format!="REGULAR") wh=1
	else wh=0
}

call fix_vector(struct%n,ft,jt,nt)

if(wd==1 && struct%file_format(1:1)!="G") {
  if(size(tempp)<product(nt)){
    write(0,*) sep_thread_num(),"=thread window size=",product(nt)," buffer size=",size(tempp)
    call seperr("window size > buffer size")
  }
}

	    if(0!=sep3d_rite(tag,struct%tag, nt,ft,jt,tempp, gp,wd,wh,wg))
	    	call erexit("trouble writing out the data \n");

deallocate(nt,ft,jt)
