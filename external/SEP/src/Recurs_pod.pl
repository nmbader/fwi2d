#!/usr/bin/perl


#if($#ARGV < 1) {
#
#die "$0: $selfdoc ";
#}

require "ctime.pl";


#MAN PAGE CONFIG
$center="SEPlib Manual Pages";
$date=&ctime(time);
$date=~s/\s+$//;
$release="6.0";



#MANUAL PAGE LOCATION
$default_man="1";
$man{"seplib_base"}=3;

$pwd=`pwd`;
chop $pwd;

@subs=(
"seplib_base",
"seplib_prog",
"vplot",
"interact");

$pod_list=join(":",@subs);
$pod_path1="--podpath=.:$pod_list";
$pod_list=join(":../",@subs);
$pod_path2="--podpath=..:../$pod_list";
$pod_base="$pwd/docs/pod";



$html_dir_base="$pwd/docs/html";
$html_base="--htmlroot=file:/usr/local/src/our/sep/html";
$html_base="--htmlroot=file:/net/koko/scr1/bob/seplib-6.3.0/html";
$html_base1="--htmlroot=.";
$html_base2="--htmlroot=.";


$man_root="$pwd/docs/man";


#LIBRARY DESCRIPTIONS
$lib_desc{"sep"}="SEPlib Base Library";
$lib_desc{"geef90"}="Fortran 90 GEE Operator Library";
$lib_desc{"supersetf90"}="Fortran 90 interface for handling sep3d datasets";
$lib_desc{"septravel"}="Functions that calculate traveltimes";
$lib_desc{"sepvelan"}="Functions that deal with velocity";
$lib_desc{"sepvelanf"}="F77 functions that deal with velocity";
$lib_desc{"sepvelanf90"}="F90 fFunctions that deal with velocity";
$lib_desc{"sepfilter"}="Functions that deal with filtering";
$lib_desc{"sepfilterf90"}="Fortran90 Functions that deal with filtering";
$lib_desc{"sep2df90"}="Fortran90 interface for handling sep2d datasets";
$lib_desc{"sepaux"}="Useful functions that don't fit anywher else";
$lib_desc{"sepauxf90"}="Useful F90 functions that don't fit anywher else";
$lib_desc{"sepfft"}="Functions that do FFTs";
$lib_desc{"sepmath"}="Math functions (especially Complex C)";
$lib_desc{"sepmathf90"}="F90 Math functions";

$lib_desc{"sep3d"}="Seplib base libary";
$lib_desc{"sepsu"}="Functions that allow Seplib and Su to interface";
$lib_desc{"sepocf90"}="Functions that allow out of core optimization";
$lib_desc{"sepweif90"}="A library for performing  wave equation migration";



#CATOGORY DESCRIPTIONS
$cat_desc{"sep_graphics"}="SEPlib Plotting programs";
$cat_desc{"graphics/vplot/util/shells"}="Useful tools for common vplot commands";
$cat_desc{"graphics/vplot/filters"}="Vplot filters";
$cat_desc{"graphics/vplot"}="Generic Vplot";
$cat_desc{"graphics/Rickmovie"}="Rickmovie";
$cat_desc{"util/cube"}="Programs that take advantage of the SEP hypercube";
$cat_desc{"util/headers"}="Programs that work on SEP3d headers";
$cat_desc{"util/info"}="Programs that provide info about SEP datasets";
$cat_desc{"util/vector"}="Programs that do vector operations on SEP datasets";
$cat_desc{"util/unix"}="Programs that perform unix-like operations";
$cat_desc{"converters"}="Programs that convert to and from seplib types";
$cat_desc{"tools"}="Useful tools that we use when programing";
$cat_desc{"seis/travel"}="Travel time and rays programs \n";
$cat_desc{"seis/filter"}="Filtering programs \n";
$cat_desc{"seis/image"}="Imaging programs \n";
$cat_desc{"seis/model"}="Modeling programs \n";
$cat_desc{"seis/velan"}="Programs that deal with velocity \n";
$cat_desc{"interact"}="Interactive SEPlib programs \n";



if($#ARGV != 0) {
	print STDERR "$0: [mode]
	
	mode=pod	
	mode=man
	mode=html
	mode=tex
";
die;
}
if($ARGV[0] eq "pod") { &create_pod;}
elsif($ARGV[0] eq "man") { &create_man;}
elsif($ARGV[0] eq "html") { &create_html;}
elsif($ARGV[0] eq "tex") { &create_tex;}
else{
	print STDERR "$0: [mode]
	
	mode=pod	
	mode=man
	mode=html
	mode=tex
";
die;
}


sub create_pod{

$ncat=0; $nlib=0;
system("rm -rf $pod_base/*");
for($i=0; $i < @subs; $i++){

@ends=("F90","r90","F90","c","pl");
my $end = join("\$|",@ends);

print "END $end \n";

$sub=$subs[$i];
open(FIND, "find $sub -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
FILE:
while ($file = <FIND>) {
  chomp $file;
  next FILE unless -T $file;
  if($file=~/$end$/){
	$pod_dir="$pod_base/$sub/";
	system("mkdir -p $pod_dir");

  
	open(IT,$file);
		$buffer="";
		$start=0;
		$look=-1;
		while ($line=<IT>){
  	if ($start == 0 && $line=~/^\s*\=head1\s+NAME/) {
			$start = 2;
			$comment="";
			$look=0;
			$cat=0; $lib=0;
		}
  	elsif ($start == 0 && $line=~/^(.*)\=head1 NAME/) {
			$test=$1;
			$comment = quotemeta $1; 
			if($test =~/^!/ || $test =~ /^#/){
 				$start = 1;
				$look=0;
			}
			else { print STDERR "bad comment $comment: $start \n";}
 	  }
		if($start==1){ $line=~s/$comment//;}
		if($look==1 && $line=~/\s*(\w+)\s+.+$/){ $base=$1;
			$pod_file="$pod_dir".$base.".pod";
			open(OUT,">$pod_file");
			$look=-1;
			if($line=~/^\s*(\S+)\s*-\s*(\S.+\S)\s*$/){$desc=$2;$line="$1 - $2\n";}
      $new=1;
      if($exist{$base}==1) {$new=0;}
      $exist{$base}=1;
		}
		elsif($look==0) {$look=1;}

    

		if($line=~/^\s*$/){ $line="\n";}
		if($line=~/^=head1\s+(\S.+\S)\s*$/){ $line="=head1 $1\n\n";}
		if($start >0 ) {$buffer.=$line;}

		if($lib==1){
			if($line=~/^\s*$/){}
			elsif($line=~/\s*B\S*<\s*(.+)\s*>/){
				$lib=2; $libra=$1;
       if($new==1){ &add_lib;}
			}
			elsif($line=~/\s*(\S*.+\S)\s*/){
				$lib=2; $libra=$1;
       if($new==1){ &add_lib;}
			}
		}

		if($cat==1){
			if($line=~/^\s*$/){}
			elsif($line=~/\s*B\S*<\s*(.+)\s*>/){
				$cat=2; $catog=$1;
       if($new==1){ &add_cat;}
			}
			elsif($line=~/\s*(\S*.+\S)\s*/){
				$cat=2; $catog=$1;
       if($new==1){ &add_cat;}
			}
		}

		if($line=~/^\s*=head1\s*CATEGORY/){ $cat=1;}
		if($line=~/^\s*=head1\s*LIBRARY/){ $lib=1;}


		if($line =~/^\=cut/){ 
			print OUT "$buffer\n";
			$buffer="";	
			close(OUT); 
			$look=-1;
			$start=0;
			if($lib==0 && $cat==0){ 
				print STDERR "$file:$base no category or library \n";
			}
			elsif($lib==1){
				print STDERR "couldn't find library name: $file:$base \n";
			}
			elsif($cat==1){
				print STDERR "couldn't find category name  $file:$base\n";
			}
#			print STDERR "$base:$file \n";
		}
	}
	if($start!=0){
		print STDERR "$file:$base ended before cut \n";
	}
}
else{
#  print STDERR "NO  $file \n";
}
}


}
&finish_it;
$dir=$pod_base;&create_make;
}

sub create_make{


open(FIND, "find $dir -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
FILE:
while ($file = <FIND>) {
  chomp $file;
	if(-d $file){
		if($file=~/$dir\/(.+)$/){	$part=$1;
			$subs.=" $part";
		}
		else{print STDERR "trouble reading dir: $file -$dir \n";}
	}
	elsif(-T $file){
		if($file=~/$dir\/(.+)\/(.+)$/){ $dirit=$1; $file=$2;
			$file_list{$dirit}.="$file ";
		}
		elsif($file=~/$dir\/(.+)$/){ $file=$1;
			$file_list{"main"}.="$file ";
		}
	}
	else{ print STDERR "unrecognized $file \n";}
}

@dirs=split(" ",$subs);

open(MAIN,">$dir/Makefile.am");
print MAIN "include \${top_srcdir}/include/SEP.aux.rules
SUBDIRS=$subs
EXTRA_DIST=".$file_list{"main"}."
\n";
close(MAIN);



foreach $directory (@dirs) {
	if($directory=~/^.+([0-9])$/){ $suf=$1;}
	open(MAIN,">$dir/$directory/Makefile.am");
$count=0;
$big=0;
if($ARGV[0] eq "man"){
	@files=split(" ",$file_list{$directory});
	$line1="include \${top_srcdir}/include/SEP.aux.rules
NUM0=";
	$line2="MN0=";
	$line3="FINAL0=";
  $line4="MN= \${MN0}";
  $line5="FINAL= \${FINAL0}";
  $line6="NUM= \${NUM0}";

	foreach $file (@files){
    $count++;
		if($file=~/.+\.mn/) {$line2.=" $file";}
		elsif($file=~/.+\.[0-9]/) {
			$line1.=" $file";
			$line3.="\${mandir}/${directory}/$file ";
			}
		else{print STDERR "ignoring $file \n";}
   if($count==30){ $count=0;$big++;
    $line1.="\nNUM$big=";
    $line2.="\nMN$big=";
    $line3.="\nFINAL$big=";
    $line4.=" \${MN$big}";
    $line5.=" \${FINAL$big}";
    $line6.=" \${NUM$big}";
   }
	}
	print MAIN  "
$line1

$line2

$line3

$line4

$line5

$line6

if HAVE_MAN_FMT

MANUA=\${MN}

\${mandir}/${directory}/%.$suf:	\${srcdir}/%.mn
	cp \${srcdir}/\$*.mn  \$@
	


else

MANUA=\${NUM}

\${mandir}/${directory}/%:	\${srcdir}/%
	cp  \${srcdir}/\$* \$@


endif

dirs:
	\$(mkinstalldirs) \${mandir}/${directory}

install-data-local: dirs \${FINAL}
	

uninstall-local:	
	\${TOUCH} \${FINAL}
	\${RM} \${FINAL}

EXTRA_DIST=\${NUM} \${MN}
\n";
}
else{
	print MAIN "include \${top_srcdir}/include/SEP.aux.rules
EXTRA_DIST=$file_list{$directory}
\n";
}

close(MAIN);
}

}

sub create_tex{
open(FIND, "find docs/pod -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
FILE:
while ($file = <FIND>) {
  chomp $file;
	print STDERR  " $file \n";

if($file=~/^pod\/(.+)\/(.+)\.pod/){
	$dir=$1; $name=$2;
	$pod_file=$file;
	

	print STDERR "on $name \n";
	$command="pod2latex   $pod_file";
	system($command) == 0
         or  print STDERR "NAME=$name \n";


}

}


}



sub create_man{


system("rm -rf $man_root/*");

print STDERR "pod_base=$pod_base\n";
open(FIND, "find $pod_base -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
$count=0;
FILE:
while ($file = <FIND>) {
  chomp $file;

if($file=~/^$pod_base\/(.+)\/(.+)\.pod/){
	print STDERR "$file =file \n";
	$dir=$1; $name=$2;
	if($man{$dir}){ $man_suf=$man{$dir};}
	else { $man_suf=$default_man;}
	$man_dir="$man_root/man$man_suf/";
	system("mkdir -p $man_dir");
	$pod_file=$file;
	$mn_file="$man_dir/$name.mn";
	$man_file="$man_dir/$name.$man_suf";

	
	$command="pod2man  $pod_file --section=$man_suf --release=$release --center=\"$center\" --date=\"$date\"  > $mn_file";

	system($command) == 0
         or  print STDERR "NAME=$name \n";

	$command="neqn  $mn_file | tbl | nroff -man > $man_file";

	system($command) == 0
         or  print STDERR "NAME=$name \n";



}

}
$dir=$man_root;&create_make;


}







sub old_stuff{


for($i=0; $i < @subs; $i++){

	$sub=$subs[$i];
	$pod_dir="pod/$sub/";
	$man_dir="man/$sub/";
	$html_dir="html/$sub/";
	$html_dir="/sepwww/pub/sep/bob/temp3/$sub/";
	$tex_dir="tex/$sub/";
	if($man{$sub}){ $man_suf=$man{$sub};}
	else { $man_suf=$default_man;}
	system("mkdir -p $pod_dir");
	system("mkdir -p $man_dir");
	system("mkdir -p $html_dir");
	system("mkdir -p $tex_dir");
	

open(FIND, "find $sub -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
FILE:
while ($file = <FIND>) {
  chomp $file;
  next FILE unless -T $file;

	open(IT,$file);
		$buffer="";
		$start=0;
		$look=-1;
		while ($line=<IT>){
  	if ($line=~/^\s*\=head1\s+NAME/) {
			$start = 2;
			$comment="";
			$look=0;
		}
  	elsif ($line=~/^(.*)\=head1 NAME/) {
			$start = 1;
			$comment = quotemeta $1; 
			$look=0;
 	  }
		if($start==1){ $line=~s/$comment//;}
		if($look==1 && $line=~/\s*(\w+)\s+.+$/){ $base=$1;
			$pod_file="$pod_dir".$base.".pod";
			$man_file="$man_dir".$base.".$man_suf";
			$html_file="$html_dir".$base.".html";
			$tex_file="$tex_dir".$base.".tex";
			open(OUT,">$pod_file");
			$look=-1;
			if($line=~/^\s*(\S+)\s*-\s*(\S.+\S)\s*$/){$line="$1 - $2\n";}
		}
		elsif($look==0) {$look=1;}

		if($line=~/^\s*$/){ $line="\n";}
		if($line=~/^=head1\s+(\S.+\S)\s*$/){ $line="=head1 $1\n\n";}
		if($start >0 ) {$buffer.=$line;}
		if($line =~/^\=cut/){ 
			print OUT "$buffer\n";
			$buffer="";	
			close(OUT); 
			$command="pod2man  $pod_file --section=$man_suf --release=$release --center=\"$center\" --date=\"$date\"  > $man_file";
			system($command);
			$command="pod2html  < $pod_file  > $html_file";
			system($command);
			$command="pod2latex   $pod_file";
			system($command);
			$look=-1;
			$start=0;
		}
	}
	close(IT);
}
}
}

sub add_cat{

if($cat_doc{$catog}){ #catogory exists

	$cat_doc{$catog}.="L<$base> - $desc\n\n";

}
else{ 
	$cat_list[$ncat]=$catog;
	$ncat++;
	if($cat_desc{$catog}){ $info=$cat_desc{$catog};}
	else {$info="UNDEFINED";$cat_desc{$catog}="UNDEFINED";}
	$cat_doc{$catog}="

=head1 NAME

$catog - $info

=head1 SYNOPSIS

A SEPlib library

=head1 DESCRIPTION

$info

=head1 PROGRAMS


L<$base> - $info



";
}

}

sub add_lib{

if($lib_doc{$libra}){ #catogory exists

	$lib_doc{$libra}.="L<$base> - $desc\n\n";

}
else{ 
	$lib_list[$nlib]=$libra;
	$nlib++;
	if($lib_desc{$libra}){ $info=$lib_desc{$libra};}
	else {$info="UNDEFINED";$lib_desc{$libra}="UNDEFINED";}
	$lib_doc{$libra}="

=head1 NAME

$libra - $info

=head1 SYNOPSIS

A SEPlib library

=head1 DESCRIPTION

$info

=head1 FUNCTIONS



L<$base> -  $desc



";
}

}

sub finish_it{


for($i=0; $i < $nlib; $i++){
	open(AAA,">$pod_base/$lib_list[$i].pod");
	print AAA $lib_doc{$lib_list[$i]};
	print AAA "=cut ";
	$lib_it.="L<$lib_list[$i]> - $lib_desc{$lib_list[$i]} \n\n";
	close(AAA);
}

for($i=0; $i < $ncat; $i++){
	$convert=$cat_list[$i];
	@ss=split("/",$convert);
	$convert=join("_",@ss);
	open(AAA,">$pod_base/$convert.pod");
	print AAA $cat_doc{$cat_list[$i]};
	print AAA "\n","=cut ";
	$cat_it.="L<$convert> - $cat_desc{$cat_list[$i]}\n\n";
	close(AAA);
}


$main="

=head1 NAME

SEPlib - SEP programs and libraries



=head1  SYNOPSIS

SEPlib utilities and libaries

=head1 DESCRIPTION

If your wanting to know more about SEPlib programs or
library functions follow the following links/ view
the following man pages

=head1 PROGRAMS


$cat_it


=head1 LIBRARIES


$lib_it


=cut

";

open(MAIN,">$pod_base/seplib.pod");
print  MAIN $main;
close(MAIN);

}

sub create_html{

print STDERR "check $html_dir_base \n";
system("rm -rf $html_dir_base/*");
open(FIND, "find $pod_base -not \\( \\( -name .svn -prune \\) -o \\( -name .deps -prune \\) -o \\( -name .libs -prune \\) \\) -type f -print |") or die "Can't run find: $!\n";
FILE:
while ($file = <FIND>) {
  chomp $file;

$run=0;
if($file=~/^$pod_base\/(.+)\/(.+)\.pod/){
	print STDERR "$file =file \n";
	$dir=$1; $name=$2;
	$html_dir="$html_dir_base/$dir";
#	$input="$dir/$name.pod";
	$input="$name.pod";
	system("mkdir -p $html_dir");
	$output="$html_dir/$name.html";
	$command="cd $pod_base/$dir; pod2html   <$input > $output $pod_path2 $html_base2";
  print STDERR "$command \n";
	system($command) == 0
         or  print STDERR "NAME=$input \n";
}
elsif($file=~/^$pod_base\/(.+)\.pod/){
	$name=$1;
	$input="$name.pod";
	$html_dir="$html_dir_base";
	$output="$html_dir/$name.html";
	$command="cd $pod_base; pod2html   <$input > $output $pod_path1 $html_base1";
  print STDERR "$command \n";
	system($command) == 0
         or  print STDERR "NAME=$input \n";
}





}

print STDERR " $html_dir \n"; 
system("cd $html_dir_base;ln -s seplib.html index.html");
$dir=$html_dir_base;&create_make;

}






