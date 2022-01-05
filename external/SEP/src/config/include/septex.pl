# library of routines, useful for tex documents

package comment;

sub removeComments {
    local ($_) = @_;
    if (s/([^\\])?(\%.*)$/$1/o) {
	$comment = $2;
    } else {
	$comment = '';
    }
    @verbs = split (/\\verb/); 
    my $firstverb = shift (@verbs);
    foreach $verb (@verbs) {
	$verb =~ s/^(.)//o;
	my $patt = $1;
	push (@verbchar,$patt);
	$patt = "\\" . $patt if ($patt =~ /\W/);
	$verb =~ s/$patt/\%/o unless ($verb =~ s/$patt/\%/o);
	$verb =~ s/^((\\%|[^\%])+[^\\])\%//o;
	push (@verb,$1);
   }
    $firstverb . join ('',@verbs);
}

sub putComments {
    my $line = shift;
    $line =~ s/\%/&Verb/eg;          # put verb back
    $line . $comment;                 #put comments back
}

# Verb stuff doesn't seem to work any more - needs deep debugging
sub Verb {
    my $verbchar = shift (@verbchar);
    join ('','\verb',$verbchar,shift (@verb),$verbchar);
}

package main;


# processes included files (recursively)
sub checkInputs {
    my (@lines) = @_;
    my $indir;
    foreach $line (@lines) { 
	$line = &comment::removeComments ($line);
	while ($line =~ /\\inputdir\{([^\}]+)\}/g) {
	    $indir = $1;
	}
	while ($line =~ /\\(sep)?(input|include)\s*\{([^\}]+)_act/g) {
            if (defined (%includes)) {
                 next unless ($includes{$3});
            }
	    my $inputfile = $3;
	    chdir ($indir) if ($indir && $1);
	    my $inputtex = $inputfile . "\.tex"; # provide suffix
	    $inputfile = $inputtex if (-e $inputtex);
	    &search_n_replace ($inputfile);
	}
    }
}

# sets directory list
sub getDirs {
    use Cwd;
    my (@indirs);
    chomp  (my $papers = `gmake PAPERS.varvalue`);
    if ($papers) {
	@indirs = split (/\s/,$papers);
    } else {
	opendir (DIR, '..');
	@indirs = grep (/^[a-z]/,readdir (DIR));
	closedir (DIR);
    }
    my $pwd = getcwd ();
    $pwd =~ s/[^\/]+$//o;
    foreach $dir (@indirs) {
	my $fuldir = ($dir =~ /^[\/\.]/)? $dir: $pwd . $dir;
	my $texfile = $fuldir . '/paper.tex';
	if ((-d $fuldir) && (-e $texfile)) {
	    push (@dirs,$dir);
	    push (@fuldirs,$fuldir);
	}
    }
}


# handler for interruption
sub handler {
    print STDERR "\n$0: aborting...\n";
    kill ('INT', $child_pid) if ($child_pid);
    exit(-1);
}

# forks to control interruption
sub syswait {
    local($_) = @_;
     if ($child_pid = fork) {
	my $status = waitpid($child_pid, 0);
	$child_pid = 0;
	return($?);
    } else {
	exec($_);
	print STDERR "@_[0]:  $!\n";
	exit($!);
    }
}

sub getBook {
    my $option = shift;
    $option = join ($option,'^','$');
    my ($book, $logfile);
    chomp (local ($_) = `gmake LATOPTS.varvalue`);
    while (/\boption=(\S+)/g) {
	$book = grep (/$option/,split (/\,/,$1));
    }
    chomp ($logfile  = `gmake LOGFILE.varvalue`);
    ($book, $logfile);
}


1; # should be the last line

