# library of routines, useful for books

# sets directory list
sub getDirs {
    my (@indirs);
    chomp  (my $papers = `gmake PAPERS.varvalue`);
    if ($papers) {
	@indirs = split (/\s/,$papers);
    } else {
	opendir (DIR, '..');
	@indirs = grep (/^[a-z]/,readdir (DIR));
	closedir (DIR);
    }
    foreach $dir (@indirs) {
	my $fuldir = '../' . $dir;
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
    my ($book, $logfile);
    chomp (local ($_) = `gmake LATOPTS.varvalue`);
    while (/\boption=(\S+)/g) {
	$book = grep (/^book$/,split (/\,/,$1));
    }
    chomp ($logfile  = `gmake LOGFILE.varvalue`);
    ($book, $logfile);
}

1; # should be the last line
