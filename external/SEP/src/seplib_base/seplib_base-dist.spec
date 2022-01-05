Summary: Geophysics software
Name: seplib-base-dist
Version: 6.0
Release: 1
Packager: Bob Clapp <bob@sep.stanford.edu>
Copyright: Stanford University
Group: Development/Libraries
URL: http://sepwww.stanford.edu/
Source: ftp://sepftp.stanford.edu/pub/sep-distr/seplib_base-6.0.tar.gz
Requires: gcc-g77
Requires: gcc
Requires: make >= 3.74
Requires: perl >= 5.001
Requires: gawk >= 3.0.3
Requires: tcsh
Requires: ld-linux.so.2
Requires:  libc.so.6
Requires: libm.so.6
Requires: ratfor90

%description
Geophysics visualization software. Useful for people studying what
goes on inside the earth.

%changelog
* Mon Nov 26 2001  Bob Clapp <bob@sep.stanford.edu>
- Created this RPM


%prep
%setup -n seplib-base-6.0

%build
csh LOCAL_BUILD
make all

%install
make install
%clean
make clean
%files
%attr(-,root,root) /usr/local/SEP
