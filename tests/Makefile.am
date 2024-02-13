

zakharov_srcdir = $(top_srcdir)/src
cpptools_srcdir = $(zakharov_srcdir)/c++-tools

zakharov_libdir = $(top_srcdir)/lib

AM_CPPFLAGS = -I$(zakharov_srcdir) -I$(cpptools_srcdir)

# Executable targets
check_PROGRAMS = Fftw
bin_PROGRAMS = Zakharov

if COMPILER_IS_COMPAQ_CXX
templates_zakharov_dir = -ptr $(zakharov_libdir)/cxx_repository
endif

Fftw_SOURCES = Fftw.cpp
Fftw_LDADD = $(zakharov_libdir)/libzakharov.la

Zakharov_SOURCES = Zakharov.cpp
Zakharov_LDFLAGS = $(templates_zakharov_dir)
Zakharov_LDADD = $(zakharov_libdir)/libzakharov.la

## clean up template repositories
clean-local:
	-rm -rf cxx_repository ti_files ii_files

## Check 
check:
	@echo Running test suite...
	@for prog in $(check_PROGRAMS); do \
		echo $$prog && (./$$prog > /dev/null 2>&1) || exit 1; \
		done
	@echo All tests passed.

