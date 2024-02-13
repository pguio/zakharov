
zakharov_srcdir = $(top_srcdir)/src
cpptools_srcdir = $(top_srcdir)/src/c++-tools

zakharov_libdir = $(top_srcdir)/lib

AM_CPPFLAGS = -I$(zakharov_srcdir) -I$(cpptools_srcdir)

zakharov_src = \
$(zakharov_srcdir)/SIconsts.cpp \
$(zakharov_srcdir)/default-mode.cpp \
$(zakharov_srcdir)/density-fluct.cpp \
$(zakharov_srcdir)/density-parab.cpp \
$(zakharov_srcdir)/diagnostics.cpp \
$(zakharov_srcdir)/efluct.cpp \
$(zakharov_srcdir)/egauss.cpp \
$(zakharov_srcdir)/mapper.cpp \
$(zakharov_srcdir)/solitons.cpp \
$(zakharov_srcdir)/spectra.cpp \
$(zakharov_srcdir)/thermal-fluct.cpp \
$(zakharov_srcdir)/zakharov.cpp \
$(zakharov_srcdir)/time-integrate-classic.cpp \
$(zakharov_srcdir)/time-integrate-pump.cpp \
$(zakharov_srcdir)/time-integrate-fixed-n.cpp \
$(zakharov_srcdir)/time-integrate-fixed-E.cpp \
$(zakharov_srcdir)/time-integrate-linear.cpp \
$(zakharov_srcdir)/time-integrate-noFp.cpp 

if GNUPLOT_ENABLED
gnuplot_interface = $(cpptools_srcdir)/gnuplot-interface.cpp
endif

cpptools_src = \
$(cpptools_srcdir)/hdf-interface.cpp \
$(cpptools_srcdir)/parser.cpp \
$(gnuplot_interface)

lib_LTLIBRARIES = libzakharov.la

if COMPILER_IS_COMPAQ_CXX
templates_zakharov_dir = -ptr $(zakharov_libdir)/cxx_repository
endif

libzakharov_la_SOURCES = $(zakharov_src) $(cpptools_src)
libzakharov_la_CXXFLAGS = $(templates_zakharov_dir)
libzakharov_la_LDFLAGS = -version-info 0:0:0

## Clean up template repositories
clean-local:
	-rm -rf cxx_repository ti_files ii_files

