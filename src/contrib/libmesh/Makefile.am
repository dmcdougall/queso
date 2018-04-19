#bindir=@bindir@/@host@_@METHOD@
#libdir=@libdir@/@host@_@METHOD@

SUBDIRS = include #src

# AUTOMAKE_OPTIONS = subdir-objects
# ACLOCAL_AMFLAGS  = -I m4
# AM_CFLAGS        = $(libmesh_CFLAGS)
# AM_LDFLAGS       = $(libmesh_LDFLAGS)

# AM_CPPFLAGS are method-independent cppflags that
# we use when compiling libmesh proper, or its utility
# codes.  These will be augmented with method-specific
# cppflags later on
# $(libmesh_contrib_INCLUDES)


# Put EIGEN here?
AM_CPPFLAGS      = -DLIBMESH_IS_COMPILING_ITSELF \
									 -I$(top_builddir)/src/contrib/libmesh/include \
                   $(libmesh_optional_INCLUDES) \
									 $(EIGEN_INCLUDE) \
                   -I$(top_builddir)/include # required for libmesh_version.h

LIBS             = $(libmesh_optional_LIBS)

# # additional files which must be included in 'make dist'
# EXTRA_DIST = reference_elements \
#              include/include_HEADERS \
#              include/rebuild_include_HEADERS.sh \
#              src/libmesh_SOURCES \
#              src/rebuild_libmesh_SOURCES.sh \
#              contrib/tecplot/README \
#              contrib/bin/strip_dup_incl_paths.pl \
#              contrib/bin/strip_dup_libs.pl \
#              CITATION
#
# nobase_data_DATA = \
#         reference_elements/2D/one_quad.xda \
#         reference_elements/2D/one_quad8.xda \
#         reference_elements/2D/one_quad9.xda \
#         reference_elements/2D/one_tri.xda \
#         reference_elements/2D/one_tri6.xda \
#         reference_elements/3D/one_hex.xda \
#         reference_elements/3D/one_hex20.xda \
#         reference_elements/3D/one_hex27.xda \
#         reference_elements/3D/one_prism.xda \
#         reference_elements/3D/one_prism15.xda \
#         reference_elements/3D/one_prism18.xda \
#         reference_elements/3D/one_pyramid.xda \
#         reference_elements/3D/one_pyramid13.xda \
#         reference_elements/3D/one_pyramid14.xda \
#         reference_elements/3D/one_tet.xda \
#         reference_elements/3D/one_tet10.xda

# Tools in the auxillary directory
# AUX_DIST = build-aux/install-sh \
#            build-aux/missing

# Support for pkgconfig
# pkgconfigdir     = $(libdir)/pkgconfig
# pkgconfig_DATA   = contrib/utils/Make.common \
#                    contrib/utils/libmesh.pc \
#                    contrib/utils/libmesh-opt.pc \
#                    contrib/utils/libmesh-dbg.pc \
#                    contrib/utils/libmesh-devel.pc \
#                    contrib/utils/libmesh-prof.pc \
#                    contrib/utils/libmesh-oprof.pc
#
# configdir        = $(sysconfdir)/libmesh
# config_DATA      = contrib/utils/Make.common \
#                    contrib/utils/libmesh.pc \
#                    contrib/utils/libmesh-opt.pc \
#                    contrib/utils/libmesh-dbg.pc \
#                    contrib/utils/libmesh-devel.pc \
#                    contrib/utils/libmesh-prof.pc \
#                    contrib/utils/libmesh-oprof.pc

# install the old make_dependencies.pl script for backwards compatibility
# contribbindir      = $(prefix)/contrib/bin
# contribbin_SCRIPTS = contrib/bin/make_dependencies.pl libtool
# EXTRA_DIST        += contrib/bin/make_dependencies.pl

# Support for libmesh-config in $(exec_prefix)/bin
# libmesh_configdir      = $(exec_prefix)/bin
# libmesh_config_SCRIPTS = contrib/bin/libmesh-config

# This rule extends the built-in make install rule, and can
# be used to:
#  when doing 'make install prefix=/other/path', update generated files
#  when doing 'make install', symlik etc/libmesh/Make.common
#  when doing 'make install', update libmesh_config
#  fix up buggy libtool postdeps under specific compiler/mpi configurations
# install-data-hook:
# 	@if (test "x$(prefix)" != "x@prefix@") ; then \
#           oldprefix="@prefix@" ; \
#           newprefix="$(prefix)" ; \
#           for genfile in `echo $(DESTDIR)$(prefix)/lib/pkgconfig/*.pc $(DESTDIR)$(prefix)/etc/libmesh/*.pc $(DESTDIR)$(prefix)/etc/libmesh/Make.common $(DESTDIR)$(prefix)/lib/pkgconfig/Make.common $(DESTDIR)$(prefix)/bin/libmesh-config` ; do \
#             echo " " ; \
#             echo " *** replacing $$oldprefix" ; \
#             echo " ***      with $$newprefix" ; \
#             echo " *** in generated file $$genfile" ; \
#             echo " " ; \
#             cd $(DESTDIR)$(prefix) && cat $$genfile | $(SED) "s,$$oldprefix,$$newprefix,g" > $$genfile.replaced ; \
#             cd $(DESTDIR)$(prefix) && mv $$genfile.replaced $$genfile ; \
#           done ; \
#         fi
# 	cd $(DESTDIR)$(prefix) && rm -f Make.common && $(LN_S) -f etc/libmesh/Make.common Make.common
# 	cat $(DESTDIR)$(libmesh_configdir)/libmesh-config | $(SED) "s/has_been_installed=no/has_been_installed=yes/g" > $(DESTDIR)$(libmesh_configdir)/libmesh-config.installed
# 	mv $(DESTDIR)$(libmesh_configdir)/libmesh-config.installed $(DESTDIR)$(libmesh_configdir)/libmesh-config && chmod +x $(DESTDIR)$(libmesh_configdir)/libmesh-config
# 	cat $(DESTDIR)$(prefix)/contrib/bin/libtool | $(SED) "s/-l //g" > $(DESTDIR)$(prefix)/contrib/bin/libtool.fixed
# 	mv $(DESTDIR)$(prefix)/contrib/bin/libtool.fixed $(DESTDIR)$(prefix)/contrib/bin/libtool && chmod +x $(DESTDIR)$(prefix)/contrib/bin/libtool
#
#
# # this rule invokes emacs on each source file to remove trailing whitespace.
# beautify:
# 	for file in `find $(srcdir)/{src,include,examples,tests} -name "*.C" -o -name "*.h" -o -name Makefile.am -type f` ; do \
#           echo $$file ; \
#           emacs -batch $$file --eval '(delete-trailing-whitespace)' -f save-buffer 2>/dev/null ; \
#         done
#
#
# # Additional files to be deleted by 'make distclean'
# # _configs.sed seems to be leftover by prefix_config
# DISTCLEANFILES =  \
#            _configs.sed \
# 	   a.out  # intel-14.0 seems to leave an a.out from configure with mpich, clean that up for distclean

# # Additional files to be deleted by 'make maintainer-clean'
# MAINTAINERCLEANFILES =  \
#          $(AUX_DIST)                              \
#          Makefile.in                              \
#          aclocal.m4                               \
#          aminclude.am                             \
#          autom4te.cache/*                         \
#          config.h                                 \
#          config.h.in                              \
#          config.log                               \
#          config.status                            \
#          config.sub                               \
#          configure                                \
#          data/Makefile.in                         \
#          src/Makefile.in                          \
#          src/physics/properties/Makefile.in       \
#          test/Makefile.in                         \
#          test/common/Makefile.in                  \
#          test/comp_ns/Makefile.in                 \
#          test/unit/Makefile.in



# SUBDIRS += contrib
noinst_LTLIBRARIES = libnv.la
include src/libmesh_SOURCES

# A convenience library to hold proper libMesh
# objects.  This will get appended with any contributed
# sources to create the final library.
# lib_LTLIBRARIES =

# if LIBMESH_DBG_MODE
  # lib_LTLIBRARIES        += libmesh_dbg.la
  # libmesh_dbg_la_SOURCES  = $(libmesh_SOURCES)
  # libmesh_dbg_la_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  # libmesh_dbg_la_CXXFLAGS = $(CXXFLAGS_DBG)
  # libmesh_dbg_la_CFLAGS   = $(CFLAGS_DBG)
  # libmesh_dbg_la_LIBADD   = contrib/libcontrib_dbg.la $(LIBS)
# endif

# if LIBMESH_DEVEL_MODE
#   lib_LTLIBRARIES          += libmesh_devel.la
#   libmesh_devel_la_SOURCES  = $(libmesh_SOURCES)
#   libmesh_devel_la_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
#   libmesh_devel_la_CXXFLAGS = $(CXXFLAGS_DEVEL)
#   libmesh_devel_la_CFLAGS   = $(CFLAGS_DEVEL)
#   libmesh_devel_la_LIBADD   = contrib/libcontrib_devel.la $(LIBS)
# endif
#
# if LIBMESH_OPT_MODE
#   lib_LTLIBRARIES        += libmesh_opt.la
#   libmesh_opt_la_SOURCES  = $(libmesh_SOURCES)
#   libmesh_opt_la_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
#   libmesh_opt_la_CXXFLAGS = $(CXXFLAGS_OPT)
#   libmesh_opt_la_CFLAGS   = $(CFLAGS_OPT)
#   libmesh_opt_la_LIBADD   = contrib/libcontrib_opt.la $(LIBS)
# endif

# if LIBMESH_PROF_MODE
#   lib_LTLIBRARIES          += libmesh_prof.la
#   libmesh_prof_la_SOURCES  = $(libmesh_SOURCES)
#   libmesh_prof_la_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
#   libmesh_prof_la_CXXFLAGS = $(CXXFLAGS_PROF)
#   libmesh_prof_la_CFLAGS   = $(CFLAGS_PROF)
#   libmesh_prof_la_LIBADD   = contrib/libcontrib_prof.la $(LIBS)
# endif
#
# if LIBMESH_OPROF_MODE
#   lib_LTLIBRARIES          += libmesh_oprof.la
#   libmesh_oprof_la_SOURCES  = $(libmesh_SOURCES)
#   libmesh_oprof_la_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
#   libmesh_oprof_la_CXXFLAGS = $(CXXFLAGS_OPROF)
#   libmesh_oprof_la_CFLAGS   = $(CFLAGS_OPROF)
#   libmesh_oprof_la_LIBADD   = contrib/libcontrib_oprof.la $(LIBS)
# endif


# if LIBMESH_ENABLE_CPPUNIT
#  SUBDIRS += tests
# endif

###########################################################


# bin_PROGRAMS = # empty, append below
# bin_SCRIPTS  = # empty, append below
# CLEANFILES   = # empty, append below
#
# ###########################################################
# # Utility programs
# opt_programs   = # empty, append below
# devel_programs = # empty, append below
# dbg_programs   = # empty, append below
# prof_programs  = # empty, append below
# oprof_programs = # empty, append below

# output_libmesh_version
# opt_programs                         += output_libmesh_version-opt
# output_libmesh_version_opt_SOURCES    = src/apps/output_libmesh_version.C
# output_libmesh_version_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
# output_libmesh_version_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
# output_libmesh_version_opt_LDADD      = libmesh_opt.la
#
# devel_programs                       += output_libmesh_version-devel
# output_libmesh_version_devel_SOURCES  = src/apps/output_libmesh_version.C
# output_libmesh_version_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
# output_libmesh_version_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
# output_libmesh_version_devel_LDADD    = libmesh_devel.la

# dbg_programs                         += output_libmesh_version-dbg
# output_libmesh_version_dbg_SOURCES    = src/apps/output_libmesh_version.C
# output_libmesh_version_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
# output_libmesh_version_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
# output_libmesh_version_dbg_LDADD      = libmesh_dbg.la

# if LIBMESH_OPT_MODE
#   bin_PROGRAMS += $(opt_programs)
# endif
#
# if LIBMESH_DEVEL_MODE
#   bin_PROGRAMS += $(devel_programs)
# endif

# if LIBMESH_DBG_MODE
  # bin_PROGRAMS += $(dbg_programs)
# endif

###########################################################
# Examples
# if LIBMESH_ENABLE_EXAMPLES
  # SUBDIRS += examples
# endif

###########################################################
# Documentation
# .PHONY: examples_doc doc

# SUBDIRS += doc

# #
# # documentation for the examples
# examples_doc:
# 	@cd $(top_builddir)/examples && $(MAKE) doc

# #
# # doxygen documentation - be sure to generate symlinks in include/libmesh!
# doc html:
# #	$(MAKE) examples_doc
# 	$(MAKE) -C doc/html doc

# # Upload the web page to github.io.
# upload:
# 	-rm -rf libMesh.github.io && git clone https://github.com/libMesh/libMesh.github.io.git
# 	cd libMesh.github.io && git rm -rf doxygen
# 	rsync -aLv --exclude=src --exclude=Makefile --exclude=.linkstamp $(top_builddir)/doc/html/ libMesh.github.io/
# 	gitcommit=`git log | head -1` && cd libMesh.github.io && git add * && git ci -m"updating web page from \"make upload\", libmesh.git at $$gitcommit"
# # uncomment below when ready
# #	cd libMesh.github.io && git push
# #	-rm -rf libMesh.github.io
#
#
# #
# # special things to do when running 'make dist'
# dist-hook:
# 	rm -rf `find $(distdir) -name libmesh_config.h`
#
#
# #
# # support 'make run_examples'
# run_examples:
# 	@$(MAKE) && cd $(top_builddir)/examples && $(MAKE) check
#
# #
# # support top-level 'make test_headers'
# test_headers:
# 	@$(MAKE) -C $(top_builddir)/include/libmesh
# 	@cd $(top_builddir)/include && $(MAKE) test_headers

# -------------------------------------------
# Optional support for code coverage analysis
# -------------------------------------------

# if CODE_COVERAGE_ENABLED
#
# lcov_dir=$(top_builddir)/doc/lcov/html
# lcov_bin=$(top_srcdir)/contrib/lcov/lcov-1.10/bin
#
# # General philosophy is to maintain code coverage for the
# # base library as generated by "make check" tests.
#
# lcov-report:
# 	@mkdir -p $(lcov_dir)
# 	$(lcov_bin)/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
# 	$(lcov_bin)/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep -v "`cd -P $(top_srcdir)/src && pwd`" | cut -d\| -f1 > $(lcov_dir)/remove
# 	$(lcov_bin)/lcov -q -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` > $(lcov_dir)/lcov.cleaned.info
# 	@rm $(lcov_dir)/remove
# 	@mv $(lcov_dir)/lcov.cleaned.info $(lcov_dir)/lcov.info
# 	$(lcov_bin)/genhtml -t "libMesh" -o $(lcov_dir) $(lcov_dir)/lcov.info
#
# lcov-reset:
# 	@rm -rf $(lcov_dir)
# 	@find . -name "*.gcda" -exec rm {} \;
# 	$(lcov_bin)/lcov --directory . --zerocounters
#
# coverage: lcov-reset check lcov-report
#
# CLEANFILES += src/apps/*.gcda src/apps/*.gcno
#
# endif
