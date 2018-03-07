#--------------------------------------------------------------------

SUBDIRS = holly_subfolder

#--------------------------------------------------------------------

# Extra files (if you want to pack them up)
EXTRA_DIST=

#--------------------------------------------------------------------

# Name of executables
# N.B.: Essentially just the filename without the .cc
check_PROGRAMS= two_d_act_inhibit one_d_act_inhibit test

# Sources for first executable
two_d_act_inhibit_SOURCES = two_d_act_inhibit.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_act_inhibit_LDADD = -L@libdir@ \
	                     -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for first executable
one_d_act_inhibit_SOURCES = one_d_act_inhibit.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_act_inhibit_LDADD = -L@libdir@  \
	                     -lgeneric $(EXTERNAL_LIBS) $(FLIBS) 

# Sources for first executable
test_SOURCES = test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
test_LDADD = -L@libdir@  \
	                     -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------------------------------------

# Essentially a copy and paste job...
AM_CPPFLAGS += -I@includedir@

