# On Windows:
#
# The DOS executables in this distribution was produced with 
# Microsoft VisualStudio 2003.NET substituting its C++ compiler 
# with the Intel C++ and Fortran compilers. 
#
# geomag:
# Generate a Project of type "Visual C++, WIN32 console"
# using the option "empty project". Add existing file 
# geomag.c into the source directory, click on the project,
# convert the project to Intel C++, and build it. 
#
# magpoint, geomdr:
# Generate a Project of type "Intel Fortran console
# application" using the option "empty project". 
# Add existing files geomag.for and magpoint.for or
# geomdr.for to the source directory and build it.
#
#
# on Unix, compile desired program using
# 
#  make geomdr   ! produces grid
#  make magpoint ! produces point values
#  make geomag   ! C version of magpoint
#

DEBUG =  -g 	 			# -g = debug, -O = optimize  
ARCH =                  # add flags, if applicable

#gnu
#CC      = gcc ${ARCH} ${DEBUG} -c
#CCLINK  = gcc ${ARCH} ${DEBUG} -lm
#FOR     = g77 ${ARCH} ${DEBUG} -c
#FORLINK = g77 ${ARCH} ${DEBUG} -lm 

#Intel
CC =      /opt/intel_cc_80/bin/icc   ${ARCH} ${DEBUG} -wd1572 -c
CCLINK =  /opt/intel_cc_80/bin/icc   ${ARCH} ${DEBUG} -static
FOR =     /opt/intel_fc_80/bin/ifort ${ARCH} ${DEBUG} ${TYPEMAP} -c
FORLINK = /opt/intel_fc_80/bin/ifort ${ARCH} ${DEBUG} -static


# Subroutine

geomag.o: geomag.for Makefile
	${FOR} geomag.for -o geomag.o

# Main programs

geomdr: geomdr.for geomag.o Makefile 
	${FOR}  geomdr.for -o geomdr.o
	${FORLINK} geomdr.o geomag.o -o geomdr 

magpoint: magpoint.for geomag.o Makefile 
	${FOR}  magpoint.for -o magpoint.o
	${FORLINK} magpoint.o geomag.o -o magpoint 

geomag: geomag.c Makefile 
	${CC}  geomag.c -o geomag_C.o
	${CCLINK}  geomag_C.o -o geomag

