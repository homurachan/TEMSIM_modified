#
# makefile to generate TEMSIM multislice package.
# assumes that FFTW already installed (see www.fftw.org)
#
# Put this file in the same directory as the TEMSIM
# C source files and type "make all" from a command line,
# to compile all of the programs.
#
#  type "make all" to compile everthing
#  type "make remove" to remove all of the compiled files.
#
# Each file name is assumed to be all lower case.
#
# convert to C++ 29-may-2012 ejk
# update for **cmd.cpp 7-may-2014 ejk
# last modified 7-may-2014 ejk
#

# define compiler with optimize flag
#CC = gcc -O
CC = g++ -O3
#DEL = del  # windows/mingw - doesn't work without .exe in file name
DEL = rm  # unix

# define libraries
MYLIBS = slicelib.o floatTIFF.o cfpix.o 
LIBS = ${MYLIBS}$
WLIBS = slicelib.o floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f

#
#  entry point to build everything
#
all:
	make atompot
	make autoslic
	make autostem
	make image
	make incostem
	make mulslice
	make probe
	make stemslic
	make sumpix


#
#  entry point to remove compiled files
#
remove:
	${DEL}$ atompot
	${DEL}$ autoslic
	${DEL}$ autostem
	${DEL}$ image
	${DEL}$ incostem
	${DEL}$ mulslice
	${DEL}$ probe
	${DEL}$ stemslic
	${DEL}$ sumpix
	${DEL}$ cfpix.o
	${DEL}$ slicelib.o
	${DEL}$ floatTIFF.o


#
#  main programs
#

atompot: atompot.cpp  ${MYLIBS}
	${CC} -o atompot atompot.cpp ${WLIBS}

autoslic: autoslic.cpp autosliccmd.cpp ${MYLIBS}
	${CC}  -o autoslic autosliccmd.cpp autoslic.cpp ${WLIBS}

autostem: autostem.cpp autostemcmd.cpp ${MYLIBS}
	${CC}  -o autostem autostemcmd.cpp autostem.cpp ${WLIBS} 

display: display.cpp  ${MYLIBS}
	${CC} -o display display.cpp ${WLIBS}

image: image.cpp  ${MYLIBS}
	${CC} -o image image.cpp ${WLIBS}

incostem: incostem.cpp incostemcmd.cpp ${MYLIBS}
	${CC}  -o incostem incostemcmd.cpp incostem.cpp probe.cpp ${WLIBS}

mulslice: mulslice.cpp ${MYLIBS}
	${CC} -o mulslice mulslice.cpp ${WLIBS}

probe: probe.cpp probecmd.cpp ${MYLIBS}
	${CC} -o probe probecmd.cpp probe.cpp ${WLIBS}

slicview: slicview.cpp ${MYLIBS}
	${CC} -o slicview slicview.cpp ${WLIBS}

stemslic: stemslic.cpp ${MYLIBS}
	${CC} -o stemslic stemslic.cpp ${WLIBS}

sumpix: sumpix.cpp ${MYLIBS}
	${CC} -o sumpix sumpix.cpp ${WLIBS}


#
# define subroutine library
#

cfpix.o: cfpix.cpp
	${CC} -c cfpix.cpp

slicelib.o: slicelib.cpp
	${CC} -c slicelib.cpp

floatTIFF.o: floatTIFF.cpp
	${CC} -c floatTIFF.cpp

