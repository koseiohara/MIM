EXEC=MIM

#VERSION=dev
VERSION=0.33r1#   ### rewrite lagmain.f90

INSTALL=~/bin


##### Fujitsu f90 #####
#FC=f90
#OPTION=-Am -X9 -Kfast -Et#                  # optimized
#OPTION=-Am -X9 -O0#                         # no optimized, safe
#OPTION=-Am -X9 -Et -O0 -Haesux#             # for debug


##### Intel ifort (for 10.0) #####
#OPTION=-warn all -O0 -assume byterecl -convert little_endian -traceback -check all -debug-parameters all
#OPTION=-warn all -O1 -assume byterecl -convert little_endian -traceback 
#OPTION=-warn all -O1 -assume byterecl -convert little_endian -traceback
#OPTION=-assume byterecl -O3 -fast -warn all -traceback -heap-arrays    # optimized
#OPTION=-assume byterecl -O2 -warn all -traceback -heap-arrays    # optimized
#OPTION=-assume byterecl -O0 -warn all -traceback -g -heap-arrays#    # no optimized
#OPTION=-assume byterecl -O0 -warn all -g -traceback -heap-arrays#    # no optimized
#OPTION=-assume byterecl -C -warn all#     # for debug

##### gfortran
FC=gfortran
FC=ifort

ifeq (${FC}, ifort)
	OPTION=-warn all -O1 -assume byterecl -convert little_endian -traceback 
else
	OPTION=-O1 -Wall -fbacktrace -fconvert=little-endian
endif


include Makefile.inc


%.o : %.f90
	${FC} -c $< ${OPTION}


all : ${EXEC}

${EXEC} : ${OBJS}
	${FC} -o $@ $^


.PHONY : run clean re

run :
	./${EXEC}

clean:
	@ rm -vf *.o *.mod ${EXEC}

re : clean all



