
#CC=/Share/home/zhangqf/lipan/gcc-6.1.0/bin/bin/g++
CC=g++

#export C_INCLUDE_PATH=/Share/home/zhangqf/lipan/gccDependance/gmp/include:/Share/home/zhangqf/lipan/gccDependance/mpc/include:/Share/home/zhangqf/lipan/gccDependance/mpfr/include:$C_INCLUDE_PATH
#export CPLUS_INCLUDE_PATH=/Share/home/zhangqf/lipan/gccDependance/gmp/include:/Share/home/zhangqf/lipan/gccDependance/mpc/include:/Share/home/zhangqf/lipan/gccDependance/mpfr/include:$CPLUS_INCLUDE_PATH
#export LD_LIBRARY_PATH=/Share/home/zhangqf/lipan/gccDependance/gmp/lib:/Share/home/zhangqf/lipan/gccDependance/mpc/lib:/Share/home/zhangqf/lipan/gccDependance/mpfr/lib:/Share/home/zhangqf/lipan/gcc-6.1.0/bin/lib64:$LD_LIBRARY_PATH

CFLAGS=-std=c++11 -Wall
ROOT=c++/GenDuplexGroup/GenDuplexGroup/

all: paris

paris: TimeFunc.o Files.o environment.o DuplexGroup.o main.o removeCR
	cd $(ROOT) && \
	$(CC) -I. main.o TimeFunc.o Files.o environment.o DuplexGroup.o -o paris && \
	rm *.o
	cd $(ROOT) && \
	mv paris ../../../bin/ && \
	mv removeCR ../../../bin/

test:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && echo "hello" `pwd`

main.o:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -c main.cpp

TimeFunc.o:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -c TimeFunc.cpp

Files.o:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -c Files.cpp

environment.o:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -c environment.cpp

DuplexGroup.o:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -c DuplexGroup.cpp

removeCR:
	cd c++/GenDuplexGroup/GenDuplexGroup/ && \
	$(CC) $(CFLAGS) -o removeCR removeCR.cpp

