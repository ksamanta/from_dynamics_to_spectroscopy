CC=g++
EXE=driver.x
CCFLAGS=-g -O3 -Wall -Wextra -pedantic -mtune=opteron -m64
SOURCES=driver.cpp spinbosoninput.cpp spinboson.cpp
GSL_LIB=/data/home/kousik/GSL/lib
GSL_INC=/data/home/kousik/GSL/include

$(EXE) : $(SOURCES)
	$(CC) $(CCFLAGS) -o $(EXE) $(SOURCES)  -L$(GSL_LIB) -I$(GSL_INC) -lgsl -lgslcblas

