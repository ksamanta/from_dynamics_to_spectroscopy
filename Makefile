CC=g++
EXE=driver.x
CCFLAGS=-g -O3 -Wall -Wextra -pedantic -mtune=opteron -m64
SOURCES=driver.cpp spinbosoninput.cpp spinboson.cpp
GSL=/data/home/kousik/GSL

$(EXE) : $(SOURCES)
	$(CC) $(CCFLAGS) -o $(EXE) $(SOURCES)  -L$(GSL)/lib -I$(GSL)/include -lgsl -lgslcblas

clean:
	rm -f $(EXE)

all: clean $(EXE)
	
