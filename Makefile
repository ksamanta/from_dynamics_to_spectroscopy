#==============================================================
# A simple GNU makefile to compile the SpinBoson code
#
# Author: Kousik Samanta (2012-11-02)
#==============================================================

CPP=g++
EXE=driver.x
OPT=-O3 -mtune=opteron -m64
WARN=-Wall -Wextra -pedantic 
DEBUG=-v -time -g -H -Q 
SOURCE=driver.cpp spinbosoninput.cpp spinboson.cpp
GSL=/data/home/kousik/GSL

$(EXE): $(SOURCE)
	$(CPP) -o $(EXE) $(SOURCE)  \
		-L$(GSL)/lib -I$(GSL)/include -lgsl -lgslcblas \
		$(OPT) $(WARN) 
	@echo " "
	@echo "---------------------------------------------"
	@echo "Compilation successful. Yay!!"
	@date

debug: 
	$(CPP) -o $(EXE) $(SOURCE)  \
		-L$(GSL)/lib -I$(GSL)/include -lgsl -lgslcblas \
		$(OPT) $(WARN) $(DEBUG)
	@echo " "
	@echo "------------------------------------"
	@echo "End of compilation for debugging."
	@date

clean:
	rm -f $(EXE)

all: clean $(EXE)
	
