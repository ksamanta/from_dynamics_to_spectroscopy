#==============================================================
# A simple GNU makefile to compile the SpinBoson code
#
# Author: Kousik Samanta (2012-11-02)
#==============================================================

CPP=/data/home/kousik/local/bin/g++
EXE=driver.x
GSL=/data/home/kousik/GSL
SOURCES=driver.cpp spinbosoninput.cpp spinboson.cpp 
OMP_FLAG=-fopenmp
GSL_FLAG=-L${GSL}/lib -I${GSL}/include -lgsl -lgslcblas -lm 
OPT_FLAG=-O3 -mtune=opteron -m64 
WARN_FLAG=-Wall -Wextra -pedantic  
DEBUG_FLAG=-v -time -g -H -Q 


${EXE}: ${SOURCES}
	${CPP} -o ${EXE} ${SOURCES} \
	${OMP_FLAG}  \
	${OPT_FLAG}  \
	${WARN_FLAG} \
	${GSL_FLAG}

	@echo " "
	@echo "---------------------------------------------"
	@echo "Compilation successful. Yay!!"
	@date


debug: 
	${CPP} -o ${EXE} ${SOURCES} \
	${OMP_FLAG}  \
	${OPT_FLAG}  \
	${WARN_FLAG} \
	${GSL_FLAG}  \
	${DEBUG_FLAG}

	@echo " "
	@echo "------------------------------------"
	@echo "End of compilation for debugging."
	@date


clean:
	rm -f ${EXE}


all: clean ${EXE}
	

