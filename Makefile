#Set directories for input files/source code (ID),
# output object files (OD), and executables (XD)
ID = ./src
OD = ./obj
XD =.

CXX=g++
#CXX=clang++

OMP=-fopenmp

OPT=-O3
WARN=-Wall -Wextra -Wpedantic -Wdouble-promotion -Wconversion

ifeq ($(CXX),clang++)
  WARN += -Wno-sign-conversion
endif

CXXFLAGS=-std=c++11 $(WARN) $(OMP) $(OPT)

detected_OS := $(shell uname -s) #will return the Operating system name
$(info )
$(info Detected operating system: $(detected_OS))
$(info )

LIBS= #-lgsl -lgslcblas #-lm

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

################################################################################

#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 transientDM combinedX \
)

#default tagret
all: checkObj checkXdir $(ALLEXES)

################################################################################
## Dependencies:

#All programs depend on these generic common headers:
COMH = $(addprefix $(ID)/, \
 ClockInfo.h /DataIO.h \
)

# Rule for files that have .cpp AND a .h file
# They depend 'only' on their own header, + generic common headers
$(OD)/%.o: $(ID)/%.cpp $(ID)/%.h $(COMH)
	$(COMP)

# Rule for files that _don't_ have a .h header. (mains)
# These also depend on the common headers
$(OD)/%.o: $(ID)/%.cpp $(COMH)
	$(COMP)

# Here: List rules for any other progs that don't fit above rules?
$(OD)/main.o: $(ID)/main.cpp $(COMH) $(ID)/ClockNetwork.h \
$(ID)/ChronoTimer.h
	$(COMP)

$(OD)/transientDM.o: $(ID)/transientDM.cpp $(COMH) $(ID)/ClockNetwork.h \
$(ID)/ChronoTimer.h
	$(COMP)

################################################################################
# Link + build all final programs

$(XD)/transientDM: $(OD)/transientDM.o $(OD)/ClockNetwork.o \
$(OD)/DMs_signalTemplates.o $(OD)/RNG_randomNumberGenerators.o
	$(LINK)

$(XD)/combinedX: $(OD)/combinedX.o
	$(LINK)

# $(XD)/main: $(OD)/main.o $(OD)/ClockNetwork.o $(OD)/DataIO.o \
# $(OD)/DMs_signalTemplates.o $(OD)/ChronoTimer.o
# 	$(LINK)

################################################################################

checkObj:
	@if [ ! -d $(OD) ]; then \
	  echo '\n ERROR: Directory: '$(OD)' doesnt exist - please create it!\n'; \
	  false; \
	else \
	  echo 'OK'; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

.PHONY: clean checkObj checkXdir
clean:
	rm -f $(ALLEXES) $(OD)/*.o
