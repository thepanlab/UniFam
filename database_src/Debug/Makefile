-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk


#  -c : not to run the linker
#  -W -Wall -pedantic : maximal sensitivity of compiler against possible problems
#  -Wno-variadic-macros : do not warn variadic macros usage (from C90 was turned on by default)
#  -fmessage-length=0 : no line wrapping for error messages, each error message appears on a single line
#  -MP : add a phony target for each dependency other than the main file, causing each to depend on nothing
#  -MMD : dependency and header file options

# compiler
CC = g++
MPICC = mpic++

ifneq ($(MAKECMSGOALS), clean)
	include $(sources:.cpp=.d)
endif
# compile time flags
CFLAGS = -O0 -g3 -W -Wno-long-long -pedantic -Wno-variadic-macros -fmessage-length=0 -MMD -MP -std=c++11
CFLAGS+=-Wall
OMPFLAG= -fopenmp


# source files and object files
# directive for source files
vpath %.cpp ../
vpath %.h ../
#SRCS := $(wildcard ../*.cpp)
SRCS := utils.cpp directoryStructure.cpp runBioProg.cpp pfclust.cpp ReadFile.cpp

# use macro to get object files
OBJS = $(SRCS:.cpp=.o)

# program names
SPLIT = split_fasta
PIPE = pipeline
TEST = testdrive

.PHONY: default split_fasta pipeline all clean

default: pipeline

split: split_fasta

pipe: pipeline

test: testdrive

all: $(SPLIT) $(PIPE)
	@echo Done
$(SPLIT): $(OBJS) split_fasta.o
	$(MPICC) $(OMPFLAG) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(SPLIT) $(OBJS) $(LIBS) split_fasta.o
	@echo $(SPLIT) compiled

$(PIPE): $(OBJS) till_hmmbuild.o
	$(MPICC) $(OMPFLAG) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(PIPE) $(OBJS) $(LIBS) till_hmmbuild.o
	@echo $(PIPE) compiled

$(TEST): $(OBJS) testdrive.o
	$(MPICC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(TEST) $(OBJS) $(LIBS) testdrive.o
	@echo $(TEST) compiled

split_fasta.o: split_fasta.cpp
	$(MPICC) $(OMPFLAG) $(CFLAGS) $(INCLUDES) -c $< -o $@

till_hmmbuild.o: till_hmmbuild.cpp
	$(MPICC) $(OMPFLAG) $(CFLAGS) $(INCLUDES) -c $< -o $@

testdrive.o: testdrive.cpp
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $< -o $@
# Either way will work, first one uses old-Fashioned suffix rules; second one uses static pattern rules
#.cpp.o:
$(OBJS): %.o:%.cpp
	$(MPICC) $(OMPFLAG) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(OBJS) $(PIPE) $(SPLIT)
