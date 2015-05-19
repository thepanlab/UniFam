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
SRCS := utils.cpp directoryStructure.cpp ReadFile.cpp

# use macro to get object files
OBJS = $(SRCS:.cpp=.o)

# program names
SPLIT = split_group
PIPE = unifam_build

.PHONY: default split_group unifam_build all clean

default: all
split: split_group
pipe: unifam_build

all: $(SPLIT) $(PIPE)
	@echo Done
$(SPLIT): $(OBJS) split_group.o
	$(CC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(SPLIT) $(OBJS) $(LIBS) split_group.o
	@echo $(SPLIT) compiled

$(PIPE): $(OBJS) unifam_build.o
	$(MPICC) $(OMPFLAG) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(PIPE) $(OBJS) $(LIBS) unifam_build.o
	@echo $(PIPE) compiled

split_group.o: split_group.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

unifam_build.o: unifam_build.cpp
	$(MPICC) $(OMPFLAG) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Either way will work, first one uses old-Fashioned suffix rules; second one uses static pattern rules
#.cpp.o:
$(OBJS): %.o:%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(OBJS) $(PIPE) $(SPLIT) unifam_build.o split_group.o
