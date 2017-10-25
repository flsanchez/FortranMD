.PHONY: default help object executable all clean
CC = gfortran

LD = $(CC)

SOURCE_F90 = $(wildcard *.f90)
OBJECTS_F90 = $(patsubst %.f90, ../bin/%_f90.o, $(SOURCE_F90))


F90FLAGS = -O3 -Wall -Wno-tabs
EXECUTABLE = ../bin/md.e

default: help

help:
	@echo "\
Options:\n\n\
  make objects:       compiler makes objects for every *.c\n\
  make executable:    compiler makes executable\n\
  make all:           build all previous\n\
  make clean:         delete output files\n\
  make help:          display this help"

objects: $(OBJECTS_F90)

executable: $(EXECUTABLE)

all: objects executable

../bin/%_f90.o: %.f90
	$(CC) $(F90FLAGS) -c $^ -o ../bin/$@

%.e: $(OBJECTS_F90)
	$(LD) $^ -o $@

clean:
	rm -rfv $(OBJECTS_F90) $(EXECUTABLE)
