# simple makefile for a Rappture-based program

CC		= g++
CFLAGS		= -g -Wall -O3

RAPPTURE_DIR	= /apps/share/rappture/20121012
INCLUDES	= -I$(RAPPTURE_DIR)/include -I/usr/include
LIBS		= -L$(RAPPTURE_DIR)/lib -lrappture -lm -L/usr/lib -lgsl -lgslcblas -lm


all: mainc

mainc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $< -o ../bin/$@ $(LIBS)

clean:
	$(RM) mainc run*.xml
