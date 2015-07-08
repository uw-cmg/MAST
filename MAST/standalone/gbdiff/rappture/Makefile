# simple makefile for a Rappture-based program

CC		= g++
CFLAGS		= -g -Wall

RAPPTURE_DIR	= /apps/share/rappture/20121012
INCLUDES	= -I$(RAPPTURE_DIR)/include
LIBS		= -L$(RAPPTURE_DIR)/lib -lrappture -lm


all: mainc

mainc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $< -o ../bin/$@ $(LIBS)

clean:
	$(RM) mainc run*.xml
