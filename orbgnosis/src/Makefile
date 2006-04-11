#
# $Id: Makefile,v 1.11 2006/04/11 11:54:20 trs137 Exp $
#
CC= g++
LD= g++
RM= /bin/rm -f
#CFLAGS= -g -Wall -pedantic
CFLAGS=-O2
#LDFLAGS= -Wall
LDFLAGS=-O2
SRCS = Lambert.cpp Vector.cpp Orbgnosis.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = orb

$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(TARGET) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
