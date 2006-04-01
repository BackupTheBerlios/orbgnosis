#
# $Id: Makefile,v 1.7 2006/04/01 02:55:38 trs137 Exp $
#
CC= g++
LD= g++
RM= /bin/rm -f
CFLAGS= -g -Wall -pedantic
#CFLAGS=-O2
LDFLAGS= -Wall
#LDFLAGS=-O2
SRCS = Lambert.cpp Vector.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = orbgnosis

$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(TARGET) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
