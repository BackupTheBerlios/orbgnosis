#
# $Id: Makefile,v 1.8 2006/04/04 23:04:33 trs137 Exp $
#
CC= g++
LD= g++
RM= /bin/rm -f
#CFLAGS= -g -Wall -pedantic
CFLAGS=-O2
#LDFLAGS= -Wall
LDFLAGS=-O2
SRCS = Lambert.cpp Vector.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = orbgnosis

$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(TARGET) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
