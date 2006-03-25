#
# $Id: Makefile,v 1.5 2006/03/25 21:52:34 trs137 Exp $
#
CC= g++
LD= g++
RM= /bin/rm -f
CFLAGS= -g -Wall -pedantic
LDFLAGS= -Wall
SRCS = Lambert.cpp Vector.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = orbgnosis

$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(TARGET) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
