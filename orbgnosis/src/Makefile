#
# $Id: Makefile,v 1.2 2006/03/20 03:52:22 trs137 Exp $
#
CC= /usr/local/bin/g++40
LD= /usr/local/bin/g++40
RM= /bin/rm -f
CFLAGS= -g -Wall -pedantic
LDFLAGS= -Wall
SRCS = Lambert.cpp Stumpff.cpp Vector.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = orbgnosis

$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(TARGET) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
