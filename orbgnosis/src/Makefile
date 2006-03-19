#
# $Id: Makefile,v 1.1 2006/03/19 23:58:43 trs137 Exp $
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

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(TARGET)
