#
# $Id: Makefile,v 1.13 2006/05/01 19:51:10 trs137 Exp $
#
# use gmake!
CC=icpc
LD=icpc
RM=/bin/rm -f
#CFLAGS=-Wall -Werror -ansi -pedantic -g
CFLAGS=-ansi -O2 -ip -xW -opt_report
SRCS = Lambert.cpp Vector.cpp Orbgnosis.cpp
OBJS:=$(SRCS:.cpp=.o)
MAIN=icc_orb
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o *.core $(MAIN)
