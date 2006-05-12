#
# $Id: Makefile,v 1.16 2006/05/12 23:51:17 trs137 Exp $
#
# *** use gmake! ***
#
#CC=icpc
#LD=icpc
CC=/usr/local/bin/g++42
LD=/usr/local/bin/g++42
RM=/bin/rm -f
#CFLAGS=-Wall -pedantic -g
#CFLAGS=-Wall -O2 -ip -static -xN
CFLAGS=-Wall -O2 -static
#REPORT=-opt_report_level min -opt_report_phase all -opt_report_file optreport
REPORT=
LDFLAGS=
SRCS = Lambert.cpp Vector.cpp Orbgnosis.cpp
OBJS:=$(SRCS:.cpp=.o)
MAIN=orb
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) $(REPORT) -c $<

clean:
	$(RM) *.o *.core $(MAIN)
