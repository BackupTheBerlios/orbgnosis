#
# $Id: Makefile,v 1.14 2006/05/05 21:50:40 trs137 Exp $
#
# *** use gmake! ***
#
CC=icpc
LD=icpc
RM=/bin/rm -f
#CFLAGS=-Wall -Werror -ansi -pedantic -g
#CFLAGS=-ansi -O2 -ip -static -xN -opt_report
CFLAGS=-O2 -ip -static -xN
REPORT=-opt_report_level max -opt_report_phase all -opt_report_file optreport
#REPORT=
LDFLAGS=
SRCS = Lambert.cpp Vector.cpp Orbgnosis.cpp
OBJS:=$(SRCS:.cpp=.o)
MAIN=icc_orb
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp Global.h
	$(CC) $(CFLAGS) $(REPORT) -c $<

clean:
	$(RM) *.o *.core $(MAIN)
