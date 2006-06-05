#
# $Id: Makefile,v 1.21 2006/06/05 16:45:40 trs137 Exp $
#
# *** use gmake! ***
#
CC=icpc
LD=icpc
#CC=/usr/local/bin/g++42
#LD=/usr/local/bin/g++42
RM=/bin/rm -f
CFLAGS=-Wall -pedantic -g -DNO_ICC_981
#CFLAGS=-Wall -O2 -ip -static -xN
#REPORT=-opt_report_level min -opt_report_phase all -opt_report_file optreport
REPORT=
LDFLAGS=
SRCS = ULambert.cpp BLambert.cpp Vector.cpp Orbgnosis.cpp
OBJS:=$(SRCS:.cpp=.o)
MAIN=orb
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp Orbgnosis.h
	$(CC) $(CFLAGS) $(REPORT) -c $<

clean:
	$(RM) *.o *.core $(MAIN)
