#
# $Id: Makefile,v 1.27 2006/06/09 00:07:12 trs137 Exp $
#
# *** use gmake! ***
#
#CC=icpc
#LD=icpc
CC=g++42
LD=g++42
RM=/bin/rm -f
#CFLAGS=-Wall -O2 -ip -static -xN
CFLAGS=-g -Wall -Wextra -Wfloat-equal -Wshadow -Wunsafe-loop-optimizations -Wredundant-decls -Winline -Weffc++ -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pedantic -DNO_ICC_981
#REPORT=-opt_report_level min -opt_report_phase all -opt_report_file optreport
REPORT=
LDFLAGS=
SRCS = Orbgnosis.cpp Vector.cpp ULambert.cpp BLambert.cpp Traj.cpp Body.cpp Sat.cpp
OBJS:=$(SRCS:.cpp=.o)
MAIN=orb
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp Orbgnosis.h
	$(CC) $(CFLAGS) $(REPORT) -c $<

clean:
	$(RM) *.o *.core $(MAIN)
