CC = g++41
LD = g++41
#CFLAGS =  -O0 -Wredundant-decls -Winline -Weffc++
#CFLAGS += -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pedantic
#CFLAGS += -Wno-multichar -Wpointer-arith -Wwrite-strings -Wconversion
#CFLAGS += -Wcast-qual -Wshadow -Wabi -Wreorder

SRCS = src/Orbgnosis.cpp src/Tour.cpp src/Graph.cpp src/Vec3.cpp
SRCS += src/nsga2/allocate.c src/nsga2/decode.c src/nsga2/fillnds.c src/nsga2/mutation.c src/nsga2/rank.c
SRCS += src/nsga2/auxiliary.c src/nsga2/display.c src/nsga2/initialize.c src/nsga2/report.c
SRCS += src/nsga2/crossover.c src/nsga2/dominance.c src/nsga2/list.c src/nsga2/problemdef.c src/nsga2/sort.c
SRCS += src/nsga2/crowddist.c src/nsga2/eval.c src/nsga2/merge.c src/nsga2/rand.c src/nsga2/tourselect.c
#SRCS += src/nsga2/src/nsga2r.c

OBJS = $(SRCS:.cpp=.o)
PROG = test
RM = /bin/rm -f

$(PROG): $(SRCS) Makefile
	$(CC) $(CFLAGS) $(SRCS) -o $(PROG)

clean:
	$(RM) *.o *.core $(PROG)
