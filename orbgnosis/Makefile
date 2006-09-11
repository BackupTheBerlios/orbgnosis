#
# This makefile is free software; you can redistribute it and/or
# modify it under the terms of the enclosed JBQ code license.
# See LICENSE.Makefile
# 
# This makefile is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# License for more details.
#
# $Id: Makefile,v 1.9 2006/09/11 01:10:49 trs137 Exp $

# The file name, ...
BINARY      := orbgnosis
VERSION     := 0.2
BUILDDATE   := $(shell date +%m%d%Y)
BUILD_TYPE  := CVS

# output directories

OBJDIR  := obj
DEPDIR  := dep
SRCDIR  := src
DOCDIR  := doc

TEMPDEPFILE := /tmp/orbgnosis_temp_deps

# compiler, linker, ...

CC := g++ -c -pipe

# These can be set from the shell to change the build type
ifeq ($(BUILD_DEBUG),true)
    CFLAGS += -g -O0 
    BUILD_TYPE = Debug
else
# cpu-type "native" only supported in gcc 4.2.
    CFLAGS += -O2
	#CFLAGS += -march=native -mtune=native
    #CFLAGS += -mtune=G4
    CFLAGS += -Wdisabled-optimization
endif

CFLAGS += -DVERSION_STRING=\"$(BUILDDATE)\"
CFLAGS += -Wredundant-decls -Winline -Weffc++
CFLAGS += -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pedantic
CFLAGS += -Wno-multichar -Wpointer-arith -Wwrite-strings -Wconversion
CFLAGS += -Wcast-qual -Wshadow -Wabi -Wreorder
# The following are not supported in gcc < 4.0
#CFLAGS += -Wstrict-null-sentinel
# The following are not supported gcc < 4.2.
#CFLAGS += -Wunsafe-loop-optimizations

# be able to add flags from the env at will without chaging the makefile
CFLAGS += $(ORBGNOSIS_ADDL_BUILD_FLAGS)

LD := g++

DEP := g++ -MM

# create the dep dir.
MAKEDEP := $(shell mkdir -p $(DEPDIR))

# external commands
TAR := tar -cjvf
MAKE:= gmake

################################################

# the engine

MAKEFILE := Makefile

FULLNAME := $(subst \ ,_,$(BINARY))-$(BUILDDATE)

BASESOURCES := $(shell cd $(SRCDIR) && ls -1 *.cpp)
SOURCES := $(addprefix $(SRCDIR)/,$(BASESOURCES))
OBJECTS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(basename $(BASESOURCES))))
DEPENDS := $(addprefix $(DEPDIR)/,$(addsuffix .d,$(basename $(BASESOURCES))))

BASEHEADERS := $(shell cd $(SRCDIR) && ls -1 *.h)
HEADERS := $(addprefix $(SRCDIR)/,$(BASEHEADERS))

.PHONY : default release sourcearchive clean all

.DELETE_ON_ERROR : $(BINARY)

default : build

build : $(BINARY)

#   rule to create the object file directory if needed
$(OBJDIR)::
	@[ -d $(OBJDIR) ] || mkdir $(OBJDIR) > /dev/null 2>&1

clean :
	@echo cleaning
	@rm -rf $(BINARY) $(OBJDIR) $(DEPDIR) $(filter-out data.tbz,$(wildcard *.tbz)) *.tbz~

all : build sourcearchive doxygen pdfmanual pdfsource

sourcearchive: $(FULLNAME)-src.tbz

pdfsource: $(FULLNAME)-listing.pdf

cvsup:
	@echo CVS up -d
	@cvs up -d

$(BINARY) : $(OBJDIR) $(OBJECTS)
	@rm -f $(TEMPDEPFILE)
	@echo linking $@
	@$(LD) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@echo compiling [$(BUILD_TYPE) Build] : $@
	@$(CC) $< $(CFLAGS) -o $@

$(DEPDIR)/%.d : $(SRCDIR)/%.cpp
	@/bin/echo generating dependencies for $<
	@/bin/echo $@ $(OBJDIR)/$(shell $(DEP) $<) > $(TEMPDEPFILE)
	@/bin/echo $(OBJDIR)/$(basename $(@F))".o : $(MAKEFILE)" | /bin/cat - $(TEMPDEPFILE) > $@

$(FULLNAME)-src.tbz : $(SOURCES) $(HEADERS) $(MAKEFILE)
	@rm -rf $@~
	@mkdir -p $@~/$(FULLNAME)-src/
	@cp -r $(SRCDIR) $(MAKEFILE) $@~/$(FULLNAME)-src/
	@cp LICENSE.orbgnosis LICENSE.Makefile $@~/$(FULLNAME)-src/
	@rm -rf $@~/$(FULLNAME)-src/src/CVS
	@cd $@~ && $(TAR) $@ $(FULLNAME)-src
	@mv -f $@~/$@ .
	@rm -rf $@~

doxygen: Doxyfile $(SOURCES) $(HEADERS)
	@echo generating documentation with Doxygen
	@doxygen Doxyfile

pdfmanual: Doxyfile $(DOCDIR)/latex/Makefile
	@echo generating documentation with Doxygen
	@doxygen Doxyfile
	@echo running pdflatex...
	@cd $(DOCDIR)/latex && $(MAKE) > /dev/null 2>&1
	@mv -f $(DOCDIR)/latex/refman.pdf .
	@mv refman.pdf refman-$(BUILDDATE).pdf
	@cd $(DOCDIR)/latex && $(MAKE) clean > /dev/null 2>&1

$(FULLNAME)-listing.pdf: $(SOURCES) $(HEADERS) $(MAKEFILE)
	@echo generating formatted postscript from source
	@a2ps --header="$(FULLNAME)" -g -2 --toc $(SRCDIR)/* > $(FULLNAME)-listing.ps
	@ps2pdf $(FULLNAME)-listing.ps
	@rm -f $(FULLNAME)-listing.ps

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),sourcearchive)
ifneq ($(MAKECMDGOALS),doxygen)
ifneq ($(MAKECMDGOALS),pdfmanual)
ifneq ($(MAKECMDGOALS),pdfsource)
-include $(DEPENDS)
endif
endif
endif
endif
endif
endif
