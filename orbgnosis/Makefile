# This makefile is free software; you can redistribute it and/or
# modify it under the terms of the enclosed JBQ code license.
# See LICENSE.Makefile
# 
# This makefile is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# License for more details.
#
# $Id: Makefile,v 1.19 2006/10/25 11:55:54 trs137 Exp $

# This makefile depends on GNU make and gcc.  Optionally, you may
# use a high performance commercial compiler by editing TARGET, but
# GNU make and gcc are always mandatory.
# Optionally, you will need doxygen, pdflatex, a2ps and ps2pdf if you
# wish to generate documentation or PDF program listings.
#
# Choose the platform and compiler. Valid choices are:
#
# TARGET := "default"     # GNU compiler, any architecture
# TARGET := "icc_p4"     # Intel compiler, Pentium 4 (arl1426.arl.psu.edu)
# TARGET := "hammer"      # Intel compiler, Opteron (hammer.aset.psu.edu)
# TARGET := "power"       # IBM VisualAge, Power5 (power.ics.psu.edu)
#
TARGET := "icc_p4"

# The file name, ...
BINARY      := orbgnosis
VERSION     := 0.9
BUILDDATE   := $(shell date +%m%d%Y)
BUILD_TYPE  := CVS
#BUILD_DEBUG := true

# output directories
OBJDIR := obj
DEPDIR := dep
SRCDIR := src
DOCDIR := doc
TEMPDEPFILE := /tmp/orbgnosis_temp_deps

# Set the compiler ...
CC := g++41 -c -pipe
CXX := g++41 -c -pipe
LD := g++41
DEP := g++41 -MM  # always use gcc to generate deps
ifeq ($(TARGET), "icc_p4")
    CC := icpc -c -pipe
    LD := icpc
endif
ifeq ($(TARGET), "hammer")
    CC := icpc -c -pipe
    LD := icpc
endif
ifeq ($(TARGET), "power")
    CC := xlC -c
    LD := xlC 
endif

# This can be set from the shell to change the build type
ifeq ($(BUILD_DEBUG),true)
    CFLAGS += -g -O0 
    BUILD_TYPE = Debug
else
    CFLAGS += -O3
    ifeq ($(TARGET), "icc_p4")
        CFLAGS += -Wall -Wcheck -w1 -Wp64
        CFLAGS += -align -ip -march=pentium4
        CFLAGS += -opt_report -opt_report_level=max
        CFLAGS += -opt_report_phase all
        CFLAGS += -opt_report_file=icc_opt_report
    endif
    ifeq ($(TARGET), "hammer")
        CFLAGS += -Wall -Wcheck -w1 -Wp64
        CFLAGS += -align -ip
        CFLAGS += -opt_report -opt_report_level=max
        CFLAGS += -opt_report_phase all
        CFLAGS += -opt_report_file=icc_opt_report
    endif
    ifeq ($(TARGET), "power")
        CFLAGS += -qwarn64
        #CFLAGS += -qformat
        #CFLAGS += -qinfo=all
        CFLAGS += -qarch=pwr5 -qtune=pwr5
        CFLAGS += -qhot -qreport
    endif
endif

ifeq ($(CC),"g++ -c -pipe")  
    CFLAGS += -Wredundant-decls -Winline -Weffc++
    CFLAGS += -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pedantic
    CFLAGS += -Wno-multichar -Wpointer-arith -Wwrite-strings -Wconversion
    CFLAGS += -Wcast-qual -Wshadow -Wabi -Wreorder
endif


CFLAGS += -DVERSION_STRING=\"$(BUILDDATE)\"

# You can add flags from the environment at will without chaging the makefile
CFLAGS += $(ORBGNOSIS_ADDL_BUILD_FLAGS)

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
BASESOURCES += $(shell cd $(SRCDIR) && ls -1 *.c)
SOURCES := $(addprefix $(SRCDIR)/,$(BASESOURCES))
OBJECTS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(basename $(BASESOURCES))))
DEPENDS := $(addprefix $(DEPDIR)/,$(addsuffix .d,$(basename $(BASESOURCES))))

BASEHEADERS := $(shell cd $(SRCDIR) && ls -1 *.h)
HEADERS := $(addprefix $(SRCDIR)/,$(BASEHEADERS))

.PHONY : default release sourcearchive clean all

.DELETE_ON_ERROR : $(BINARY)

default : build

build : $(BINARY)
	@echo Completed build for $(FULLNAME).  Please read through the output above
	@echo for any errors or warnings.
	@echo CC = $(CC)
	@echo CXX = $(CXX)
	@echo CFLAGS = $(CFLAGS)

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

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo compiling C file [$(BUILD_TYPE) Build] : $@
	@$(CC) $< $(CFLAGS) -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@echo compiling C++ file [$(BUILD_TYPE) Build] : $@
	@$(CC) $< $(CFLAGS) -o $@

$(DEPDIR)/%.d : $(SRCDIR)/%.c
	@echo generating C dependencies for $<
	@echo $@ $(OBJDIR)/$(shell $(DEP) $<) > $(TEMPDEPFILE)
	@echo $(OBJDIR)/$(basename $(@F))".o : $(MAKEFILE)" | /bin/cat - $(TEMPDEPFILE) > $@

$(DEPDIR)/%.d : $(SRCDIR)/%.cpp
	@echo generating C++ dependencies for $<
	@echo $@ $(OBJDIR)/$(shell $(DEP) $<) > $(TEMPDEPFILE)
	@echo $(OBJDIR)/$(basename $(@F))".o : $(MAKEFILE)" | /bin/cat - $(TEMPDEPFILE) > $@

$(FULLNAME)-src.tbz : $(SOURCES) $(HEADERS) $(MAKEFILE)
	@rm -rf $@~
	@mkdir -p $@~/$(FULLNAME)-src/
	@cp -r $(SRCDIR) data $(MAKEFILE) $@~/$(FULLNAME)-src/
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
