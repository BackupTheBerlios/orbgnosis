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
# $Id: Makefile,v 1.1 2006/06/12 20:11:19 trs137 Exp $

# The file name, ...
BINARY      := orbgnosis
VERSION     := 0.2
BUILDDATE   := $(shell date +%b_%d_%Y)
BUILD_TYPE  := CVS

# output directories

OBJDIR  := obj
DEPDIR  := dep
SRCDIR  := src
DATADIR := data

TEMPDEPFILE := /tmp/orbgnosis_temp_deps

# compiler, linker, ...

CC := g++42 -c -pipe

# These can be set from the shell to change the build type
ifeq ($(BUILD_DEBUG),true)
    CFLAGS += -g -O0
    BUILD_TYPE = Debug
else
    CFLAGS += -O2
endif

CFLAGS := -DVERSION_STRING=\"$(VERSION)\"
CFLAGS += -march=pentium4 -mtune=pentium4
CFLAGS += -Wredundant-decls -Winline -Weffc++ -Wstrict-null-sentinel
CFLAGS += -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pedantic
CFLAGS += -Wno-multichar -Wpointer-arith -Wwrite-strings -Wconversion
CFLAGS += -Wcast-qual -Wshadow

# be able to add flags from the env at will without chaging the makefile
CFLAGS += $(ORBGNOSIS_ADDL_BUILD_FLAGS)

LD := g++42

DEP := g++42 -MM

# create the dep dir.
MAKEDEP := $(shell mkdir -p $(DEPDIR))

################################################

# the engine

MAKEFILE := Makefile

FULLNAME := $(subst \ ,_,$(BINARY))-$(VERSION)

BASESOURCES := $(shell cd $(SRCDIR) && ls -1 *.cpp)
SOURCES := $(addprefix $(SRCDIR)/,$(BASESOURCES))
OBJECTS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(basename $(BASESOURCES))))
DEPENDS := $(addprefix $(DEPDIR)/,$(addsuffix .d,$(basename $(BASESOURCES))))

BASEHEADERS := $(shell cd $(SRCDIR) && ls -1 *.h)
HEADERS := $(addprefix $(SRCDIR)/,$(BASEHEADERS))

.PHONY : default release clean binarchive sourcearchive all

.DELETE_ON_ERROR : $(BINARY)

default : build

build : $(BINARY)

#   rule to create the object file directory if needed
$(OBJDIR)::
	@[ -d $(OBJDIR) ] || mkdir $(OBJDIR) > /dev/null 2>&1

clean :
	@echo cleaning
	@rm -rf $(BINARY) $(OBJDIR) $(DEPDIR) $(filter-out data.zip,$(wildcard *.zip)) *.zip~

all : sourcearchive binarchive

sourcearchive : $(FULLNAME)-src.zip

binarchive : $(FULLNAME).zip

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

$(FULLNAME).zip : $(BINARY) $(DATA) $(EXTRADATA) $(MAKEFILE)
	@rm -rf $@~
	@mkdir -p $@~/$(FULLNAME)/
ifneq "$(EXTRADATA)" ""
	@unzip $(EXTRADATA) -d $@~/$(FULLNAME)/
endif
	@rm -rf $@~/$(FULLNAME)/CVS $@~/$(FULLNAME)/protocol/CVS $@~/$(FULLNAME)/help/CVS $@~/$(FULLNAME)/scripts/CVS
	@rm -rf $@~/$(FULLNAME)/CVS $@~/$(FULLNAME)/help/content/CVS $@~/$(FULLNAME)/help/images/CVS
	@cd $@~ && $(ZIP) $@ $(FULLNAME)
	@mv -f $@~/$@ .
	@rm -rf $@~

$(FULLNAME)-src.zip : $(SOURCES) $(HEADERS) $(DATA) $(EXTRADATA) $(MAKEFILE)
	@rm -rf $@~
	@mkdir -p $@~/$(FULLNAME)-src/
	@rm -rf $@~/$(FULLNAME)-src/src/CVS $@~/$(FULLNAME)-src/data/CVS
	@rm -rf $@~/$(FULLNAME)-src/CVS $@~/$(FULLNAME)-src/data/protocol/CVS $@~/$(FULLNAME)-src/data/help/CVS $@~/$(FULLNAME)-src/data/scripts/CVS
	@rm -rf $@~/$(FULLNAME)-src/CVS $@~/$(FULLNAME)-src/data/help/content/CVS $@~/$(FULLNAME)-src/data/help/images/CVS
	@cd $@~ && $(ZIP) $@ $(FULLNAME)-src
	@mv -f $@~/$@ .
	@rm -rf $@~

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),sourcearchive)
-include $(DEPENDS)
endif
endif
endif
