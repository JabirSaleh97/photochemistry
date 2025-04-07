#############
### Targets #
#############

DEST		= /opt/ciopt/bin

LIBRARY		= 

PROGBASE	= CIOpt

PROGRAM		= $(PROGBASE).$(VERNO).e

TARGET		= $(PROGRAM)

# PRESET VARIABLES
#################################
### Compilers, Programs and Flags
#################################

### C Compiler
CC              = gcc
CFLAGS		= 

### Fortran Compiler
#FC		= pgf90
FC		= ifort -align
FFLAGS		= 

### Linker/Loader
#LD		= pgf90
LD		= ifort
LDFLAGS		= 

### Shell Command
SHELL		= /bin/sh

### Installation Program
INSTALL		= cp

### Print command
PRINT		= pr

.f.o:
	$(FC) -c $(FFLAGS) $*.f

#################################
### Sources Headers and Objects #
#################################

SRCS		= ciopt.f \
		  ReadOutput.f \
		  WriteInput.f \
		  powellb.f \
		  runmp.f \
		  func.f \
		  writexyz.f \
		  EvalCons.f \
		  svn_version.f \
		  sdrive.f \
		  lbfgs.f

SYSHDRS		=

OBJS		= ciopt.o \
		  ReadOutput.o \
		  WriteInput.o \
		  powellb.o \
		  runmp.o \
		  func.o \
		  writexyz.o \
		  EvalCons.o \
		  svn_version.o \
		  sdrive.o \
		  lbfgs.o

MAKEFILE	= Makefile

#######################
### Program Libraries #
#######################

LIBS		=
########################
### Additional Targets #
########################

all:		$(TARGET)


$(PROGRAM):	$(OBJS) $(MAKEFILE)
		@rm -f $(PROGBASE).e
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@ln -s $(PROGRAM) $(PROGBASE).e
		@echo "done"

docs:		
		robodoc --src . --doc doc --multidoc --html --index --nodesc 

svn_version.f:	
		@echo    '      subroutine svn_version(revision)'      > svn_version.f
		@echo    '      character*20 revision               ' >> svn_version.f
		@echo -n '      revision = "'                         >> svn_version.f
		@svnversion -n .                                      >> svn_version.f
		@echo    '"'                                          >> svn_version.f
		@echo    '      return  '                             >> svn_version.f
		@echo    '      end  '                                >> svn_version.f

install:	$(PROGRAM)
		$(INSTALL) $(PROGRAM) $(DEST)
		@rm -f $(DEST)/$(PROGBASE).e
		@ln -s $(DEST)/$(PROGRAM) $(DEST)/$(PROGBASE).e






