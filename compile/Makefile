# -- Makefile for sugar compiler (standalone and MEX file)
#
# To make the standalone file, type "make" or "make standalone"
# To make the mex file, type "make mex"
#
# Compile flags are:
#  USE_XT -- use the Xt toolkit for display
#  MEMPOOL_CLEAR_DEBUG -- clear memory after freeing, for debug purposes
#  HAS_LAPACK -- if BLAS and LAPACK are available
#  __SUGAR_MEX -- if a mex file is being compiled

MAKE=make

CC=gcc
CFLAGS=-g -Wall -DUSE_XT -DMEMPOOL_CLEAR_DEBUG -DHAS_LAPACK -Wall -pedantic
DRAWLIBS=  -L/usr/X11R6/lib/ -lXaw -lXt -lX11
LALIBS=    -llapack -lblas -lg2c
CLIBS= $(DRAWLIBS) $(LALIBS) -lnsl -lm

CMEX=/usr/sww/pkg/matlab-6.0/bin/mex
#CMEX=/usr/sww/pkg/matlab-5.3/bin/cmex
CMEX=mex
MEXFLAGS=-g -DMEMPOOL_CLEAR_DEBUG -D__SUGAR_MEX 
MEXLIBS=


# -----------------------------------------------------------------

OBJS = sugar_c.o \
       mempool.o shash.o symbol.o \
       sugar.tab.o lex.yy.o \
       parse.o codegen.o codestack.o evalexpr.o intrinsic.o \
       codeprint.o position.o varassign.o modelmgr.o \
       paramprint.o writegeom.o assembly.o dcsolve.o \
       xtnetdraw.o models.o

MXOBJS = sugar_c.o \
       mempool.o shash.o symbol.o \
       sugar.tab.o lex.yy.o \
       parse.o codegen.o codestack.o evalexpr.o intrinsic.o \
       codeprint.o codemx.o position.o varassign.o modelmgr.o \
       paramprint.o writegeom.o assembly.o

# -----------------------------------------------------------------

standalone:
	$(MAKE) sugar-c CMPL="$(CC)" CMPLFLAGS="$(CFLAGS)" LIBS="$(CLIBS)"

mex:
	$(MAKE) sugarmex CMPL="$(CMEX)" CMPLFLAGS="$(MEXFLAGS)" LIBS="$(MEXLIBS)"

# -----------------------------------------------------------------

.SUFFIXES : .c

.c.o:
	$(CMPL) -c $(CMPLFLAGS) $<

sugar-c: $(OBJS)
	$(CMPL) -o sugar-c $(OBJS) $(LIBS)

sugarmex: $(MXOBJS)
	$(CMPL) $(MXOBJS) 

sugar.tab.c: sugar.y
	bison -v -d sugar.y

lex.yy.c: sugar.lex
	flex sugar.lex
	cat lex.yy.c | awk '!/unistd.h/' > lex.yy-win.c

clean:
	rm -f *.o

clean-all:
	rm -f *.o sugar.tab.c sugar.tab.h sugar.output lex.yy.c \
              lex.yy-win.c sugar-c

count:
	wc -l codegen.[ch] codemx.[ch] codeprint.[ch] codestack.[ch] \
	evalexpr.[ch] intrinsic.[ch] \
        mempool.[ch] parse.[ch] position.[ch] shash.[ch] \
	symbol.[ch] sugar.lex sugar.y \
	xtnetdraw.[ch] modelmgr.[ch] models.[ch] varassign.[ch] \
        paramprint.[ch] \
	writegeom.[ch] dcsolve.[ch] assembly.[ch]

mxcopy:
	cp *mex* /home/dbindel/work/sugar/src-2.0/analysis

sugar_c.o: sugar_c.c
mempool.o: mempool.c
shash.o: shash.c
parse.o: parse.c parse.h symbol.h
symbol.o: symbol.c symbol.h
codegen.o: codegen.c codegen.h codestack.h evalexpr.h
position.o: position.c 
codeprint.o: codeprint.c codeprint.h codestack.h evalexpr.h
codestack.o: codestack.c codestack.h evalexpr.h
evalexpr.o: evalexpr.c evalexpr.h
sugar.tab.o: sugar.tab.c parse.h
lex.yy.o: lex.yy.c sugar.tab.c
lex.yy-win.o: lex.yy-win.c sugar.tab.c


