CC := gcc

N := 32# default
SUBN := 32#default

ifeq ($(N),64)
	SUBN := $(N)
endif

DEFINE_N = -DNUM="($(N)U)" -DSNUM="($(SUBN)U)"

ifneq ($(N), 32)
	ifeq ($(SUB_MULT), 1)
		DEFINE_N := $(DEFINE_N) -DUSE_SUBMAT_MULT
	endif
endif

DB_SUFFIX :=_db

SRCDIR := src
ODIR := $(SRCDIR)\obj
IDIR := include
BINDIR := bin

LUTDIR := $(SRCDIR)/lut
DDIR := data

ifeq ($(N), $(SUBN))
	FNAME_SUFFIX := _$(N)
else
	FNAME_SUFFIX := _$(N)_$(SUBN)
	ifeq ($(SUB_MULT), 1)
		SUBMULT_SUFFIX := _submult
	endif
endif


LUT_OUT := $(BINDIR)/lut_$(SUBN).exe
BASE_OUT_FNAME := $(BINDIR)/xorF$(SUBMULT_SUFFIX)
OUT := $(BASE_OUT_FNAME)$(FNAME_SUFFIX).exe
OUT_DB := $(BASE_OUT_FNAME)$(DB_SUFFIX)$(FNAME_SUFFIX).exe

WARN_FLAGS := -Wall -Wextra -pedantic
DBFLAGS := -ggdb3 -gdwarf-5
FASTFLAGS := -O2
LFLAGS := -mavx2 -I$(IDIR)

CFLAGS := $(WARN_FLAGS) $(LFLAGS)

BASE_OBJ_FNAME := $(ODIR)/%$(SUBMULT_SUFFIX)
XORF_SRC := $(wildcard $(SRCDIR)/*.c)
XORF_OBJ := $(XORF_SRC:$(SRCDIR)/%.c=$(BASE_OBJ_FNAME)$(FNAME_SUFFIX).o)
XORF_DB_OBJ := $(XORF_SRC:$(SRCDIR)/%.c=$(BASE_OBJ_FNAME)$(DB_SUFFIX)$(FNAME_SUFFIX).o)

IFILES := $(wildcard $(IDIR)/*.h)

DFILES_ := lut_32.data lut_64.data
DFILES := $(DFILES_:%.data=$(DDIR)/%.data)
ifeq ($(SUBN),32)
	DFILE = $(word 1,$(DFILES))
else
	DFILE = $(word 2,$(DFILES))
endif

ifeq ($(OS),Windows_NT)
	RM = cmd /C del /Q /F
	RRM = - cmd /C rmdir /Q /S
else
	RM = rm -f
	RRM = rm -f -r
endif

default: fast

run db: | $(DFILE)

run: $(OUT)
	./$^

plain: $(OUT) run

fast: comp_fast run
comp_fast: set_fast $(OUT)
set_fast:
	$(eval CFLAGS += $(FASTFLAGS))

run_db: $(OUT_DB)
	gdb ./$<
db: comp_db run_db
comp_db: set_db $(OUT_DB)
set_db:
	$(eval CFLAGS += $(DBFLAGS))

$(OUT): $(XORF_OBJ)
$(OUT_DB): $(XORF_DB_OBJ)

$(OUT) $(OUT_DB): | $(BINDIR)
	$(CC) $(DEFINE_N) $^ -o $@ $(CFLAGS)
pre:
	$(CC) $(DEFINE_N) $(SRCDIR)/square_matrix.c -o mat_pre.c -E $(CFLAGS)

$(BASE_OBJ_FNAME)$(FNAME_SUFFIX).o $(BASE_OBJ_FNAME)$(DB_SUFFIX)$(FNAME_SUFFIX).o:$(SRCDIR)/%.c Makefile $(IFILES)| $(ODIR)
	$(CC) $(DEFINE_N) $< -c -o $@ $(CFLAGS)


# ---------------- lut -------------------

lut: $(DFILE)

$(DFILE): $(LUT_OUT) Makefile include/constants.h
	./$<

$(LUT_OUT): $(LUTDIR)/generate_lut.c | data
	$(CC) $(DEFINE_N) $^ -o $@ $(CFLAGS)


$(DDIR) $(ODIR) $(BINDIR):
	- mkdir "$@"

clean: 
	$(RRM) $(DDIR)
	$(RRM) $(ODIR)
	$(RRM) $(BINDIR)