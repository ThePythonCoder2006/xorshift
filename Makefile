CC := gcc

N := 32# default
SUBN := 32#default

ifeq ($(N),64)
	SUBN := $(N)
endif

DEFINE_N := -DNUM="($(N)U)" -DSNUM="($(SUBN)U)"

DB_SUFFIX :=_db

SRCDIR := src
ODIR := $(SRCDIR)\obj
IDIR := include
BINDIR := bin

LUTDIR := $(SRCDIR)/lut
DDIR := data

OUT := $(BINDIR)/xorF_$(N).exe
LUT_OUT := $(BINDIR)/lut_$(N).exe
OUT_DB := $(BINDIR)/xorF$(DB_SUFFIX)_$(N).exe

WARN_FLAGS := -Wall -Wextra -pedantic
DBFLAGS := -ggdb3 -gdwarf-5
FASTFLAGS := -O2
LFLAGS := -mavx2 -I$(IDIR)

CFLAGS := $(WARN_FLAGS) $(LFLAGS)

XORF_SRC := $(wildcard $(SRCDIR)/*.c)
XORF_OBJ := $(XORF_SRC:$(SRCDIR)/%.c=$(ODIR)/%_$(N).o)
XORF_DB_OBJ := $(XORF_SRC:$(SRCDIR)/%.c=$(ODIR)/%$(DB_SUFFIX)_$(N).o)

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

$(OUT): $(XORF_OBJ) | $(BINDIR)
	$(CC) $(DEFINE_N) $^ -o $@ $(CFLAGS)
pre:
	$(CC) $(DEFINE_N) $(SRCDIR)/square_matrix.c -o mat_pre.c -E $(CFLAGS)

$(ODIR)/%_$(N).o:$(SRCDIR)/%.c Makefile $(IFILES)| $(ODIR)
	$(CC) $(DEFINE_N) $< -c -o $@ $(CFLAGS)

# ---------------- db --------------------

db: $(OUT_DB)
	gdb ./$<
comp_db: $(OUT_DB)

$(OUT_DB): $(XORF_DB_OBJ) | $(BINDIR)
	$(CC) $(DEFINE_N) $^ -o $@ $(CFLAGS) $(DBFLAGS)

$(ODIR)/%$(DB_SUFFIX)_$(N).o:$(SRCDIR)/%.c Makefile $(IFILES)| $(ODIR)
	$(CC) $(DEFINE_N) $< -c -o $@ $(CFLAGS) $(DBFLAGS)

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