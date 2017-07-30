
CFLAGS = $(shell gsl-config --cflags)
CFLAGS += -O3 -I include
LIBS   = $(shell gsl-config --libs)
CC     = g++
BDIR   = bin
ODIR   = obj
SDIR   = subroutines
EDIR   = EXEs
OS    := $(shell uname)

ifeq ($(OS),Linux)
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
else
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
endif

_EXECUTABLES = rd_gen_reweight \
	       rd_gen_reweightPBC \
	       rd_AB_reweight_cellPBC_3Drep \
	       rd_self_reweight_cellPBC \
	       rd_trap_rewgt_3D \
	       rd_AB_reweight_cellPBC_3Dskip \
	       rd_self_reweight_cellPBC_skip \
	       rd_AB_reweight_cellPBC_2Dcomsig \
	       rd_gen_reweight_FULLCELL \
	       rd_gen_reweight_FULLCELL_AVOID \
	       rd_gen_reweight_FULLCELL_AVOID_subs \
	       rd_gen_reweight_FULLCELL_complex \
	       rd_gen_reweight_FULLCELL_AVOID_MF \
	       rd_gen_noreweight_FULLCELL_AVOID \
               coordsgen \
	       rd_Vr_irrev_reweightPBC_Force \
	       rd_self_reweight_cellPBC_3Dirr \
	       rd_self_reweight_cellPBC_3Dirr_abs \
	       rd_trap_trueGF_3D \
	       rd_trap_rewgt_2D \
	       rd_trap_trueGF_2D \
	       rd_AB_reweight_cellPBC_2Dirr \
	       rd_self_reweight_cellPBC_2Dirr \
	       rd_Vr_irrev_reweightPBC_2fit_dens \
	       rd_AB_reweight_cellPBC_2Dskip \
	       rd_self_reweight_cellPBC_2Dskip \
	       rd_self_reweight_cellPBC_3Dirr_absNEWINPUTS \
	       rd_self_reweight_cellPBC_3DirrNEWINPUTS \


OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

EXECUTABLES = $(patsubst %,$(BDIR)/%,$(_EXECUTABLES))

_SOURCES = ${_EXECUTABLES:=.cpp}
SOURCES = $(patsubst %,$(EDIR)/%,$(_SOURCES))

all: dirs $(EXECUTABLES)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) $(CFLAGS2) -c $< -o $@  
	@echo "------------"

$(EXECUTABLES): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS)
	@echo "------------"

dirs: 
	mkdir -p bin
	mkdir -p obj

clean:
	rm -f $(ODIR)/*.o $(EXECUTABLES)


