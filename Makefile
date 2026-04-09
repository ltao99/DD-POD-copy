# ======================= Paths ============================
PATHLIB       = /clonetroop/dasilva/libs
PATHCOMPILER  = /clonetroop/dasilva/libs
PATHMUMPS     = $(PATHLIB)
PATHMETIS     = $(PATHLIB)
PATHMPI       = $(PATHCOMPILER)/bin

MPI_BIN       = $(PATHMPI)
MPI_LIB       = $(PATHLIB)/lib
MPI_INCLUDE   = $(PATHLIB)/include

LPORDDIR      = $(PATHMUMPS)/lib
IPORD         = -I$(PATHMUMPS)/include
LPORD         = -L$(LPORDDIR) -lpord

LMETISDIR     = $(PATHMETIS)/lib
IMETIS        = -I$(PATHMETIS)/include
LMETIS        = -L$(LMETISDIR) -lparmetis -lmetis

# If using SCOTCH, define here (comment out if not used)
LSCOTCH       = -lscotch
LORDERINGS    = $(LMETIS) $(LPORD) $(LSCOTCH)

PATHLASC      = $(PATHLIB)
LIB_SPBLAS    = $(PATHLIB)

# Intel MKL path
MKLROOT ?= /clonetroop/software/intel/oneapi/mkl/latest

# ======================= Compiler and Flags ==============
FC       = $(MPI_BIN)/mpif90
MPIRUN   = $(MPI_BIN)/mpirun
FFLAGS   = -O3 -march=native -funroll-loops -I$(MKLROOT)/include
INCLUDES = -I$(MPI_INCLUDE) -I$(LIB_SPBLAS)/include $(IPORD) $(IMETIS)
LIBDIRS  = -L$(MPI_LIB) -L$(PATHLASC)/lib -L$(PATHCOMPILER)/lib -Wl,-rpath,$(MPI_LIB)

LIBS = -Wl,--start-group \
  -lzmumps -lmumps_common -lpord \
  -L$(MKLROOT)/lib/intel64 \
  $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
  $(MKLROOT)/lib/intel64/libmkl_sequential.a \
  $(MKLROOT)/lib/intel64/libmkl_core.a \
  $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
  -lpthread -lm \
  -lparmetis -lmetis \
  -lmpi \
  -Wl,--end-group


# ======================= Source Files =====================
VPATH = src
SRC = \
  sizes.f90 \
  MOD_Constants.f90 \
  commons_3d.f90 \
  set_sub_coordinates.f90 \
  read_basis.f90 \
  string_utils.f90 \
  model.f90 \
  MOD_Quicksort.f90 \
  readata.f90 \
  set_nodes.f90 \
  msource.f90 \
  set_bet_del_alp.f90 \
  inguess.f90 \
  matrix_handling.f90 \
  mat_coef_global.f90 \
  rhs.f90 \
  solver_MUMPS.f90 \
  transm.f90 \
  pri_fields_newz0.f90 \
  r-b-proc.f90 \
  calc_error.f90 \
  l2_utils.f90 \
  output.f90 \
  max3d-par.f90

OBJ = $(SRC:.f90=.o)
EXE = max3d-par

# ======================= Build Rules ======================
all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBDIRS) -o $@ $^ $(LIBS)

%.o: src/%.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

# ======================= Run Rule =========================
NP ?= 18

run: $(EXE)
	OMP_NUM_THREADS=1 LD_LIBRARY_PATH=$(MPI_LIB):$$LD_LIBRARY_PATH $(MPIRUN) -np $(NP) ./$(EXE)

# ======================= Clean Rule =======================
clean:
	rm -f *.o *.mod *.i *.xyz fort.* $(EXE)






