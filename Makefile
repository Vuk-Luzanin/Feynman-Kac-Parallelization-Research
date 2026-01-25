# environment variables set:
# export OMP_WAIT_POLICY=ACTIVE

# TODO(vuk): add profiler 

# all compiled code will be stored in gen directory
BUILD_DIR = result

# per-strategy build directories
SEQUENTIAL_BUILD_DIR = $(BUILD_DIR)/feynman_sequential
OMP_BUILD_DIR = $(BUILD_DIR)/feynman_omp
PTHREADS_BUILD_DIR = $(BUILD_DIR)/feynman_pthreads
# applications build directory
APPLICATIONS_BUILD_DIR = $(BUILD_DIR)/feynman_applications
APPLICATION_HEAT_EQUATION_BUILD_DIR = $(APPLICATIONS_BUILD_DIR)/heat_equation
APPLICATION_GIRSANOV_BUILD_DIR = $(APPLICATIONS_BUILD_DIR)/girsanov_importance_sampling


# where all source code is located
SOURCE_DIR = src
SEQUENTIAL_DIR = $(SOURCE_DIR)/base
OMP_DIR_SRC = $(SOURCE_DIR)/OpenMP
PTHREADS_DIR_SRC = $(SOURCE_DIR)/Pthreads

APPLICATIONS_DIR_SRC = $(SOURCE_DIR)/applications
APPLICATIONS_HEAT_EQUATION_DIR_SRC = $(APPLICATIONS_DIR_SRC)/heat_equation
APPLICATIONS_GIRSANOV_DIR_SRC = $(APPLICATIONS_DIR_SRC)/girsanov_importance_sampling

# on my machine it is gcc-14 (regular gcc can be used instead)
OMPCC = gcc -fopenmp
# Ofast -> O3 + -ffast-math	
CC_FLAGS = -Ofast
# optimization while linking
CC_FLAGS += -flto
# -march=native -> march (machine architecture - to be native) - finds characteristics of my cpu and uses all its instruction
# -ftree-vectorize -> enables vectorization (use of SIMD instructions)
CC_FLAGS += -march=native
# -funroll-loops -> unroll the loops
CC_FLAGS += -funroll-loops
# ignore uninitialized variable warnings
CC_FLAGS += -Wno-maybe-uninitialized

# preuredjuej redosled i organizaciju ugnezdenih petlji
CC_FLAGS += -floop-interchange -floop-block -floop-strip-mine
CC_FLAGS += -fprefetch-loop-arrays
CC_FLAGS += -Wall -Wextra 
# include path so util.h can be found -> to know where to find .h files
CC_FLAGS += -I$(SOURCE_DIR)  
LIBS = -lm

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif


# $(^) stands for everything written after: for ex. $(BUILD_DIR)/prime: -> all dependencies
# -o defines name of output file
# $(@) stands for target -> written before : -> $(BUILD_DIR)/prime

DIM = 1d 2d 3d

# all is defined as main target when running make
all: simulations apps

# build all feynman-kac implementations
simulations: \
  $(foreach d,$(DIM),$(SEQUENTIAL_BUILD_DIR)/feynman_sequential_$(d)) \
  $(foreach d,$(DIM),$(OMP_BUILD_DIR)/feynman_omp_$(d)) \
  $(foreach d,$(DIM),$(PTHREADS_BUILD_DIR)/feynman_pthreads_$(d)) 

# sequential
$(SEQUENTIAL_BUILD_DIR)/feynman_sequential_%: $(SEQUENTIAL_DIR)/feynman_sequential_%.c | $(SEQUENTIAL_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS)

# OpenMP
$(OMP_BUILD_DIR)/feynman_omp_%: $(OMP_DIR_SRC)/feynman_omp_%.c $(SOURCE_DIR)/util.c | $(OMP_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS)

# Pthreads
$(PTHREADS_BUILD_DIR)/feynman_pthreads_%: $(PTHREADS_DIR_SRC)/feynman_pthreads_%.c $(SOURCE_DIR)/util.c | $(PTHREADS_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS) -lpthread

# buiild applications
apps: \
  $(APPLICATION_HEAT_EQUATION_BUILD_DIR)/heat_equation_sequential \
  $(APPLICATION_HEAT_EQUATION_BUILD_DIR)/heat_equation_omp \
  $(APPLICATION_GIRSANOV_BUILD_DIR)/girsanov_importance_sampling_sequential \
  $(APPLICATION_GIRSANOV_BUILD_DIR)/girsanov_importance_sampling_omp

# Applications: Heat Equation - sequential
$(APPLICATION_HEAT_EQUATION_BUILD_DIR)/heat_equation_sequential: $(APPLICATIONS_HEAT_EQUATION_DIR_SRC)/heat_equation_sequential.c $(SOURCE_DIR)/util.c | $(APPLICATION_HEAT_EQUATION_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS) 

# Applications: Heat Equation - omp
$(APPLICATION_HEAT_EQUATION_BUILD_DIR)/heat_equation_omp: $(APPLICATIONS_HEAT_EQUATION_DIR_SRC)/heat_equation_omp.c $(SOURCE_DIR)/util.c | $(APPLICATION_HEAT_EQUATION_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS) 

# Applications: Girsanov - sequential
$(APPLICATION_GIRSANOV_BUILD_DIR)/girsanov_importance_sampling_sequential: $(APPLICATIONS_GIRSANOV_DIR_SRC)/girsanov_importance_sampling_sequential.c $(SOURCE_DIR)/util.c | $(APPLICATION_GIRSANOV_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS) 

# Applications: Girsanov - omp
$(APPLICATION_GIRSANOV_BUILD_DIR)/girsanov_importance_sampling_omp: $(APPLICATIONS_GIRSANOV_DIR_SRC)/girsanov_importance_sampling_omp.c $(SOURCE_DIR)/util.c | $(APPLICATION_GIRSANOV_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $^ -o $@ $(LIBS) 


# create needed directories
$(SEQUENTIAL_BUILD_DIR) $(OMP_BUILD_DIR) $(PTHREADS_BUILD_DIR) $(APPLICATIONS_BUILD_DIR) $(APPLICATION_HEAT_EQUATION_BUILD_DIR) $(APPLICATION_GIRSANOV_BUILD_DIR):
	mkdir -p $@

# removing gen directory
clean:
	rm -rf $(BUILD_DIR)

# work with profiler:
#  add  -fprofile-generate to CC_FLAGS
# compile
# this will generate gmon.out file
# to see that file run:
# - gprof ./result/TEST/feynman_pthreads_1d  gmon.out > analysis.txt
# then, delete -fprofile-generate flag, and add -fprofile-use instead
