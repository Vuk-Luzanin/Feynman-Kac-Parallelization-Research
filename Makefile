
# environment variables set:
# export OMP_WAIT_POLICY=ACTIVE

# TODO(vuk): add profiler 

# all compiled code will be stored in gen directory
BUILD_DIR = result

# per-strategy build directories
OMP_BUILD_DIR = $(BUILD_DIR)/feynman_omp
PTHREADS_BUILD_DIR = $(BUILD_DIR)/feynman_pthreads


# where all source code is located
SOURCE_DIR = src
OMP_DIR = $(SOURCE_DIR)/OpenMP
PTHREADS_DIR = $(SOURCE_DIR)/Pthreads

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

# all is defined as main target when running make
all: $(OMP_BUILD_DIR)/feynman_omp_1d $(OMP_BUILD_DIR)/feynman_omp_2d $(OMP_BUILD_DIR)/feynman_omp_3d \
	 $(PTHREADS_BUILD_DIR)/feynman_pthreads_1d $(PTHREADS_BUILD_DIR)/feynman_pthreads_2d $(PTHREADS_BUILD_DIR)/feynman_pthreads_3d

# OpenMP
$(OMP_BUILD_DIR)/feynman_omp_1d: $(OMP_DIR)/feynman_omp_1d.c $(SOURCE_DIR)/util.c | $(OMP_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(OMP_BUILD_DIR)/feynman_omp_2d: $(OMP_DIR)/feynman_omp_2d.c $(SOURCE_DIR)/util.c | $(OMP_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(OMP_BUILD_DIR)/feynman_omp_3d: $(OMP_DIR)/feynman_omp_3d.c $(SOURCE_DIR)/util.c | $(OMP_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)


#Pthreads
$(PTHREADS_BUILD_DIR)/feynman_pthreads_1d: $(PTHREADS_DIR)/feynman_pthreads_1d.c $(SOURCE_DIR)/util.c | $(PTHREADS_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS) -lpthread

$(PTHREADS_BUILD_DIR)/feynman_pthreads_2d: $(PTHREADS_DIR)/feynman_pthreads_2d.c $(SOURCE_DIR)/util.c | $(PTHREADS_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS) -lpthread

$(PTHREADS_BUILD_DIR)/feynman_pthreads_3d: $(PTHREADS_DIR)/feynman_pthreads_3d.c $(SOURCE_DIR)/util.c | $(PTHREADS_BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS) -lpthread

# create needed directories
$(OMP_BUILD_DIR):
	mkdir -p $@

$(PTHREADS_BUILD_DIR):
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

