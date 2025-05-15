# all compiled code will be stored in gen directory
BUILD_DIR = gen
# where all source code is located
SOURCE_DIR = src

# on my machine it is gcc-14 (regular gcc can be used instead)
OMPCC = gcc -fopenmp
CC_FLAGS = -O3
CC_FLAGS += -Wall -Wextra
LIBS = -lm

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif


# $(^) stands for everything written after: for ex. $(BUILD_DIR)/prime: -> all dependencies
# -o defines name of output file
# $(@) stands for target -> written before : -> $(BUILD_DIR)/prime

# all is defined as main target when running make
all: $(BUILD_DIR)/feynman_omp_1d $(BUILD_DIR)/feynman_omp_2d $(BUILD_DIR)/feynman_omp_3d \
	 $(BUILD_DIR)/feynman_pthreads_1d $(BUILD_DIR)/feynman_pthreads_3d

# OpenMP
$(BUILD_DIR)/feynman_omp_1d: $(SOURCE_DIR)/feynman_omp_1d.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/feynman_omp_2d: $(SOURCE_DIR)/feynman_omp_2d.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/feynman_omp_3d: $(SOURCE_DIR)/feynman_omp_3d.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)


#Pthreads
$(BUILD_DIR)/feynman_pthreads_1d: $(SOURCE_DIR)/feynman_pthreads_1d.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS) -lpthread

$(BUILD_DIR)/feynman_pthreads_3d: $(SOURCE_DIR)/feynman_pthreads_3d.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS) -lpthread

$(BUILD_DIR):		# when running this, execute this command: this command will run if $(BUILD_DIR) does not exist, -p adds parent directories in the path of the new one
	mkdir -p $(BUILD_DIR)

# removing gen directory
clean:
	rm -rf $(BUILD_DIR)
