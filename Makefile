# Compiler and Linker settings
CXX = g++
# Add include path for the sdsl library to compiler flags
CXXFLAGS = -std=c++17 -O3 -fopenmp -I./libs/sdsl/include
# Add library path for the sdsl library to linker flags
LDFLAGS = -L./libs/sdsl/lib
# Specify the libraries to link against, including both divsufsort versions
LDLIBS = -lsdsl -ldivsufsort -ldivsufsort64 -lstdc++fs -lz -lrt -pthread

# Project structure and file definitions
BUILD_DIR = build
TARGET_EXEC = dp_cost
EXECUTABLE = $(BUILD_DIR)/$(TARGET_EXEC)

# Automatically find all .cpp files in the current directory
SRCS = $(wildcard *.cpp)
# Create a corresponding list of object files to be placed in the build directory
OBJS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRCS))

# SDSL submodule configuration
SDSL_DIR = libs/sdsl
SDSL_LIB_FILE = $(SDSL_DIR)/lib/libsdsl.a

# --- Makefile Targets ---

# The default target, called when you just run 'make'
.PHONY: all
all: $(EXECUTABLE)

# Rule to link the final executable file.
# After a successful link, the intermediate object files are removed.
$(EXECUTABLE): $(OBJS) | $(BUILD_DIR)
	@echo "--- Linking executable: $(EXECUTABLE) ---"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
	@echo "--- Removing intermediate object file(s) ---"
	@rm $(OBJS)

# A rule to create the build directory. This is an order-only prerequisite
# for rules that write into the build directory.
$(BUILD_DIR):
	@echo "--- Creating build directory if it does not exist ---"
	@mkdir -p $(BUILD_DIR)

# Pattern rule to compile .cpp source files into .o object files.
# This rule now has an order-only prerequisite on the sdsl library file and the build directory.
$(BUILD_DIR)/%.o: %.cpp | $(SDSL_LIB_FILE) $(BUILD_DIR)
	@echo "--- Compiling $< into $@ ---"
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to build the sdsl-lite dependency.
# NOTE: The sdsl-lite library is built as a single unit (`libsdsl.a`) using its own build system.
# It is not designed to be built piece-by-piece. Attempting to manually compile only the files
# related to csa_wt, construct, etc., would be complex and error-prone due to deep internal
# dependencies within the library.
# The good news is that this is a one-time build cost. More importantly, when linking your
# final executable, the linker will automatically include ONLY the parts of the library your
# code actually uses, keeping the final file size optimized.
$(SDSL_LIB_FILE):
	@echo "--- Building dependency: SDSL library ---"
	cd $(SDSL_DIR) && ./install.sh .

# Target to clean up the project's build files and the submodule's build files.
.PHONY: clean
clean:
	@echo "--- Cleaning project and SDSL build files ---"
	rm -rf $(BUILD_DIR)
	rm -rf $(SDSL_DIR)/build

