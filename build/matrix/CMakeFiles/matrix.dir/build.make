# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/public/cmake/2.8.9/bin/cmake

# The command to remove a file.
RM = /usr/public/cmake/2.8.9/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/public/cmake/2.8.9/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yeba/B_ITensor

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yeba/B_ITensor/build

# Include any dependencies generated for this target.
include matrix/CMakeFiles/matrix.dir/depend.make

# Include the progress variables for this target.
include matrix/CMakeFiles/matrix.dir/progress.make

# Include the compile flags for this target's objects.
include matrix/CMakeFiles/matrix.dir/flags.make

matrix/CMakeFiles/matrix.dir/matrix.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/matrix.cc.o: ../matrix/matrix.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/matrix.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/matrix.cc.o -c /home/yeba/B_ITensor/matrix/matrix.cc

matrix/CMakeFiles/matrix.dir/matrix.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/matrix.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/matrix.cc > CMakeFiles/matrix.dir/matrix.cc.i

matrix/CMakeFiles/matrix.dir/matrix.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/matrix.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/matrix.cc -o CMakeFiles/matrix.dir/matrix.cc.s

matrix/CMakeFiles/matrix.dir/matrix.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/matrix.cc.o.requires

matrix/CMakeFiles/matrix.dir/matrix.cc.o.provides: matrix/CMakeFiles/matrix.dir/matrix.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/matrix.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/matrix.cc.o.provides

matrix/CMakeFiles/matrix.dir/matrix.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/matrix.cc.o

matrix/CMakeFiles/matrix.dir/utility.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/utility.cc.o: ../matrix/utility.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/utility.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/utility.cc.o -c /home/yeba/B_ITensor/matrix/utility.cc

matrix/CMakeFiles/matrix.dir/utility.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/utility.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/utility.cc > CMakeFiles/matrix.dir/utility.cc.i

matrix/CMakeFiles/matrix.dir/utility.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/utility.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/utility.cc -o CMakeFiles/matrix.dir/utility.cc.s

matrix/CMakeFiles/matrix.dir/utility.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/utility.cc.o.requires

matrix/CMakeFiles/matrix.dir/utility.cc.o.provides: matrix/CMakeFiles/matrix.dir/utility.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/utility.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/utility.cc.o.provides

matrix/CMakeFiles/matrix.dir/utility.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/utility.cc.o

matrix/CMakeFiles/matrix.dir/sparse.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/sparse.cc.o: ../matrix/sparse.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/sparse.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/sparse.cc.o -c /home/yeba/B_ITensor/matrix/sparse.cc

matrix/CMakeFiles/matrix.dir/sparse.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/sparse.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/sparse.cc > CMakeFiles/matrix.dir/sparse.cc.i

matrix/CMakeFiles/matrix.dir/sparse.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/sparse.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/sparse.cc -o CMakeFiles/matrix.dir/sparse.cc.s

matrix/CMakeFiles/matrix.dir/sparse.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/sparse.cc.o.requires

matrix/CMakeFiles/matrix.dir/sparse.cc.o.provides: matrix/CMakeFiles/matrix.dir/sparse.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/sparse.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/sparse.cc.o.provides

matrix/CMakeFiles/matrix.dir/sparse.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/sparse.cc.o

matrix/CMakeFiles/matrix.dir/david.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/david.cc.o: ../matrix/david.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/david.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/david.cc.o -c /home/yeba/B_ITensor/matrix/david.cc

matrix/CMakeFiles/matrix.dir/david.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/david.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/david.cc > CMakeFiles/matrix.dir/david.cc.i

matrix/CMakeFiles/matrix.dir/david.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/david.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/david.cc -o CMakeFiles/matrix.dir/david.cc.s

matrix/CMakeFiles/matrix.dir/david.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/david.cc.o.requires

matrix/CMakeFiles/matrix.dir/david.cc.o.provides: matrix/CMakeFiles/matrix.dir/david.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/david.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/david.cc.o.provides

matrix/CMakeFiles/matrix.dir/david.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/david.cc.o

matrix/CMakeFiles/matrix.dir/sparseref.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/sparseref.cc.o: ../matrix/sparseref.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/sparseref.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/sparseref.cc.o -c /home/yeba/B_ITensor/matrix/sparseref.cc

matrix/CMakeFiles/matrix.dir/sparseref.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/sparseref.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/sparseref.cc > CMakeFiles/matrix.dir/sparseref.cc.i

matrix/CMakeFiles/matrix.dir/sparseref.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/sparseref.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/sparseref.cc -o CMakeFiles/matrix.dir/sparseref.cc.s

matrix/CMakeFiles/matrix.dir/sparseref.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/sparseref.cc.o.requires

matrix/CMakeFiles/matrix.dir/sparseref.cc.o.provides: matrix/CMakeFiles/matrix.dir/sparseref.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/sparseref.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/sparseref.cc.o.provides

matrix/CMakeFiles/matrix.dir/sparseref.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/sparseref.cc.o

matrix/CMakeFiles/matrix.dir/hpsortir.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/hpsortir.cc.o: ../matrix/hpsortir.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/hpsortir.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/hpsortir.cc.o -c /home/yeba/B_ITensor/matrix/hpsortir.cc

matrix/CMakeFiles/matrix.dir/hpsortir.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/hpsortir.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/hpsortir.cc > CMakeFiles/matrix.dir/hpsortir.cc.i

matrix/CMakeFiles/matrix.dir/hpsortir.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/hpsortir.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/hpsortir.cc -o CMakeFiles/matrix.dir/hpsortir.cc.s

matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.requires

matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.provides: matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.provides

matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/hpsortir.cc.o

matrix/CMakeFiles/matrix.dir/daxpy.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/daxpy.cc.o: ../matrix/daxpy.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/daxpy.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/daxpy.cc.o -c /home/yeba/B_ITensor/matrix/daxpy.cc

matrix/CMakeFiles/matrix.dir/daxpy.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/daxpy.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/daxpy.cc > CMakeFiles/matrix.dir/daxpy.cc.i

matrix/CMakeFiles/matrix.dir/daxpy.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/daxpy.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/daxpy.cc -o CMakeFiles/matrix.dir/daxpy.cc.s

matrix/CMakeFiles/matrix.dir/daxpy.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/daxpy.cc.o.requires

matrix/CMakeFiles/matrix.dir/daxpy.cc.o.provides: matrix/CMakeFiles/matrix.dir/daxpy.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/daxpy.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/daxpy.cc.o.provides

matrix/CMakeFiles/matrix.dir/daxpy.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/daxpy.cc.o

matrix/CMakeFiles/matrix.dir/matrixref.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/matrixref.cc.o: ../matrix/matrixref.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/matrixref.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/matrixref.cc.o -c /home/yeba/B_ITensor/matrix/matrixref.cc

matrix/CMakeFiles/matrix.dir/matrixref.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/matrixref.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/matrixref.cc > CMakeFiles/matrix.dir/matrixref.cc.i

matrix/CMakeFiles/matrix.dir/matrixref.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/matrixref.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/matrixref.cc -o CMakeFiles/matrix.dir/matrixref.cc.s

matrix/CMakeFiles/matrix.dir/matrixref.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/matrixref.cc.o.requires

matrix/CMakeFiles/matrix.dir/matrixref.cc.o.provides: matrix/CMakeFiles/matrix.dir/matrixref.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/matrixref.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/matrixref.cc.o.provides

matrix/CMakeFiles/matrix.dir/matrixref.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/matrixref.cc.o

matrix/CMakeFiles/matrix.dir/storelink.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/storelink.cc.o: ../matrix/storelink.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/storelink.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/storelink.cc.o -c /home/yeba/B_ITensor/matrix/storelink.cc

matrix/CMakeFiles/matrix.dir/storelink.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/storelink.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/storelink.cc > CMakeFiles/matrix.dir/storelink.cc.i

matrix/CMakeFiles/matrix.dir/storelink.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/storelink.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/storelink.cc -o CMakeFiles/matrix.dir/storelink.cc.s

matrix/CMakeFiles/matrix.dir/storelink.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/storelink.cc.o.requires

matrix/CMakeFiles/matrix.dir/storelink.cc.o.provides: matrix/CMakeFiles/matrix.dir/storelink.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/storelink.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/storelink.cc.o.provides

matrix/CMakeFiles/matrix.dir/storelink.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/storelink.cc.o

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o: ../matrix/conjugate_gradient.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/conjugate_gradient.cc.o -c /home/yeba/B_ITensor/matrix/conjugate_gradient.cc

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/conjugate_gradient.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/conjugate_gradient.cc > CMakeFiles/matrix.dir/conjugate_gradient.cc.i

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/conjugate_gradient.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/conjugate_gradient.cc -o CMakeFiles/matrix.dir/conjugate_gradient.cc.s

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.requires

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.provides: matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.provides

matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o

matrix/CMakeFiles/matrix.dir/dgemm.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/dgemm.cc.o: ../matrix/dgemm.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/dgemm.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/dgemm.cc.o -c /home/yeba/B_ITensor/matrix/dgemm.cc

matrix/CMakeFiles/matrix.dir/dgemm.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/dgemm.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/dgemm.cc > CMakeFiles/matrix.dir/dgemm.cc.i

matrix/CMakeFiles/matrix.dir/dgemm.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/dgemm.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/dgemm.cc -o CMakeFiles/matrix.dir/dgemm.cc.s

matrix/CMakeFiles/matrix.dir/dgemm.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/dgemm.cc.o.requires

matrix/CMakeFiles/matrix.dir/dgemm.cc.o.provides: matrix/CMakeFiles/matrix.dir/dgemm.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/dgemm.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/dgemm.cc.o.provides

matrix/CMakeFiles/matrix.dir/dgemm.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/dgemm.cc.o

matrix/CMakeFiles/matrix.dir/svd.cc.o: matrix/CMakeFiles/matrix.dir/flags.make
matrix/CMakeFiles/matrix.dir/svd.cc.o: ../matrix/svd.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/yeba/B_ITensor/build/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object matrix/CMakeFiles/matrix.dir/svd.cc.o"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix.dir/svd.cc.o -c /home/yeba/B_ITensor/matrix/svd.cc

matrix/CMakeFiles/matrix.dir/svd.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix.dir/svd.cc.i"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/yeba/B_ITensor/matrix/svd.cc > CMakeFiles/matrix.dir/svd.cc.i

matrix/CMakeFiles/matrix.dir/svd.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix.dir/svd.cc.s"
	cd /home/yeba/B_ITensor/build/matrix && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/yeba/B_ITensor/matrix/svd.cc -o CMakeFiles/matrix.dir/svd.cc.s

matrix/CMakeFiles/matrix.dir/svd.cc.o.requires:
.PHONY : matrix/CMakeFiles/matrix.dir/svd.cc.o.requires

matrix/CMakeFiles/matrix.dir/svd.cc.o.provides: matrix/CMakeFiles/matrix.dir/svd.cc.o.requires
	$(MAKE) -f matrix/CMakeFiles/matrix.dir/build.make matrix/CMakeFiles/matrix.dir/svd.cc.o.provides.build
.PHONY : matrix/CMakeFiles/matrix.dir/svd.cc.o.provides

matrix/CMakeFiles/matrix.dir/svd.cc.o.provides.build: matrix/CMakeFiles/matrix.dir/svd.cc.o

# Object files for target matrix
matrix_OBJECTS = \
"CMakeFiles/matrix.dir/matrix.cc.o" \
"CMakeFiles/matrix.dir/utility.cc.o" \
"CMakeFiles/matrix.dir/sparse.cc.o" \
"CMakeFiles/matrix.dir/david.cc.o" \
"CMakeFiles/matrix.dir/sparseref.cc.o" \
"CMakeFiles/matrix.dir/hpsortir.cc.o" \
"CMakeFiles/matrix.dir/daxpy.cc.o" \
"CMakeFiles/matrix.dir/matrixref.cc.o" \
"CMakeFiles/matrix.dir/storelink.cc.o" \
"CMakeFiles/matrix.dir/conjugate_gradient.cc.o" \
"CMakeFiles/matrix.dir/dgemm.cc.o" \
"CMakeFiles/matrix.dir/svd.cc.o"

# External object files for target matrix
matrix_EXTERNAL_OBJECTS =

matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/matrix.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/utility.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/sparse.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/david.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/sparseref.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/hpsortir.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/daxpy.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/matrixref.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/storelink.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/dgemm.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/svd.cc.o
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/build.make
matrix/libmatrix.a: matrix/CMakeFiles/matrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libmatrix.a"
	cd /home/yeba/B_ITensor/build/matrix && $(CMAKE_COMMAND) -P CMakeFiles/matrix.dir/cmake_clean_target.cmake
	cd /home/yeba/B_ITensor/build/matrix && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
matrix/CMakeFiles/matrix.dir/build: matrix/libmatrix.a
.PHONY : matrix/CMakeFiles/matrix.dir/build

matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/matrix.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/utility.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/sparse.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/david.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/sparseref.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/hpsortir.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/daxpy.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/matrixref.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/storelink.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/conjugate_gradient.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/dgemm.cc.o.requires
matrix/CMakeFiles/matrix.dir/requires: matrix/CMakeFiles/matrix.dir/svd.cc.o.requires
.PHONY : matrix/CMakeFiles/matrix.dir/requires

matrix/CMakeFiles/matrix.dir/clean:
	cd /home/yeba/B_ITensor/build/matrix && $(CMAKE_COMMAND) -P CMakeFiles/matrix.dir/cmake_clean.cmake
.PHONY : matrix/CMakeFiles/matrix.dir/clean

matrix/CMakeFiles/matrix.dir/depend:
	cd /home/yeba/B_ITensor/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yeba/B_ITensor /home/yeba/B_ITensor/matrix /home/yeba/B_ITensor/build /home/yeba/B_ITensor/build/matrix /home/yeba/B_ITensor/build/matrix/CMakeFiles/matrix.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : matrix/CMakeFiles/matrix.dir/depend
