# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.10.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.10.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build

# Include any dependencies generated for this target.
include CMakeFiles/convolve_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/convolve_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/convolve_test.dir/flags.make

CMakeFiles/convolve_test.dir/convolve_test.cxx.o: CMakeFiles/convolve_test.dir/flags.make
CMakeFiles/convolve_test.dir/convolve_test.cxx.o: ../convolve_test.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/convolve_test.dir/convolve_test.cxx.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/convolve_test.dir/convolve_test.cxx.o -c /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/convolve_test.cxx

CMakeFiles/convolve_test.dir/convolve_test.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convolve_test.dir/convolve_test.cxx.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/convolve_test.cxx > CMakeFiles/convolve_test.dir/convolve_test.cxx.i

CMakeFiles/convolve_test.dir/convolve_test.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convolve_test.dir/convolve_test.cxx.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/convolve_test.cxx -o CMakeFiles/convolve_test.dir/convolve_test.cxx.s

CMakeFiles/convolve_test.dir/convolve_test.cxx.o.requires:

.PHONY : CMakeFiles/convolve_test.dir/convolve_test.cxx.o.requires

CMakeFiles/convolve_test.dir/convolve_test.cxx.o.provides: CMakeFiles/convolve_test.dir/convolve_test.cxx.o.requires
	$(MAKE) -f CMakeFiles/convolve_test.dir/build.make CMakeFiles/convolve_test.dir/convolve_test.cxx.o.provides.build
.PHONY : CMakeFiles/convolve_test.dir/convolve_test.cxx.o.provides

CMakeFiles/convolve_test.dir/convolve_test.cxx.o.provides.build: CMakeFiles/convolve_test.dir/convolve_test.cxx.o


CMakeFiles/convolve_test.dir/opencl_utils.cxx.o: CMakeFiles/convolve_test.dir/flags.make
CMakeFiles/convolve_test.dir/opencl_utils.cxx.o: ../opencl_utils.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/convolve_test.dir/opencl_utils.cxx.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/convolve_test.dir/opencl_utils.cxx.o -c /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/opencl_utils.cxx

CMakeFiles/convolve_test.dir/opencl_utils.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convolve_test.dir/opencl_utils.cxx.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/opencl_utils.cxx > CMakeFiles/convolve_test.dir/opencl_utils.cxx.i

CMakeFiles/convolve_test.dir/opencl_utils.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convolve_test.dir/opencl_utils.cxx.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/opencl_utils.cxx -o CMakeFiles/convolve_test.dir/opencl_utils.cxx.s

CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.requires:

.PHONY : CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.requires

CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.provides: CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.requires
	$(MAKE) -f CMakeFiles/convolve_test.dir/build.make CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.provides.build
.PHONY : CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.provides

CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.provides.build: CMakeFiles/convolve_test.dir/opencl_utils.cxx.o


# Object files for target convolve_test
convolve_test_OBJECTS = \
"CMakeFiles/convolve_test.dir/convolve_test.cxx.o" \
"CMakeFiles/convolve_test.dir/opencl_utils.cxx.o"

# External object files for target convolve_test
convolve_test_EXTERNAL_OBJECTS =

convolve_test: CMakeFiles/convolve_test.dir/convolve_test.cxx.o
convolve_test: CMakeFiles/convolve_test.dir/opencl_utils.cxx.o
convolve_test: CMakeFiles/convolve_test.dir/build.make
convolve_test: CMakeFiles/convolve_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable convolve_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/convolve_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/convolve_test.dir/build: convolve_test

.PHONY : CMakeFiles/convolve_test.dir/build

CMakeFiles/convolve_test.dir/requires: CMakeFiles/convolve_test.dir/convolve_test.cxx.o.requires
CMakeFiles/convolve_test.dir/requires: CMakeFiles/convolve_test.dir/opencl_utils.cxx.o.requires

.PHONY : CMakeFiles/convolve_test.dir/requires

CMakeFiles/convolve_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/convolve_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/convolve_test.dir/clean

CMakeFiles/convolve_test.dir/depend:
	cd /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build /Users/tdesell/Dropbox/webpages/home_page/files/opencl_matrix_add/build/CMakeFiles/convolve_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/convolve_test.dir/depend

