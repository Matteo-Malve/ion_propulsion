# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dealii/shared-folder/ion_propulsion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dealii/shared-folder/ion_propulsion

# Utility rule file for cleanmesh.

# Include any custom commands dependencies for this target.
include CMakeFiles/cleanmesh.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cleanmesh.dir/progress.make

CMakeFiles/cleanmesh:
	cd /home/dealii/shared-folder/mesh_storage && rm *.vtu

cleanmesh: CMakeFiles/cleanmesh
cleanmesh: CMakeFiles/cleanmesh.dir/build.make
.PHONY : cleanmesh

# Rule to build all files generated by this target.
CMakeFiles/cleanmesh.dir/build: cleanmesh
.PHONY : CMakeFiles/cleanmesh.dir/build

CMakeFiles/cleanmesh.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cleanmesh.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cleanmesh.dir/clean

CMakeFiles/cleanmesh.dir/depend:
	cd /home/dealii/shared-folder/ion_propulsion && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dealii/shared-folder/ion_propulsion /home/dealii/shared-folder/ion_propulsion /home/dealii/shared-folder/ion_propulsion /home/dealii/shared-folder/ion_propulsion /home/dealii/shared-folder/ion_propulsion/CMakeFiles/cleanmesh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cleanmesh.dir/depend

