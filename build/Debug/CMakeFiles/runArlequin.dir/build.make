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
CMAKE_SOURCE_DIR = /home/patricia/Arlequin10

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/patricia/Arlequin10/build/Debug

# Include any dependencies generated for this target.
include CMakeFiles/runArlequin.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/runArlequin.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/runArlequin.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/runArlequin.dir/flags.make

CMakeFiles/runArlequin.dir/main.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/main.cpp.o: ../../main.cpp
CMakeFiles/runArlequin.dir/main.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/runArlequin.dir/main.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/main.cpp.o -MF CMakeFiles/runArlequin.dir/main.cpp.o.d -o CMakeFiles/runArlequin.dir/main.cpp.o -c /home/patricia/Arlequin10/main.cpp

CMakeFiles/runArlequin.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/main.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/main.cpp > CMakeFiles/runArlequin.dir/main.cpp.i

CMakeFiles/runArlequin.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/main.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/main.cpp -o CMakeFiles/runArlequin.dir/main.cpp.s

CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o: ../../src/Arlequin.cpp
CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o -MF CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o.d -o CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o -c /home/patricia/Arlequin10/src/Arlequin.cpp

CMakeFiles/runArlequin.dir/src/Arlequin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/Arlequin.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/Arlequin.cpp > CMakeFiles/runArlequin.dir/src/Arlequin.cpp.i

CMakeFiles/runArlequin.dir/src/Arlequin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/Arlequin.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/Arlequin.cpp -o CMakeFiles/runArlequin.dir/src/Arlequin.cpp.s

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o: ../../src/BoundaryIntegrationQuadrature.cpp
CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o -MF CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o.d -o CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o -c /home/patricia/Arlequin10/src/BoundaryIntegrationQuadrature.cpp

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/BoundaryIntegrationQuadrature.cpp > CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.i

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/BoundaryIntegrationQuadrature.cpp -o CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.s

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o: ../../src/BoundaryIntegrationQuadratureIso.cpp
CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o -MF CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o.d -o CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o -c /home/patricia/Arlequin10/src/BoundaryIntegrationQuadratureIso.cpp

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/BoundaryIntegrationQuadratureIso.cpp > CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.i

CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/BoundaryIntegrationQuadratureIso.cpp -o CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.s

CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o: ../../src/BoundaryShapeFunction.cpp
CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o -MF CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o.d -o CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o -c /home/patricia/Arlequin10/src/BoundaryShapeFunction.cpp

CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/BoundaryShapeFunction.cpp > CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.i

CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/BoundaryShapeFunction.cpp -o CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.s

CMakeFiles/runArlequin.dir/src/Element.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/Element.cpp.o: ../../src/Element.cpp
CMakeFiles/runArlequin.dir/src/Element.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/runArlequin.dir/src/Element.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/Element.cpp.o -MF CMakeFiles/runArlequin.dir/src/Element.cpp.o.d -o CMakeFiles/runArlequin.dir/src/Element.cpp.o -c /home/patricia/Arlequin10/src/Element.cpp

CMakeFiles/runArlequin.dir/src/Element.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/Element.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/Element.cpp > CMakeFiles/runArlequin.dir/src/Element.cpp.i

CMakeFiles/runArlequin.dir/src/Element.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/Element.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/Element.cpp -o CMakeFiles/runArlequin.dir/src/Element.cpp.s

CMakeFiles/runArlequin.dir/src/FluidData.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/FluidData.cpp.o: ../../src/FluidData.cpp
CMakeFiles/runArlequin.dir/src/FluidData.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/runArlequin.dir/src/FluidData.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/FluidData.cpp.o -MF CMakeFiles/runArlequin.dir/src/FluidData.cpp.o.d -o CMakeFiles/runArlequin.dir/src/FluidData.cpp.o -c /home/patricia/Arlequin10/src/FluidData.cpp

CMakeFiles/runArlequin.dir/src/FluidData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/FluidData.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/FluidData.cpp > CMakeFiles/runArlequin.dir/src/FluidData.cpp.i

CMakeFiles/runArlequin.dir/src/FluidData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/FluidData.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/FluidData.cpp -o CMakeFiles/runArlequin.dir/src/FluidData.cpp.s

CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o: ../../src/IntegrationQuadrature.cpp
CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o -MF CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o.d -o CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o -c /home/patricia/Arlequin10/src/IntegrationQuadrature.cpp

CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/IntegrationQuadrature.cpp > CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.i

CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/IntegrationQuadrature.cpp -o CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.s

CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o: ../../src/QuadraticShapeFunction.cpp
CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o -MF CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o.d -o CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o -c /home/patricia/Arlequin10/src/QuadraticShapeFunction.cpp

CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/QuadraticShapeFunction.cpp > CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.i

CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/QuadraticShapeFunction.cpp -o CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.s

CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/flags.make
CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o: ../../src/SpecialIntegrationQuadrature.cpp
CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o: CMakeFiles/runArlequin.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o -MF CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o.d -o CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o -c /home/patricia/Arlequin10/src/SpecialIntegrationQuadrature.cpp

CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.i"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/patricia/Arlequin10/src/SpecialIntegrationQuadrature.cpp > CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.i

CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.s"
	/usr/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/patricia/Arlequin10/src/SpecialIntegrationQuadrature.cpp -o CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.s

# Object files for target runArlequin
runArlequin_OBJECTS = \
"CMakeFiles/runArlequin.dir/main.cpp.o" \
"CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o" \
"CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o" \
"CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o" \
"CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o" \
"CMakeFiles/runArlequin.dir/src/Element.cpp.o" \
"CMakeFiles/runArlequin.dir/src/FluidData.cpp.o" \
"CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o" \
"CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o" \
"CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o"

# External object files for target runArlequin
runArlequin_EXTERNAL_OBJECTS =

runArlequin: CMakeFiles/runArlequin.dir/main.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/Arlequin.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadrature.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/BoundaryIntegrationQuadratureIso.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/BoundaryShapeFunction.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/Element.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/FluidData.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/IntegrationQuadrature.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/QuadraticShapeFunction.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/src/SpecialIntegrationQuadrature.cpp.o
runArlequin: CMakeFiles/runArlequin.dir/build.make
runArlequin: /home/patricia/mpich-3.4.2/install/lib/libmpicxx.so
runArlequin: /home/patricia/mpich-3.4.2/install/lib/libmpi.so
runArlequin: /home/patricia/petsc-3.18.0/arch-linux2-c-debug/lib/libpetsc.so
runArlequin: CMakeFiles/runArlequin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/patricia/Arlequin10/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable runArlequin"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/runArlequin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/runArlequin.dir/build: runArlequin
.PHONY : CMakeFiles/runArlequin.dir/build

CMakeFiles/runArlequin.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/runArlequin.dir/cmake_clean.cmake
.PHONY : CMakeFiles/runArlequin.dir/clean

CMakeFiles/runArlequin.dir/depend:
	cd /home/patricia/Arlequin10/build/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/patricia/Arlequin10 /home/patricia/Arlequin10 /home/patricia/Arlequin10/build/Debug /home/patricia/Arlequin10/build/Debug /home/patricia/Arlequin10/build/Debug/CMakeFiles/runArlequin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/runArlequin.dir/depend

