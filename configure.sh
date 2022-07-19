#!/usr/bin/env bash

set -e # Instructs the shell to exit if a command fails.

source env.sh # source this from your run script too
              # determine which environment directory is to be used for
              # the current host (i.e. allows to load specific modules
              # on ollie)
              
mkdir build || true # make sure not to commit this to svn or git
                    # compiled code will be saved to build directory
cd build
cmake .. # generates a build system
         # .. is the source tree (directory which must contain
         # CMakeLists.txt file)
		 # not required when re-compiling (if build tree with
		 # CMakeCache.txt file exists)
		 
make install -j`nproc --all` # -j: set number of CPU cores/threads used
                             # nproc --all: counts number of available
                             # CPU cores/threads
