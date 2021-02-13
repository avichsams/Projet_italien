# ---
# CMake implementation for TP program
# Source set[compiler] scripts to force a specific compiler.
# Luca Sciacovelli (luca.sciacovelli@ensam.eu)
# ---

cmake_minimum_required(VERSION 3.0)
project(program)
enable_language(Fortran)

# default installation
get_filename_component(default_prefix ".." ABSOLUTE)
set( CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
     "Choose the installation directory; by default it installs in the main directory."
     FORCE)

# Find the modules included with HPCFD
set( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake ${CMAKE_MODULE_PATH} )

# ------------------------------------------------------------------------------
# Setup compiler flags
# Optional flags

# set( CMAKE_VERBOSE_MAKEFILE ON ) # Verbose CMAKE

# force debug flag in debug build
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message( "-- Build type is " ${CMAKE_BUILD_TYPE} " ...forcing OPT_VERBOSE ON" )
  set( OPT_VERBOSE ON CACHE BOOL "docstring" FORCE)
endif( )

# Check target optional flags conflicts
list( APPEND TARGET_OPTIONS ${FLAG_SANDY-IVY} ${FLAG_HASWELL-BROADWELL} ${FLAG_SKYLAKE-COFFELAKE} )
set( index 0 )
foreach( out_loop ${TARGET_OPTIONS} )
  if( ${out_loop} )
    # eliminate self-conflict
    list( REMOVE_AT TARGET_OPTIONS "${index}" )
    # check all the others
    foreach( in_loop ${TARGET_OPTIONS} )
      if ( ${in_loop} )
        message(FATAL_ERROR "Can't have optimization flags ON for different targets at same time!")
      endif( )
    endforeach( )
  endif( )
  # counter
  MATH(EXPR index "${index}+1")
endforeach( )

# echo platform
message ("-- Detected platform is: " ${CMAKE_SYSTEM_NAME} ", " ${CMAKE_SYSTEM})

# ------------------------------------------------------------------------------
# INTEL
# ------------------------------------------------------------------------------

# FFLAGS depend on the compiler and or platform
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -r8 -w -fp-model precise -xhost")

   # ----------------------------------------------------------------------------
   set (CMAKE_Fortran_FLAGS_DEBUG  "-O0 -g -r8 -fpp -C -fpe0 -init=snan \
                                    -debug all -gen-interface -traceback \
                                    -check arg_temp_created -check bounds \
                                    -warn all -fstack-protector-all -diag-disable 8889")
   # Apple system specific flags for ifort
   if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
       set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xhost ")
       set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -no_pie")
   endif ()

# ------------------------------------------------------------------------------
# GNU
# ------------------------------------------------------------------------------
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fdefault-real-8 -fdefault-double-8 \
                                    -cpp -funroll-all-loops -march=native -ffree-line-length-512 ")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fdefault-real-8 -fdefault-double-8 -cpp \
                                    -fcheck=all -Wall -Wextra  -ffree-line-length-512 \
                                    -ffpe-trap=invalid,zero,overflow,underflow,denormal ")
else ()
  message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message(FATAL_ERROR "This compiler is unknown, it needs to be added in CMakeLists.txt")
endif ()

# set executable name
set( EXE_NAME program.x )

# generate source list
set( EXE_SRC mod_precision.f90
             mod_mode.f90
             mod_filter.f90
             mod_deriv.f90
             mod_grid.f90
             mod_fluid.f90
             mod_work.f90
             mod_time.f90
             compute_lad.f90
             compute_thermo.f90
             derivative.f90
             deriv_coeff.f90
             deriv_04_exp.f90
             deriv_10_exp.f90
             deltat_fromCFL.f90
             filter_coeff.f90
             filter_exp.f90
             filter_shock_jameson.f90
             filter_shock_bogey.f90
             rhs_inviscid.f90
             rhs_viscous.f90
             runge_kutta4.f90
             runge_kutta6.f90
             runge_kutta2.f90
             setup_init_cond.f90
             update_dirichbc.f90
             update_neumannbc.f90
             update_ghost.f90
             update_var.f90
             write_solution.f90
             program.f90
             )

set( EXECUTABLES ${EXECUTABLES} "${EXE_NAME}" )
add_executable( ${EXE_NAME} ${EXE_SRC} )
install( TARGETS ${EXECUTABLES} RUNTIME DESTINATION "bin" )