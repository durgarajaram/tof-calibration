cmake_minimum_required(VERSION 2.0)
project(tofCalib)
include_directories(. "$ENV{MAUS_ROOT_DIR}/third_party/install/include/python2.7/"
                      "$ENV{ROOTSYS}/include/" 
                      "$ENV{MAUS_ROOT_DIR}/third_party/install/include/"
                      "$ENV{MAUS_ROOT_DIR}/"
                      "$ENV{MAUS_ROOT_DIR}/src/legacy/"
                      "$ENV{MAUS_ROOT_DIR}/src/common_cpp/")

link_directories(. "$ENV{ROOTSYS}/lib"
                   "$ENV{MAUS_ROOT_DIR}/third_party/install/lib/"
                   "$ENV{MAUS_ROOT_DIR}/build")

set(CMAKE_CXX_FLAGS "-D_STEPIV_DATA")
add_executable (tofCalib calib.cpp
                         MiceTofCalib.cc
                         TofTriggerCalib.cc
                         TofT0Calib.cc
                         TofTWCalib.cc)

target_link_libraries(tofCalib Core 
                               Cint 
                               RIO 
                               Net 
                               Hist 
                               Graf 
                               Graf3d 
                               Gpad 
                               Tree 
                               Rint 
                               Postscript 
                               Matrix 
                               Physics 
                               MathCore 
                               Thread 
                               Gui 
                               Minuit 
                               Spectrum
                               CLHEP
                               MausCpp)

