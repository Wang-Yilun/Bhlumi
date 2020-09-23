cmake_minimum_required (VERSION 2.8.11)
project (BHLumi)
enable_testing()

message("######################################################")
message("#                                                    #")
message("#       BHLUMI software V1.0                         #")
message("#       S.Jadach                                     #")
message("#                                                    #")
message("#       Documentation:                               #")
message("#       https://github.com/KrakowHEPSoft/BHLUMI      #")
message("#       GPL 3 license                                #")
message("#                                                    #")
message("######################################################")



list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
#set(CMAKE_MACOSX_RPATH TRUE) # use @rpath keyword for the macOS rpath

###########################
#
#    Compile seetings
#
###########################
# Add -fPIC for 64 bit systems
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()



###########################
#
#    Packages search
#
###########################
# looking for yaml-cpp
find_package(Yaml-cpp)
# looking for boost
find_package(Boost REQUIRED)
# GSL
find_package(GSL REQUIRED)
# ROOT:
#find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})
include_directories(${CMAKE_SOURCE_DIR}/include)

###########################
#
#    Libraries
#
###########################
file(GLOB BHLUMI_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc)
add_library(BHLUMI_static ${BHLUMI_SOURCES})
add_library(BHLUMI_shared SHARED ${BHLUMI_SOURCES})
set_target_properties(BHLUMI_static PROPERTIES OUTPUT_NAME HEPLike)
set_target_properties(BHLUMI_shared PROPERTIES OUTPUT_NAME HEPLike)
#target_link_libraries(BHLUMI_shared PRIVATE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
target_link_libraries(BHLUMI_shared PRIVATE ${Boost_LIBRARIES} ${GSL_LIBRARIES} )   

###########################
#
#    Executables
#
###########################

if(EXECUTABLE)
  

  add_executable(Example main/Example.cc)
  target_link_libraries(Example BHLUMI_shared yaml-cpp  ${Boost_LIBRARIES} )    
  add_test(BLExample Example)


  
endif()
