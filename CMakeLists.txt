# FasTC
# Copyright (c) 2012 University of North Carolina at Chapel Hill. All rights reserved.
#
# Permission to use, copy, modify, and distribute this software and its documentation for educational, 
# research, and non-profit purposes, without fee, and without a written agreement is hereby granted, 
# provided that the above copyright notice, this paragraph, and the following four paragraphs appear 
# in all copies.
#
# Permission to incorporate this software into commercial products may be obtained by contacting the 
# authors or the Office of Technology Development at the University of North Carolina at Chapel Hill <otd@unc.edu>.
#
# This software program and documentation are copyrighted by the University of North Carolina at Chapel Hill. 
# The software program and documentation are supplied "as is," without any accompanying services from the 
# University of North Carolina at Chapel Hill or the authors. The University of North Carolina at Chapel Hill 
# and the authors do not warrant that the operation of the program will be uninterrupted or error-free. The 
# end-user understands that the program was developed for research purposes and is advised not to rely 
# exclusively on the program for any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL OR THE AUTHORS BE LIABLE TO ANY PARTY FOR 
# DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE 
# USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL OR THE 
# AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL AND THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, 
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY 
# STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY 
# OF NORTH CAROLINA AT CHAPEL HILL AND THE AUTHORS HAVE NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
# ENHANCEMENTS, OR MODIFICATIONS.
#
# Please send all BUG REPORTS to <pavel@cs.unc.edu>.
#
# The authors may be contacted via:
#
# Pavel Krajcevski
# Dept of Computer Science
# 201 S Columbia St
# Frederick P. Brooks, Jr. Computer Science Bldg
# Chapel Hill, NC 27599-3175
# USA
# 
# <http://gamma.cs.unc.edu/FasTC/>

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
PROJECT(GenTC)

OPTION(TREAT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors. We use the highest warnings levels for compilers." OFF)

IF(MSVC)
	SET(MSVC_INSTALL_PATH "${PROJECT_SOURCE_DIR}/Windows")
	SET(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${MSVC_INSTALL_PATH}")

	IF(MSVC10)
		SET(MSVC_VERSION_STRING vc100)
	ELSEIF(MSVC11)
		SET(MSVC_VERSION_STRING vc110)
	ELSEIF(MSVC90)
		SET(MSVC_VERSION_STRING vc90)
	ELSEIF(MSVC80)
		SET(MSVC_VERSION_STRING vc80)
	ENDIF()

	# !FIXME! Actually detect compiler architecture version....
	IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
		SET(MSVC_ARCHITECTURE_STRING x64)
	ELSE()
		SET(MSVC_ARCHITECTURE_STRING x86)
	ENDIF()

	SET(MSVC_LIB_DIR "${MSVC_INSTALL_PATH}/lib/${MSVC_ARCHITECTURE_STRING}/${MSVC_VERSION_STRING}")
	SET(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH};${MSVC_LIB_DIR}")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

ELSEIF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fms-extensions")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fms-extensions")

ENDIF(MSVC)

IF(TREAT_WARNINGS_AS_ERRORS)
  IF(MSVC)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /WX")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX")
  ELSEIF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  ENDIF(MSVC)
ENDIF(TREAT_WARNINGS_AS_ERRORS)

SET(FASTC_DIRECTORY "" CACHE FILEPATH "Path to the FasTC directory")

IF(FASTC_DIRECTORY STREQUAL "")
  MESSAGE(FATAL_ERROR "FasTC not found.")
ELSE()
  SET(FASTC_BINDIR "${CMAKE_CURRENT_BINARY_DIR}/build_FasTC/")
  FILE(MAKE_DIRECTORY ${FASTC_BINDIR})
  ADD_SUBDIRECTORY(${FASTC_DIRECTORY} ${FASTC_BINDIR})
ENDIF()

INCLUDE_DIRECTORIES(${FASTC_DIRECTORY}/BPTCEncoder/include)
INCLUDE_DIRECTORIES(${FASTC_BINDIR}/BPTCEncoder/include)
INCLUDE_DIRECTORIES(${FASTC_DIRECTORY}/Base/include)
INCLUDE_DIRECTORIES(${FASTC_DIRECTORY}/Core/include)
INCLUDE_DIRECTORIES(${FASTC_DIRECTORY}/IO/include)
INCLUDE_DIRECTORIES(${FASTC_BINDIR}/Base/include)
INCLUDE_DIRECTORIES(${FASTC_BINDIR}/IO/include)

SET(SOURCES
  "main.cpp"
  "SLIC.cpp"
  "Partition.cpp")

SET(HEADERS
  "SLIC.h"
  "Partition.h"
  "VPTree.h")

ADD_EXECUTABLE(sc ${SOURCES} ${HEADERS})
ADD_EXECUTABLE(vptree_test "VPTreeTest.cpp" "VPTree.h")

IF( MSVC )
  SET_TARGET_PROPERTIES(sc PROPERTIES LINK_FLAGS "/LTCG")
ENDIF()

TARGET_LINK_LIBRARIES( sc FasTCBase )
TARGET_LINK_LIBRARIES( sc FasTCIO )
TARGET_LINK_LIBRARIES( sc FasTCCore )

TARGET_LINK_LIBRARIES( vptree_test FasTCCore )
