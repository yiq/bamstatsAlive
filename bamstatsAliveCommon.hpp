/**
 * @file bamstatsAliveCommon.hpp
 * This is the common header for the entire project.
 *
 * This header file will be compiled into precompiled header, and
 * automatically be included in all source code's compilation.
 *
 * @author Yi Qiao
 */

#ifndef BAMSTATSALIVECOMMON_H
#define BAMSTATSALIVECOMMON_H

// BamTools
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <src/utils/bamtools_pileup_engine.h>

// Standard C/C++ Headers
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>

// Commonly used STL classes
#include <string>
#include <vector>
#include <map>
#include <algorithm>

// Utilize assertions
#include <assert.h>

// Jansson JSON manipulation library
#include <jansson.h>


// Logging facility for Debug
#ifdef RELEASE
namespace BamstatsAlive {
	namespace Log {
		static std::ostream bitBucket(0);
	}
}
#define LOGS (Log::bitBucket)
#endif

#ifdef DEBUG
#define LOGS (std::cerr)
#endif

#endif
