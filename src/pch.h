#pragma once
#include <vector>
#include <string>
using namespace std;

#ifdef HAVE_OPENCV
#	include <opencv2/core/core.hpp>
#endif

#if defined _MSC_VER
#	include<emmintrin.h> 
#	define USE_SSE
#endif

#ifndef _OPENMP
#	include<omp.h>
#	define USE_OMP
#endif