#pragma once
#include "matrix.h"
#include "matop.h"
#include "matrix.inl.hpp"
#include "matop.inl.hpp"

namespace wcv {
	typedef _Scalar<int> Scalar4i;
	typedef _Scalar<float> Scalar4f;
	typedef _Scalar<double> Scalar4d;

	typedef _Size<int> Size4i;
	typedef _Size<float> Size4f;
	typedef _Size<double> Size4d;

	typedef Range_<int> Range4i;
	typedef Range_<float> Range4f;
	typedef Range_<double> Range4d;

	typedef Matrix_<uchar> Image8u;
	typedef Image8u Image;
	typedef Matrix_<short> Image16s;
	typedef Matrix_<unsigned short> Image16u;
	typedef Matrix_<int> Image32s;
	typedef Matrix_<unsigned int> Image32u;
	typedef Matrix_<float> Image32f;
	typedef Matrix_<double> Image64f;
}