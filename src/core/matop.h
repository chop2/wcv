#pragma once
#include "matrix.h"

namespace wcv {

	typedef Matrix_<uchar> Mat8u;
	typedef Matrix_<short> Mat16s;
	typedef Matrix_<unsigned short> Mat16u;
	typedef Matrix_<int>	Mat32s;
	typedef Matrix_<unsigned int> Mat32u;
	typedef Matrix_<float> Mat32f;
	typedef Matrix_<double> Mat64f;


	/*********************************generate matrix function******************************************/
	template<typename _Tp>
	/**@brief operate add*/
	void add(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst);

	template<typename _Tp>
	/**@brief operate add*/
	void add(const Matrix_<_Tp>& s1, const Scalar4d scalar, Matrix_<_Tp>& dst);

	template<typename _Tp>
	/**@brief operate subtract*/
	void sub(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst);

	template<typename _Tp>
	/**@brief operate subtract*/
	void sub(const Matrix_<_Tp>& s1, const Scalar4d scalar, Matrix_<_Tp>& dst);
	
	template<typename _Tp>
	Scalar4d sum(const Matrix_<_Tp>& s);

	template<typename _Tp>
	/**@brief make image extend by input kernel
	e.g. image by 7x7 boarder in memory:
	|***************************************************|
	|***************************************************|
	|***************************************************|
	|***xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx***|
	|***xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx***|
						...
	|***xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx***|
	|***************************************************|
	|***************************************************|
	|***************************************************|
	"x" means data in raw image
	"*" means boarder data
	@param src - input matrix data
	@param kSize - input kernel size
	@param type - boarder type,see eBoarderType
		"*" will be filled by zeros,if use NONE.	-> 000|123
		"*" will be same with raw data, if use SAME -> 123|123
		"*" will be filled mirror data,if use MIRROR -> 321|123
	@param dst - output image with boarder
	*/
	void copymakeBoarder(const Matrix_<_Tp>& src, Size4i kSize, eBoarderType type, Matrix_<_Tp>& dst);

	template<typename _Tp>
	void templateOp(const Matrix_<_Tp>& src, const Matrix_<_Tp>& kernel, Matrix_<_Tp>& dst, eBoarderType type);

	template<typename _Tp>
	void split(const Matrix_<_Tp>& src, std::vector<Matrix_<_Tp> >& mvs);

	template<typename _Tp>
	void merge(const std::vector<Matrix_<_Tp> >& mvs, Matrix_<_Tp>& dst);

	template<typename _Tp = double>
	Matrix_<double> getGaussianKernel1D(int size, double sigma);

	template<typename _Tp = double>
	Matrix_<double> getGaussianKernel2D(const Size4i& kSize, float sita);

	template<typename _Tp = double>
	Matrix_<double> getDefaultGaussianKernel2D_3x3();
};