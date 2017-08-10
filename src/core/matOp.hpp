#pragma once
#include "core_base.hpp"

namespace wcv {

	#define Matrix_ Matrix_
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
	void add(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.height == s2.height &&
			s1.width == s2.width &&
			s1.nchannels == s2.nchannels);
		if (dst.empty())
			dst.create(s1.height, s1.width, s1.nchannels);
		int i = 0;
		for (; i < s1.totalSizes() - 8;i += 8) {
			*(dst.data + i + 0) = *(s1.data + i + 0) + *(s2.data + i + 0);
			*(dst.data + i + 1) = *(s1.data + i + 1) + *(s2.data + i + 1);
			*(dst.data + i + 2) = *(s1.data + i + 2) + *(s2.data + i + 2);
			*(dst.data + i + 3) = *(s1.data + i + 3) + *(s2.data + i + 3);
			*(dst.data + i + 4) = *(s1.data + i + 4) + *(s2.data + i + 4);
			*(dst.data + i + 5) = *(s1.data + i + 5) + *(s2.data + i + 5);
			*(dst.data + i + 6) = *(s1.data + i + 6) + *(s2.data + i + 6);
			*(dst.data + i + 7) = *(s1.data + i + 7) + *(s2.data + i + 7);
		}

		for (; i < s1.totalSizes(); i++) {
			*(dst.data + i) = *(s1.data + i) + *(s2.data + i);
		}
	}

	template<typename _Tp>
	/**@brief operate add*/
	void add(const Matrix_<_Tp>& s1, const Scalar4d scalar, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.height == s2.height &&
			s1.width == s2.width &&
			s1.nchannels == s2.nchannels);
		if (dst.empty())
			dst.create(s1.height, s1.width, s1.nchannels);
		if (s1.nchannels == 1) {
			int i = 0;
			for (; i < s1.totalSizes(); i++) {
				*(dst.data + i + 0) = *(s1.data + i + 0) + scalar[0];
				*(dst.data + i + 1) = *(s1.data + i + 1) + scalar[0];
				*(dst.data + i + 2) = *(s1.data + i + 2) + scalar[0];
				*(dst.data + i + 3) = *(s1.data + i + 3) + scalar[0];
				*(dst.data + i + 4) = *(s1.data + i + 4) + scalar[0];
				*(dst.data + i + 5) = *(s1.data + i + 5) + scalar[0];
				*(dst.data + i + 6) = *(s1.data + i + 6) + scalar[0];
				*(dst.data + i + 7) = *(s1.data + i + 7) + scalar[0];
			}

			for (; i < s1.totalSizes(); i++) {
				*(dst.data + i) = *(s1.data + i) + *(s2.data + i);
			}
		} else if (s1.nchannels == 3) {
			for (size_t i = 0; i < s1.height; i++)	{
				for (size_t j = 0; j < s1.width; j++) {
					*(dst.data + i*dst.step() + j*dst.nchannels + 0) += scalar[0];
					*(dst.data + i*dst.step() + j*dst.nchannels + 1) += scalar[1];
					*(dst.data + i*dst.step() + j*dst.nchannels + 2) += scalar[2];
				}
			}
		}
	}

	template<typename _Tp>
	/**@brief operate subtract*/
	void sub(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.height == s2.height &&
			s1.width == s2.width &&
			s1.nchannels == s2.nchannels);
		if (dst.empty())
			dst.create(s1.height, s1.width, s1.nchannels);
		int i = 0;
		for (; i < s1.totalSizes() - 8; i += 8) {
			*(dst.data + i + 0) = *(s1.data + i + 0) - *(s2.data + i + 0);
			*(dst.data + i + 1) = *(s1.data + i + 1) - *(s2.data + i + 1);
			*(dst.data + i + 2) = *(s1.data + i + 2) - *(s2.data + i + 2);
			*(dst.data + i + 3) = *(s1.data + i + 3) - *(s2.data + i + 3);
			*(dst.data + i + 4) = *(s1.data + i + 4) - *(s2.data + i + 4);
			*(dst.data + i + 5) = *(s1.data + i + 5) - *(s2.data + i + 5);
			*(dst.data + i + 6) = *(s1.data + i + 6) - *(s2.data + i + 6);
			*(dst.data + i + 7) = *(s1.data + i + 7) - *(s2.data + i + 7);
		}

		for (; i < s1.totalSizes(); i++) {
			*(dst.data + i) = *(s1.data + i) - *(s2.data + i);
		}
	}

	template<typename _Tp>
	/**@brief operate subtract*/
	void sub(const Matrix_<_Tp>& s1, const Scalar4d scalar, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.height == s2.height &&
			s1.width == s2.width &&
			s1.nchannels == s2.nchannels);
		if (dst.empty())
			dst.create(s1.height, s1.width, s1.nchannels);
		if (s1.nchannels == 1) {
			int i = 0;
			for (; i < s1.totalSizes(); i++) {
				*(dst.data + i + 0) = *(s1.data + i + 0) - scalar[0];
				*(dst.data + i + 1) = *(s1.data + i + 1) - scalar[0];
				*(dst.data + i + 2) = *(s1.data + i + 2) - scalar[0];
				*(dst.data + i + 3) = *(s1.data + i + 3) - scalar[0];
				*(dst.data + i + 4) = *(s1.data + i + 4) - scalar[0];
				*(dst.data + i + 5) = *(s1.data + i + 5) - scalar[0];
				*(dst.data + i + 6) = *(s1.data + i + 6) - scalar[0];
				*(dst.data + i + 7) = *(s1.data + i + 7) - scalar[0];
			}

			for (; i < s1.totalSizes(); i++) {
				*(dst.data + i) = *(s1.data + i) - *(s2.data + i);
			}
		}
		else if (s1.nchannels == 3) {
			for (size_t i = 0; i < s1.height; i++) {
				for (size_t j = 0; j < s1.width; j++) {
					*(dst.data + i*dst.step() + j*dst.nchannels + 0) -= scalar[0];
					*(dst.data + i*dst.step() + j*dst.nchannels + 1) -= scalar[1];
					*(dst.data + i*dst.step() + j*dst.nchannels + 2) -= scalar[2];
				}
			}
		}
	}

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
	void copymakeBoarder(const Matrix_<_Tp>& src,Size4i kSize, eBoarderType type, Matrix_<_Tp>& dst) {
		assert(src.checkValid());
		assert(kSize.width % 2 != 0 && kSize.height % 2 != 0);
		int newheight = src.height + kSize.height - 1;
		int newwidth = src.width + kSize.width - 1;
		if (dst.empty()) {
			dst.create(newheight, newwidth, src.nchannels,0);
		}

		//! copy raw data
		_Tp* ptrSrc = NULL, *ptrDst = NULL;
		for (size_t i = kSize.height/2; i <= dst.height - kSize.height/2; i++) {
			ptrSrc = src.data + (i - kSize.height / 2)*src.step();
			ptrDst = dst.data + i*dst.step();
			int step = dst.step();
			size_t offsLeft = 0;
			size_t offsRight = (kSize.width / 2 * src.nchannels + src.step());
			size_t cpySize = kSize.width / 2 * src.nchannels;
			if (type == SAME) {	
				//left side
				memcpy(ptrDst, ptrSrc, cpySize);
				//mid side
				memcpy(ptrDst + cpySize, ptrSrc, step * sizeof(_Tp));
				//right side
				memcpy(ptrDst + offsRight, ptrSrc + src.step() - cpySize, cpySize);
			} else if (type == MIRROR) {
				_Tp* ptrD, *ptrS;
				ptrD = ptrDst;
				int cpys = cpySize / src.nchannels;
				//left side
				while (cpys--) {
					ptrS = ptrSrc + cpys*src.nchannels;
					for (size_t k = 0; k < src.nchannels; k++) {
						*(ptrD++) = *(ptrS++);
					}
				}
				//mid side
				memcpy(ptrDst + cpySize, ptrSrc, step * sizeof(_Tp));
				ptrD = ptrDst + offsRight;
				cpys = 0;
				//right side
				while (cpys++ < cpySize / src.nchannels) {
					ptrS = ptrSrc +(src.width - cpys)*src.nchannels;
					for (size_t k = 0; k < src.nchannels; k++) {
						*(ptrD++) = *(ptrS++);
					}
				}
			}
		}

		size_t cpySize = kSize.height / 2 * dst.step() * sizeof(_Tp);
		//! copy top-bottom boarder
		if (type == SAME) {
			_Tp* ptrDT = dst.data;
			_Tp* ptrST = dst.data + cpySize;
			_Tp* ptrDB = dst.data + dst.height * dst.step() - cpySize;
			_Tp* ptrSB = dst.data + (dst.height - kSize.height) * dst.step();
			memcpy(ptrDT, ptrST, cpySize);
			memcpy(ptrDB, ptrSB, cpySize);
		} else if (type == MIRROR) {
			for (size_t i = 0; i < kSize.height / 2; i++) {
				cpySize = dst.step() * sizeof(_Tp);
				_Tp* ptrDT = dst.data + i * dst.step();
				_Tp* ptrST = dst.data + ((kSize.height -1) - i - 1) * dst.step();
				_Tp* ptrDB = dst.data + (dst.height - 1 - i) * dst.step();
				_Tp* ptrSB = dst.data + (dst.height - (kSize.height - 1) + i) * dst.step();
				memcpy(ptrDT, ptrST, cpySize);
				memcpy(ptrDB, ptrSB, cpySize);
			}
		}
	}

	template<typename _Tp>
	void templateOp(const Matrix_<_Tp>& src,const Mat32f& kernel,Matrix_<_Tp>& dst,eBoarderType type) {
		assert(src.nchannels == 1 && kernel.nchannels == 1);
		Matrix_<_Tp> src_bd;
		copymakeBoarder(src, kernel.size(), type, src_bd);
		if (!dst.empty()) dst.release();
		dst.create(src.height, src.width, src.nchannels, 0);

		for (size_t i = kernel.width / 2; i < src_bd.height - kernel.width / 2; i++) {
			for (size_t j = kernel.width / 2; j < src_bd.width - kernel.width / 2; j++) {
				double sum = 0;
				for (int m = -kernel.height/2; m < kernel.height / 2; m++) {
					for (int n = -kernel.width / 2; n < kernel.width / 2; n++) {
						sum += (double)MAT_ELEM_S(src_bd, (i + m), (j + n)) * 
							(double)MAT_ELEM_S(kernel, (m + kernel.height / 2), (n + kernel.width / 2));
					}
				}
				sum /= kernel.size().area();
				int i0 = i - kernel.height / 2;
				int j0 = j - kernel.width / 2;
				MAT_ELEM_S(dst,i0,j0) = static_cast<_Tp>(sum);
			}
		}
	}
};