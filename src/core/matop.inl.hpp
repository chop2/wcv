#pragma once
#include "matop.h"

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
	void add(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.rows == s2.rows &&
			s1.cols == s2.cols &&
			s1.channels == s2.channels);
		if (dst.empty())
			dst.create(s1.rows, s1.cols, s1.channels);
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
		assert(s1.rows == s2.rows &&
			s1.cols == s2.cols &&
			s1.channels == s2.channels);
		if (dst.empty())
			dst.create(s1.rows, s1.cols, s1.channels);
		if (s1.channels == 1) {
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
		} else if (s1.channels == 3) {
			for (size_t i = 0; i < s1.rows; i++)	{
				for (size_t j = 0; j < s1.cols; j++) {
					*(dst.data + i*dst.step() + j*dst.channels + 0) += scalar[0];
					*(dst.data + i*dst.step() + j*dst.channels + 1) += scalar[1];
					*(dst.data + i*dst.step() + j*dst.channels + 2) += scalar[2];
				}
			}
		}
	}

	template<typename _Tp>
	/**@brief operate subtract*/
	void sub(const Matrix_<_Tp>& s1, const Matrix_<_Tp>& s2, Matrix_<_Tp>& dst) {
		assert(s1.checkValid() && s2.checkValid());
		assert(s1.rows == s2.rows &&
			s1.cols == s2.cols &&
			s1.channels == s2.channels);
		if (dst.empty())
			dst.create(s1.rows, s1.cols, s1.channels);
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
		assert(s1.rows == s2.rows &&
			s1.cols == s2.cols &&
			s1.channels == s2.channels);
		if (dst.empty())
			dst.create(s1.rows, s1.cols, s1.channels);
		if (s1.channels == 1) {
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
		else if (s1.channels == 3) {
			for (size_t i = 0; i < s1.rows; i++) {
				for (size_t j = 0; j < s1.cols; j++) {
					*(dst.data + i*dst.step() + j*dst.channels + 0) -= scalar[0];
					*(dst.data + i*dst.step() + j*dst.channels + 1) -= scalar[1];
					*(dst.data + i*dst.step() + j*dst.channels + 2) -= scalar[2];
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
		assert(kSize.cols % 2 != 0 && kSize.rows % 2 != 0);
		int newheight = src.rows + kSize.rows - 1;
		int newwidth = src.cols + kSize.cols - 1;
		if (dst.empty()) {
			dst.create(newheight, newwidth, src.channels,0);
		}

		//! copy raw data
		_Tp* ptrSrc = NULL, *ptrDst = NULL;
		for (size_t i = kSize.rows/2; i <= dst.rows - kSize.rows/2; i++) {
			ptrSrc = src.data + (i - kSize.rows / 2)*src.step();
			ptrDst = dst.data + i*dst.step();
			int step = dst.step();
			size_t offsLeft = 0;
			size_t offsRight = (kSize.cols / 2 * src.channels + src.step());
			size_t cpySize = kSize.cols / 2 * src.channels;
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
				int cpys = cpySize / src.channels;
				//left side
				while (cpys--) {
					ptrS = ptrSrc + cpys*src.channels;
					for (size_t k = 0; k < src.channels; k++) {
						*(ptrD++) = *(ptrS++);
					}
				}
				//mid side
				memcpy(ptrDst + cpySize, ptrSrc, step * sizeof(_Tp));
				ptrD = ptrDst + offsRight;
				cpys = 0;
				//right side
				while (cpys++ < cpySize / src.channels) {
					ptrS = ptrSrc +(src.cols - cpys)*src.channels;
					for (size_t k = 0; k < src.channels; k++) {
						*(ptrD++) = *(ptrS++);
					}
				}
			}
		}

		size_t cpySize = kSize.rows / 2 * dst.step() * sizeof(_Tp);
		//! copy top-bottom boarder
		if (type == SAME) {
			_Tp* ptrDT = dst.data;
			_Tp* ptrST = dst.data + cpySize;
			_Tp* ptrDB = dst.data + dst.rows * dst.step() - cpySize;
			_Tp* ptrSB = dst.data + (dst.rows - kSize.rows) * dst.step();
			memcpy(ptrDT, ptrST, cpySize);
			memcpy(ptrDB, ptrSB, cpySize);
		} else if (type == MIRROR) {
			for (size_t i = 0; i < kSize.rows / 2; i++) {
				cpySize = dst.step() * sizeof(_Tp);
				_Tp* ptrDT = dst.data + i * dst.step();
				_Tp* ptrST = dst.data + ((kSize.rows -1) - i - 1) * dst.step();
				_Tp* ptrDB = dst.data + (dst.rows - 1 - i) * dst.step();
				_Tp* ptrSB = dst.data + (dst.rows - (kSize.rows - 1) + i) * dst.step();
				memcpy(ptrDT, ptrST, cpySize);
				memcpy(ptrDB, ptrSB, cpySize);
			}
		}
	}

	template<typename _Tp>
	void templateOp(const Matrix_<_Tp>& src,const Mat32f& kernel,Matrix_<_Tp>& dst,eBoarderType type) {
		assert(src.channels == 1 && kernel.channels == 1);
		Matrix_<_Tp> src_bd;
		copymakeBoarder(src, kernel.size(), type, src_bd);
		if (!dst.empty()) dst.release();
		dst.create(src.rows, src.cols, src.channels, 0);

		for (size_t i = kernel.cols / 2; i < src_bd.rows - kernel.cols / 2; i++) {
			for (size_t j = kernel.cols / 2; j < src_bd.cols - kernel.cols / 2; j++) {
				double sum = 0;
				double sumk = 0;
				for (int m = -kernel.rows/2; m < kernel.rows / 2; m++) {
					for (int n = -kernel.cols / 2; n < kernel.cols / 2; n++) {
						double v0 = (double)MAT_ELEM_S(src_bd, (i + m), (j + n));
						double v1 = (double)MAT_ELEM_S(kernel, (m + kernel.rows / 2), (n + kernel.cols / 2));
						sum += v0 * v1;							
						sumk += v1;
					}
				}
				sum /= (sumk+0.5);
				int i0 = i - kernel.rows / 2;
				int j0 = j - kernel.cols / 2;
				MAT_ELEM_S(dst,i0,j0) = static_cast<_Tp>(sum);
			}
		}
	}

	template<typename _Tp>
	void split(const Matrix_<_Tp>& src, std::vector<Matrix_<_Tp> >& mvs) {
		assert(src.checkValid());
		if (src.channels == 1) {
			mvs.push_back(src);
			return;
		}

		for (size_t c = 0; c < src.channels; c++) {
			Matrix_<_Tp> c_mat = Matrix_<_Tp>(src.rows,
				src.cols, 1);
			for (size_t i = 0; i < src.rows; i++) {
				for (size_t j = 0; j < src.cols; j++) {
					MAT_ELEM_S(c_mat,i,j) = MAT_ELEM_M(src, i, j, c);
				}
			}
			mvs.push_back(c_mat);
		}
	}

	template<typename _Tp>
	void merge(const std::vector<Matrix_<_Tp> >& mvs, Matrix_<_Tp>& dst) {
		if (mvs.empty()) return;
		int r = mvs[0].rows;
		int c = mvs[0].cols;
		int cn = (int)mvs.size();
		if (!dst.empty()) dst.release();
		dst.create(r, c, cn, 0);
#		ifdef USE_OMP
#		pragma omp parallel for
#		endif
		for (size_t k = 0; k < cn; k++) {
			for (size_t i = 0; i < r; i++) {
				for (size_t j = 0; j < c; j++) {
					MAT_ELEM_M(dst, i, j, k) = MAT_ELEM_S(mvs[k], i, j);
				}
			}
		}
	}

	template<typename _Tp>
	Matrix_<float> getGaussianKernel2D(const Size4i& kSize, float sita) {
		assert(kSize.height % 2 != 0 &&
			kSize.width % 2 != 0);
		Matrix_<float> kernel(kSize.height, kSize.width, 1);
		int i0 = kSize.height / 2;
		int j0 = kSize.width / 2;
		float A = 1. / (sigma * sqrt(2 * PI));
		for (size_t i = 0; i < kSize.height; i++) {
			for (size_t j = 0; j < kSize.width; j++) {
				kernel.at(i, j) = A * exp(-((i - i0)*(i - i0) / (2 * sigma*sigma) +
					(j - j0)*(j - j0) / (2 * sigma*sigma)));
			}
		}
		return kernel;
	}

	template<typename _Tp1, typename _Tp2>
	void cvt_img_fmt(const Matrix_<_Tp1>& src, Matrix_<_Tp2>& dst) {
		int i = 0;
		if (!dst.checkValid())
			dst.create(src.rows, src.cols, src.channels);
		for (; i < src.totalSizes() - 8; i += 8) {
			dst.data[i + 0] = (_Tp2)(src.data[i + 0]);
			dst.data[i + 1] = (_Tp2)(src.data[i + 1]);
			dst.data[i + 2] = (_Tp2)(src.data[i + 2]);
			dst.data[i + 3] = (_Tp2)(src.data[i + 3]);
			dst.data[i + 4] = (_Tp2)(src.data[i + 4]);
			dst.data[i + 5] = (_Tp2)(src.data[i + 5]);
			dst.data[i + 6] = (_Tp2)(src.data[i + 6]);
			dst.data[i + 7] = (_Tp2)(src.data[i + 7]);
		}

		for (; i < src.totalSizes(); i++) {
			dst.data[i] = (_Tp2)(src.data[i]);
		}
	}
};