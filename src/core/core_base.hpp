#pragma once
#include "../utils/utility.h"
#include "../pch.h"
#include <algorithm>
#include <vector>
#include <assert.h>
#include <iomanip>
namespace wcv {

#define IMAGE_8UC1	0x00000801
#define IMAGE_8UC2	0x00000802
#define IMAGE_8UC3	0x00000803
#define IMAGE_8UC4	0x00000804

#define IMAGE_16SC1	0x00000f01
#define IMAGE_16SC2	0x00000f02
#define IMAGE_16SC3	0x00000f03

#define IMAGE_16UC1	0x00000f01
#define IMAGE_16UC2	0x00000f02
#define IMAGE_16UC3	0x00000f03

#define IMAGE_32SC1	0x00002001
#define IMAGE_32SC2	0x00002002
#define IMAGE_32SC3	0x00002003

#define IMAGE_32FC1	0x00002001
#define IMAGE_32FC2	0x00002002
#define IMAGE_32FC3	0x00002003

#define IMAGE_64FC1	0x00004001
#define IMAGE_64FC2	0x00004002
#define IMAGE_64FC3	0x00004003

#define IMAGE_MAX_CN	8

#define IMAGE_TYPE_MASK  0x000000ff
#define IMAGE_ELEM_SIZE(type) (type >> 8 & IMAGE_TYPE_MASK) 	
#define IMAGE_ELEM_CN(type) (type & IMAGE_TYPE_MASK)
	
//	/**@brief get pixel element from single image*/
//#define MAT_ELEM_S(mat,r,c)			\
//		(*(mat.data + r*mat.step() + c))
//
//#define MAT_ELEM_M(mat,r,c,k)		\
//		(*(mat.data + r*mat.step() + c*mat.channels + k))
//
//	/**@brief	
//		"*" will be filled by zeros,if use NONE.	-> 000|123
//		"*" will be same with raw data, if use SAME -> 123|123
//		"*" will be filled mirror data,if use MIRROR -> 321|123	
//	*/
//	enum eBoarderType {
//		NONE = 0,
//		SAME,
//		MIRROR
//	};
//
//	template<typename _Tp>
//	class _Scalar
//	{
//	public:
//		_Scalar<_Tp>(_Tp v0, _Tp v1, _Tp v2, _Tp v3) {
//			_data[0] = v0, _data[1] = v1;
//			_data[2] = v2, _data[3] = v3;
//		};
//		_Scalar<_Tp>(_Tp v0, _Tp v1, _Tp v2) {
//			_data[0] = v0, _data[1] = v1;
//			_data[2] = v2, _data[3] = 0;
//		};
//		_Scalar<_Tp>(_Tp v0, _Tp v1) {
//			_data[0] = v0, _data[1] = v1;
//			_data[2] = 0, _data[3] = 0;
//		};
//		_Scalar<_Tp>(_Tp v0) {
//			_data[0] = v0, _data[1] = 0;
//			_data[2] = 0, _data[3] = 0;
//		};
//		_Scalar<_Tp>() {
//			memset(_data, 0, 4 * sizeof(_Tp));
//		}
//
//		static _Scalar<_Tp> all(_Tp val) {
//			_Scalar<_Tp> s;
//			memset(s._data, val, 4 * sizeof(_Tp));
//			return s;
//		}
//
//		~_Scalar<_Tp>() {};
//
//		_Tp& operator[](int idx) {
//			return _data[idx];
//		}
//
//	private:
//		_Tp _data[4];
//	};
//
//	typedef _Scalar<int> Scalar4i;
//	typedef _Scalar<float> Scalar4f;
//	typedef _Scalar<double> Scalar4d;
//
//	template<typename _Tp>
//	class _Size
//	{
//	public:
//		_Size() :rows(0), cols(0) {};
//		_Size(int w, int h) :cols(w), rows(h) {};
//		~_Size() {};
//
//		_Tp area() const { return rows*cols; }
//		_Tp rows;
//		_Tp cols;
//	};
//
//	typedef _Size<int> Size4i;
//	typedef _Size<float> Size4f;
//	typedef _Size<double> Size4d;
//
//	
//	template <typename _Tp = int>
//	/**@brief Range from "start" to "end",not include end
//	that's to say: [start,end)
//	*/
//	class Range_
//	{
//	public:
//		Range_() :start(0), end(0) {};
//		Range_(int start, int end) :
//			start(start), end(end) {}
//		~Range_() {};
//	public:
//		_Tp start;
//		_Tp end;
//	};
//
//	typedef Range_<int> Range4i;
//	typedef Range_<float> Range4f;
//	typedef Range_<double> Range4d;

	
//	template<typename _Tp>
//	class Matrix_
//	{
//	public:
//		_Tp* data;
//		int cols;
//		int rows;
//		int channels;
//		Matrix_() :data(nullptr),
//			cols(0), rows(0), channels(0) {};
//		Matrix_(int h, int w, int c, void* data) :
//			rows(h), cols(w), channels(c) , data((_Tp*)data) {};
//		Matrix_(int h, int w, int c) {
//			create(h, w, c);
//		};
//		Matrix_(const Matrix_& rhs) :
//			rows(rhs.rows), cols(rhs.cols),
//			channels(rhs.channels),data(nullptr) {
//			create(rows, cols, channels);
//			memcpy(this->data, rhs.data, rhs.totalBytes());
//		};
//
//		~Matrix_() {
//			release();
//		};
//
//		void create(int h, int w, int c, int val = 205) {
//			rows = h;
//			cols = w;
//			channels = c;
//			this->data = (_Tp*)fastAlloc(w*h*c * sizeof(_Tp));
//			memset(this->data, val, w*h*c * sizeof(_Tp));
//		}
//
//		void release() {
//			rows = 0;
//			cols = 0;
//			channels = 0;
//			if (this->data) {
//				fastFree(this->data);
//				this->data = nullptr;
//			}
//		}
//
//		size_t totalSizes() const {
//			return size_t(rows * cols * channels);
//		};
//
//		size_t totalBytes() const {
//			return size_t(rows * cols * channels * sizeof(_Tp));
//		}
//
//		int step() const {
//			return cols * channels;
//		};
//
//		bool checkValid() const {
//			bool bCheck = true;
//			if ((rows < 0 || cols < 0) ||
//				(rows == 0 && cols == 0) ||
//				data == nullptr)
//				bCheck = false;
//			return bCheck;
//		};
//
//		bool empty() const {
//			return (rows*cols*channels == 0) || (data == nullptr);
//		};
//
//		Size4i size() const { return Size4i(cols, rows); };
//
//		std::string toString(bool singleLine = true) const {
//			std::stringstream buff;
//			buff << "[";
//			for (size_t i = 0; i < rows; i++) {
//				for (size_t j = 0; j < cols; j++) {
//					for (size_t k = 0; k < channels; k++) {
//						if (j < cols - 1)
//							buff << (float)MAT_ELEM_M((*this), i, j, k) << ",";
//						else
//							buff << (float)MAT_ELEM_M((*this), i, j, k) << ";";
//					}
//				}
//				if (!singleLine && i < rows - 1)
//					buff << endl;
//				else if (i == rows - 1)
//					buff << "]";
//			}
//			return buff.str();
//		};
//
//		friend std::ostream& operator <<(std::ostream& os, const Matrix_& rhs) {
//			os.setf(ios::fixed);
//			os.precision(8);
//			os << "[";
//			for (size_t i = 0; i < rhs.rows; i++) {
//				for (size_t j = 0; j < rhs.cols; j++) {
//					for (size_t k = 0; k < rhs.channels; k++) {
//						if (j < rhs.cols - 1)
//							os << MAT_ELEM_M(rhs, i, j, k) << ",";
//						else
//							os << MAT_ELEM_M(rhs, i, j, k) << ";";
//					}
//				}
//				if (i < rhs.rows - 1)
//					os << endl;
//				else
//					os << "]";
//			}
//			return os;
//		}
//
//		Matrix_ subMat(const Range4i& rowRange, const Range4i& colRange) {
//			assert(rowRange.end >= rowRange.start &&
//				rowRange.end <= rows && rowRange.start >= 0);
//			assert(colRange.end >= colRange.start &&
//				colRange.end <= cols && colRange.start >= 0);
//			Matrix_ submat(rowRange.end - rowRange.start,
//				colRange.end - colRange.start,channels);
//			size_t cpySize = (colRange.end - colRange.start)*channels;
//			for (size_t i = rowRange.start; i < rowRange.end; i++) {
//				uchar* ptrS = this->data + i*step() + colRange.start*channels;
//				uchar* ptrD = submat.data + i*submat.step();
//				memcpy(ptrD, ptrS, cpySize * sizeof(_Tp));
//			}
//			return submat;
//		}
//
//		_Scalar<_Tp> dot(const Matrix_& m) {
//			assert(rows == m.rows &&
//				cols == m.cols && 
//				channels == m.channels);
//			int i = 0;
//			_Scalar<_Tp> res = _Scalar<_Tp>::all(0);
//			if (channels == 1) {
//#			ifdef USE_OMP
//#			pragma omp parallel for reduction(+:res)
//#			endif
//				for (; i < totalSizes() - 8; i += 8) {
//					res[0] += *(data + i + 0) * *(m.data + i + 0);
//					res[0] += *(data + i + 1) * *(m.data + i + 1);
//					res[0] += *(data + i + 2) * *(m.data + i + 2);
//					res[0] += *(data + i + 3) * *(m.data + i + 3);
//					res[0] += *(data + i + 4) * *(m.data + i + 4);
//					res[0] += *(data + i + 5) * *(m.data + i + 5);
//					res[0] += *(data + i + 6) * *(m.data + i + 6);
//					res[0] += *(data + i + 7) * *(m.data + i + 7);
//				}
//
//#			ifdef USE_OMP
//#			pragma omp parallel for reduction(+:res)
//#			endif
//				for (; i < totalSizes(); i++) {
//					res[0] += *(data + i) * *(m.data + i);
//				}
//			} else {
//				std::vector<Matrix_<_Tp> > mvs;
//				split(m, mvs);
//			}
//			return res;
//		}
//
//#ifdef HAVE_OPENCV
//		Matrix_(const cv::Mat& src, bool bCopy = true) {
//			CV_Assert(!src.empty());
//			int h = src.rows;
//			int w = src.cols;
//			int c = src.channels();
//			this->rows = h;
//			this->cols = w;
//			this->channels = c;
//			if (bCopy) {
//				if (this->data)
//					delete[] this->data;
//				this->data = (uchar*)fastAlloc(w*h*c * sizeof(uchar));
//				memset(this->data, 0, w*h*c * sizeof(uchar));
//			}
//			else {
//				this->data = (uchar*)src.ptr<uchar>();
//			}
//		};
//
//		cv::Mat to_cvmat() {
//			cv::Mat dst;
//			if (channels == 1) {
//				dst = cv::Mat(rows, cols, CV_8UC1);
//			}
//			else if (channels == 3) {
//				dst = cv::Mat(rows, cols, CV_8UC3);
//			}
//			else {
//				fprintf(stderr, "unsupport format image.");
//				return cv::Mat();
//			}
//			memcpy(dst.ptr<uchar>(), this->data, this->totalSizes() * sizeof(uchar));
//			return dst;
//		};
//
//		void from_cvmat(const cv::Mat& src) {
//			CV_Assert(!src.empty());
//			this->rows = src.rows;
//			this->cols = src.cols;
//			this->channels = src.channels();
//			if (data)
//				delete[] data;
//			size_t size = rows*cols*channels;
//			data = (uchar*)fastAlloc(size * sizeof(uchar));
//			memcpy(data, src.ptr<uchar>(), size * sizeof(uchar));
//		};
//#endif
//
//		Matrix_& operator =(const Matrix_& rhs) {
//			if (this == &rhs) {
//				return *this;
//			}
//			release();
//			create(rhs.rows, rhs.cols, rhs.channels);
//			size_t total = this->totalBytes();
//			memcpy(this->data, rhs.data, total);
//			return *this;
//		}
//	};


#if 0

	class Image
	{
	public:
		uchar* data;
		int cols;
		int rows;
		//int channels;
		int type;
		Image() :data(nullptr),
			cols(0), rows(0), /*channels(0),*/ type(0) {};
		Image(int h, int w, int type, unsigned char* data) :
			rows(h), cols(w), type(type)/* channels(c),*/, data(data) {};
		Image(int h, int w, int type) {
			create(h, w, type);
		};
		~Image() {
			release();
		};

		template<typename _Tp>
		_Tp* ptr() const {
			return (_Tp*)this->data;
		}

		template<typename _Tp>
		_Tp* ptr(int i0) const {
			return (_Tp*)(this->data + i0 * this->step());
		}

		template<typename _Tp>
		_Tp* ptr(int i0, int i1) const {
			return (_Tp*)(this->data + i0 * this->step() + i1*this->channels);
		}

		void create(int h, int w, int type, int val = 205) {
			rows = h, cols = w;
			int c = IMAGE_ELEM_CN(type);
			int bit = IMAGE_ELEM_SIZE(type);
			//! chack max channels and bit is power of 2
			assert(c < IMAGE_MAX_CN && (bit & (bit - 1) == 0));
			int elem_size = bit / sizeof(uchar);
			int total = h * w * c * elem_size * sizeof(uchar);
			this->data = (uchar*)fastAlloc(total);
			memset(this->data, val, total);
		}

		void release() {
			rows = 0;
			cols = 0;
			//channels = 0;
			if (this->data) {
				fastFree(this->data);
				this->data = nullptr;
			}
		}

		int channels() const {
			return IMAGE_ELEM_CN(this->type);
		}

		size_t totalSizes() const {
			return size_t(rows * cols * channels() * IMAGE_ELEM_SIZE(type) / sizeof(uchar));
		};

		size_t totalSizes() const {
			return size_t(rows * cols * channels());
		}

		int step() const {
			return cols * channels();
		};

		bool checkValid() const {
			bool bCheck = true;
			if ((rows < 0 || cols < 0) ||
				(rows == 0 && cols == 0) ||
				data == nullptr)
				bCheck = false;
			return bCheck;
		};

		bool empty() const {
			return (rows*cols*channels() == 0) || (data == nullptr);
		};

		Size4i size() const { return Size4i(cols, rows); };

#ifdef HAVE_OPENCV
		Image(const cv::Mat& src, bool bCopy = true) {
			CV_Assert(!src.empty());

			int h = src.rows;
			int w = src.cols;
			//int c = src.channels();
			this->rows = h;
			this->cols = w;

			//this->channels = c;
			cv::Mat s = src;
			size_t total = CV_ELEM_SIZE(src.flags);
			for (int i = s.dims - 1; i >= 0; i--) {
				if (step) {
					if (0 && s.step[i] != CV_AUTOSTEP) {
						CV_Assert(total <= src.step[i]);
						total = s.step[i];
					}else
						s.step[i] = total;
				}
				total *= s.size[i];
			}

			if (bCopy) {
				if (this->data)
					delete[] this->data;
				this->data = (uchar*)fastAlloc(total);
				memset(this->data, 0, total);
			}
			else {
				this->data = (uchar*)src.ptr<uchar>();
			}
		};

		cv::Mat to_cvmat() {
			cv::Mat dst;
			if (channels() == 1) {
				dst = cv::Mat(rows, cols, CV_8UC1);
			}
			else if (channels() == 3) {
				dst = cv::Mat(rows, cols, CV_8UC3);
			}
			else {
				fprintf(stderr, "unsupport format image.");
				return cv::Mat();
			}
			memcpy(dst.ptr<uchar>(), this->data, this->totalSizes() * sizeof(uchar));
			return dst;
		};

		void from_cvmat(const cv::Mat& src) {
			CV_Assert(!src.empty());
			this->rows = src.rows;
			this->cols = src.cols;
			//this->channels = src.channels();
			if (data)
				delete[] data;
			size_t size = rows*cols*channels();
			cv::Mat s = src;
			size_t total = CV_ELEM_SIZE(src.flags);
			for (int i = s.dims - 1; i >= 0; i--) {
				if (step) {
					if (0 && s.step[i] != CV_AUTOSTEP) {
						CV_Assert(total <= src.step[i]);
						total = s.step[i];
					}
					else
						s.step[i] = total;
				}
				total *= s.size[i];
			}

			data = (uchar*)fastAlloc(total/*size * sizeof(uchar)*/);
			memcpy(data, src.ptr<uchar>(), total/*size * sizeof(uchar)*/);
		};
#endif
	};
#endif // 0


	//template<typename _Tp1, typename _Tp2>
	//void cvt_img_fmt(const Matrix_<_Tp1>& src, Matrix_<_Tp2>& dst) {
	//	int i = 0;
	//	if (!dst.checkValid())
	//		dst.create(src.rows, src.cols, src.channels);
	//	for (; i < src.totalSizes() - 8; i += 8) {
	//		dst.data[i + 0] = (_Tp2)(src.data[i + 0]);
	//		dst.data[i + 1] = (_Tp2)(src.data[i + 1]);
	//		dst.data[i + 2] = (_Tp2)(src.data[i + 2]);
	//		dst.data[i + 3] = (_Tp2)(src.data[i + 3]);
	//		dst.data[i + 4] = (_Tp2)(src.data[i + 4]);
	//		dst.data[i + 5] = (_Tp2)(src.data[i + 5]);
	//		dst.data[i + 6] = (_Tp2)(src.data[i + 6]);
	//		dst.data[i + 7] = (_Tp2)(src.data[i + 7]);
	//	}

	//	for (; i < src.totalSizes(); i++) {
	//		dst.data[i] = (_Tp2)(src.data[i]);
	//	}
	//}
}