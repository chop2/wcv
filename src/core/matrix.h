#pragma once
#include "../pch.h"

namespace wcv {
	/**@brief get pixel element from single image*/
#define MAT_ELEM_S(mat,r,c)			\
		(*(mat.data + r*mat.step() + c))

#define MAT_ELEM_M(mat,r,c,k)		\
		(*(mat.data + r*mat.step() + c*mat.channels + k))

	/**@brief
	"*" will be filled by zeros,if use NONE.	-> 000|123
	"*" will be same with raw data, if use SAME -> 123|123
	"*" will be filled mirror data,if use MIRROR -> 321|123
	*/
	enum eBoarderType {
		NONE = 0,
		SAME,
		MIRROR
	};

	template<typename _Tp>
	class _Scalar
	{
	public:
		_Scalar<_Tp>(_Tp v0, _Tp v1, _Tp v2, _Tp v3) {
			_data[0] = v0, _data[1] = v1;
			_data[2] = v2, _data[3] = v3;
		};
		_Scalar<_Tp>(_Tp v0, _Tp v1, _Tp v2) {
			_data[0] = v0, _data[1] = v1;
			_data[2] = v2, _data[3] = 0;
		};
		_Scalar<_Tp>(_Tp v0, _Tp v1) {
			_data[0] = v0, _data[1] = v1;
			_data[2] = 0, _data[3] = 0;
		};
		_Scalar<_Tp>(_Tp v0) {
			_data[0] = v0, _data[1] = 0;
			_data[2] = 0, _data[3] = 0;
		};
		_Scalar<_Tp>() {
			memset(_data, 0, 4 * sizeof(_Tp));
		}

		static _Scalar<_Tp> all(_Tp val) {
			_Scalar<_Tp> s;
			memset(s._data, val, 4 * sizeof(_Tp));
			return s;
		}

		~_Scalar<_Tp>() {};

		_Tp& operator[](int idx) {
			return _data[idx];
		}

	private:
		_Tp _data[4];
	};

	typedef _Scalar<int> Scalar4i;
	typedef _Scalar<float> Scalar4f;
	typedef _Scalar<double> Scalar4d;

	template<typename _Tp>
	class _Size
	{
	public:
		_Size() :height(0), width(0) {};
		_Size(int w, int h) :width(w), height(h) {};
		~_Size() {};

		_Tp area() const { return height*width; }
		_Tp height;
		_Tp width;
	};

	typedef _Size<int> Size4i;
	typedef _Size<float> Size4f;
	typedef _Size<double> Size4d;


	template <typename _Tp = int>
	/**@brief Range from "start" to "end",not include end
	that's to say: [start,end)
	*/
	class Range_
	{
	public:
		Range_() :start(0), end(0) {};
		Range_(int start, int end) :
			start(start), end(end) {}
		~Range_() {};
	public:
		_Tp start;
		_Tp end;
	};

	typedef Range_<int> Range4i;
	typedef Range_<float> Range4f;
	typedef Range_<double> Range4d;

	template<typename _Tp>
	class Matrix_
	{
	public:
		Matrix_();
		Matrix_(int h, int w, int c, void* data);
		Matrix_(int h, int w, int c);
		Matrix_(const Matrix_& rhs);
		~Matrix_();

		void create(int h, int w, int c, int val = 205);

		void release();

		size_t totalSizes() const;

		size_t totalBytes() const;

		int step() const;

		bool checkValid() const;

		bool empty() const;

		Size4i size() const;
		
		/**@brief access 1D element at data + offset(i0)*/
		_Tp& at(int i0) const;
		/**@brief access 2D element at data + i0*step + i1*channels */
		_Tp& at(int i0, int i1) const;
		/**@brief access 2D element with channel, address: data + i0*step + i1*channels + i2*/
		_Tp& at(int i0, int i1, int i2) const;

		/**@brief return i'th row address*/
		_Tp* ptr(int i0);

		std::string toString(bool singleLine = true) const;

		friend std::ostream& operator <<(std::ostream& os, const Matrix_& rhs);

		Matrix_ subMat(const Range4i& rowRange, const Range4i& colRange);

		_Scalar<_Tp> dot(const Matrix_& m);

		Matrix_ hadamardProduct(const Matrix_& m);

#ifdef HAVE_OPENCV
		Matrix_(const cv::Mat& src, bool bCopy = true);

		cv::Mat to_cvmat();

		void from_cvmat(const cv::Mat& src);
#endif
		Matrix_& operator =(const Matrix_& rhs);

	public:
		_Tp* data;
		int cols;
		int rows;
		int channels;
	};

	template<typename _Tp1, typename _Tp2>
	void cvt_img_fmt(const Matrix_<_Tp1>& src, Matrix_<_Tp2>& dst);
}