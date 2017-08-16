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

		_Scalar<_Tp>(const _Scalar<_Tp>& rhs) {
			memcpy(_data, rhs._data, 4 * sizeof(_Tp));
		}

		static _Scalar<_Tp> all(_Tp val) {
			_Scalar<_Tp> s;
			//memset(s._data, val, 4 * sizeof(_Tp));
			s._data[0] = val;
			s._data[1] = val;
			s._data[2] = val;
			s._data[3] = val;
			return s;
		}
		friend std::ostream& operator <<(std::ostream& os, const _Scalar& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "[";
			os << rhs._data[0] << ",";
			os << rhs._data[1] << ",";
			os << rhs._data[2] << ",";
			os << rhs._data[3];
			os << "]";
			return os;
		};
		~_Scalar<_Tp>() {};

		_Tp& operator[](int idx) {
			return _data[idx];
		}

		_Scalar& operator =(const _Scalar& rhs) {
			if (this == &rhs) {
				return *this;
			}
			memcpy(_data, rhs._data, 4 * sizeof(_Tp));
			return *this;
		};

	public:
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
		friend std::ostream& operator <<(std::ostream& os, const _Size& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "(";
			os << rhs.width << "," << rhs.height;
			os << ")";
			return os;
		};
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

		friend std::ostream& operator <<(std::ostream& os, const Range_& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "[";
			os << rhs.start << "," << rhs.end;
			os << ")";
			return os;
		};
	public:
		_Tp start;
		_Tp end;
	};

	typedef Range_<int> Range4i;
	typedef Range_<float> Range4f;
	typedef Range_<double> Range4d;

	template<typename _Tp>
	class Rect_ {
	public:
		Rect_() {};
		Rect_(_Tp x, _Tp y, _Tp w, _Tp h) :
			x(x), y(y), width(w), height(h) {};
		friend std::ostream& operator <<(std::ostream& os, const Rect_& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "(";
			os << rhs.x << "," << rhs.y << "," << rhs.width << "," << rhs.height << end;
			os << ")";
			return os;
		};
	public:
		_Tp x;
		_Tp y;
		_Tp width;
		_Tp height;
	};

	typedef Rect_<int> Rect4i;
	typedef Rect_<float> Rect4f;
	typedef Rect_<double> Rect4d;

	template<typename _Tp>
	class Point_ {
	public:
		Point_() {};
		Point_(_Tp x, _Tp y) :x(x), y(y) {};
		friend std::ostream& operator <<(std::ostream& os, const Point_& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "(";
			os << rhs.x << "," << rhs.y << end;
			os << ")";
			return os;
		}
	public:
		_Tp x;
		_Tp y;
	};

	typedef Point_<int> Point4i;
	typedef Point_<float> Point4f;
	typedef Point_<double> Point4d;

	template<typename _Tp>
	class Matrix_
	{
	public:
		Matrix_();
		Matrix_(int h, int w, int c, void* data);
		Matrix_(int h, int w, int c);
		Matrix_(const Matrix_& rhs);
		~Matrix_();

		/**@brief create a matrix*/
		void create(int h, int w, int c, _Tp val = 205);

		/**@brief release matrix*/
		void release();

		/**@brief total size of matrix , h*w*c */
		size_t totalSizes() const;

		/**@brief total byte of matrix , h*w*c*sizeof(_Tp) */
		size_t totalBytes() const;

		/**@brief step of matrix , w*c */
		int step() const;

		/**@brief check matrix */
		bool checkValid() const;

		/**@brief check matrix empty or not */
		bool empty() const;

		/**@brief size of matrix , (w,h) */
		Size4i size() const;
		
		/**@brief access 1D element at data + offset(i0)*/
		_Tp& at(int i0) const;
		/**@brief access 2D element at data + i0*step + i1*channels */
		_Tp& at(int i0, int i1) const;
		/**@brief access 2D element with channel, address: data + i0*step + i1*channels + i2*/
		_Tp& at(int i0, int i1, int i2) const;

		/**@brief return i'th row address*/
		_Tp* ptr(int i0) const;

		/**@brief convert matrix to string */
		std::string toString(bool singleLine = true) const;

		/**@brief return a sub matrix */
		Matrix_ subMat(const Range4i& rowRange, const Range4i& colRange);

		/**@brief return a sub matrix */
		Matrix_ subMat(const Rect4i& roi);

		/**@brief matrix operator*/
		_Scalar<_Tp> dot(const Matrix_& m);

		/**@brief matrix operator*/
		Matrix_ hadamardProduct(const Matrix_& m);

#ifdef HAVE_OPENCV
		/**@brief interface of opencv Mat*/
		Matrix_(const cv::Mat& src, bool bCopy = true);

		/**@brief convert Matrix to opencv Mat*/
		cv::Mat to_cvmat();

		/**@brief convert cv::Mat to Matrix*/
		void from_cvmat(const cv::Mat& src);
#endif
		/**@brief stream output a matrix*/
		friend std::ostream& operator <<(std::ostream& os, const Matrix_& rhs) {
			os.setf(ios::fixed);
			os.precision(8);
			os << "[";
			for (size_t i = 0; i < rhs.rows; i++) {
				for (size_t j = 0; j < rhs.cols; j++) {
					for (size_t k = 0; k < rhs.channels; k++) {
						if (j < rhs.cols - 1) {
							if (typeid(_Tp) == typeid(uchar))
								os << (int)MAT_ELEM_M(rhs, i, j, k) << ",";
							else
								os << MAT_ELEM_M(rhs, i, j, k) << ",";
						} else {
							if (typeid(_Tp) == typeid(uchar))
								os << (int)MAT_ELEM_M(rhs, i, j, k) << ";";
							else
								os << MAT_ELEM_M(rhs, i, j, k) << ";";
						}
					}
				}
				if (i < rhs.rows - 1)
					os << endl;
				else
					os << "]";
			}
			return os;
		};
		/**@brief matrix operator*/
		Matrix_& operator =(const Matrix_& rhs);

		/**@brief matrix operator*/
		bool operator ==(const Matrix_& rhs);

		/**@brief matrix operator*/
		bool operator !=(const Matrix_& rhs);

		/**@brief matrix operator*/
		Matrix_ operator +(const Matrix_& rhs);

		/**@brief matrix operator*/
		Matrix_ operator +(const Scalar4d& scalar);

		/**@brief matrix operator*/
		Matrix_ operator -(const Matrix_& rhs);

		/**@brief matrix operator*/
		Matrix_ operator -(const Scalar4d& scalar);

		/**@brief matrix operator*/
		Matrix_ operator /(const Scalar4d& scalar);

		/**@brief matrix operator*/
		Matrix_ operator *(const Scalar4d& scalar);

		/**@brief matrix operator*/
		Matrix_ operator *(const Matrix_& rhs);

		/**@brief matrix transpose*/
		Matrix_ t();

		/**@brief return a row of matrix */
		Matrix_ row(int row);

		/**@brief return a col of matrix */
		Matrix_ col(int col);

		/**@brief return a row range of matrix */
		Matrix_ rowRange(const Range4i& rowRange);

		/**@brief return a col range of matrix */
		Matrix_ colRange(const Range4i& colRange);

		/**@brief copy self to other matrix */
		void copyTo(Matrix_& m, Rect4i roi) const;

		template<typename _Tp2>
		/**@brief convert to other type*/
		void convertTo(Matrix_<_Tp2>& dst) const;

		Matrix_ clone() const;

		/**@brief matrix concatenate.
		function will concat self and m on axes
		@param m -input matrix
		@param axes - concat axes. axes only support 2D,
			that's to say,axes == 0 => concat in row
			axes == 1 => concat in col
		*/
		Matrix_ concatenate(const Matrix_& m , int axes = 0);

		/**@brief return a matrix which filled by 0*/
		static Matrix_ zeros(int rows, int cols, int channels = 1);

		/**@brief return a matrix which filled by 1*/
		static Matrix_ ones(int rows, int cols, int channels = 1);

		/**@brief return a eye matrix*/
		static Matrix_ eyes(int rows, int cols);
		/**@brief shuffle matrix data with row
			    | 1 1 1 |      | 2 2 2 |    | 3 3 3 |
		shuffle(| 2 2 2 |)  => | 1 1 1 | or | 2 2 2 |
			    | 3 3 3 |      | 3 3 3 |    | 1 1 1 |
		*/
		static Matrix_ shuffle(const Matrix_& m);

		/**@brief return a random matrix*/
		static Matrix_ rand(int rows, int cols, int channels);

		/**@brief return a random matrix*/
		static Matrix_ randn(int rows, int cols, int channels, float esp = 1.f);
	public:
		_Tp* data;
		int cols;
		int rows;
		int channels;
	};

	template<typename _Tp1, typename _Tp2>
	void cvt_img_fmt(const Matrix_<_Tp1>& src, Matrix_<_Tp2>& dst);
}