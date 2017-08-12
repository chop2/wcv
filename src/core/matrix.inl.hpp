#pragma once
#include "matrix.h"
#include <random>
#include <chrono>

namespace wcv {

	template<typename _Tp>
	double _randn(double min, double max) {
		//return min + (max - min)*std::rand() / (RAND_MAX + 1.0);
		std::random_device rd;
		std::mt19937 gen(rd());
		std::normal_distribution<> d(0, 1);
		return min + d(gen) * (max - min);
	}

	template<typename _Tp>
	double _rand(double min,double max) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(min,max);
		return dis(gen);
	}

	template<typename _Tp>
	Matrix_<_Tp>::Matrix_() :data(nullptr),
		cols(0), rows(0), channels(0) {};

	template<typename _Tp>
	Matrix_<_Tp>::Matrix_(int h, int w, int c, void* data) :
		rows(h), cols(w), channels(c), data((_Tp*)data) {};;

	template<typename _Tp>
	Matrix_<_Tp>::Matrix_(int h, int w, int c) {
		create(h, w, c);
	};

	template<typename _Tp>
	Matrix_<_Tp>::Matrix_(const Matrix_& rhs) :
		rows(rhs.rows), cols(rhs.cols),
		channels(rhs.channels), data(nullptr) {
		create(rows, cols, channels);
		memcpy(this->data, rhs.data, rhs.totalBytes());
	};

	template<typename _Tp>
	Matrix_<_Tp>::~Matrix_() {
		release();
	};

	template<typename _Tp>
	void Matrix_<_Tp>::create(int h, int w, int c, int val) {
		int i = 0;
		rows = h;
		cols = w;
		channels = c;
		this->data = (_Tp*)fastAlloc(w*h*c * sizeof(_Tp));
		// if type of _Tp is float,double,memset not work correctly
		//memset(this->data, val, w*h*c * sizeof(_Tp));
		if (totalSizes() > 8) {
			for (; i < totalSizes() - 8; i += 8) {
				*(this->data + i + 0) = val;
				*(this->data + i + 1) = val;
				*(this->data + i + 2) = val;
				*(this->data + i + 3) = val;
				*(this->data + i + 4) = val;
				*(this->data + i + 5) = val;
				*(this->data + i + 6) = val;
				*(this->data + i + 7) = val;
			}
		}
		for (; i < totalSizes(); i++) {
			*(this->data + i) = val;
		}
	};

	template<typename _Tp>
	void Matrix_<_Tp>::release() {
		rows = 0,cols = 0;
		channels = 0;
		if (this->data) {
			fastFree(this->data);
			this->data = nullptr;
		}
	};

	template<typename _Tp>
	size_t Matrix_<_Tp>::totalSizes() const {
		return size_t(rows * cols * channels);
	};

	template<typename _Tp>
	size_t Matrix_<_Tp>::totalBytes() const {
		return size_t(rows * cols * channels * sizeof(_Tp));
	};

	template<typename _Tp>
	int Matrix_<_Tp>::step() const {
		return cols * channels;
	};

	template<typename _Tp>
	bool Matrix_<_Tp>::checkValid() const {
		bool bCheck = true;
		if ((rows < 0 || cols < 0) ||
			(rows == 0 && cols == 0) ||
			data == nullptr)
			bCheck = false;
		return bCheck;
	};

	template<typename _Tp>
	bool Matrix_<_Tp>::empty() const {
		return (rows*cols*channels == 0) || (data == nullptr);
	};

	template<typename _Tp>
	Size4i Matrix_<_Tp>::size() const { return Size4i(cols, rows); }

	template<typename _Tp>
	inline _Tp & Matrix_<_Tp>::at(int i0) const	{
		return (_Tp)*(this->data + i0);
	}

	template<typename _Tp>
	inline _Tp & Matrix_<_Tp>::at(int i0, int i1) const	{
		return (_Tp)*(this->data + i0 * step() + i1*channels);
	}

	template<typename _Tp>
	inline _Tp & Matrix_<_Tp>::at(int i0, int i1, int i2) const	{
		return (_Tp)*(this->data + i0 * step() + i1*channels + i2);
	}

	template<typename _Tp>
	inline _Tp * Matrix_<_Tp>::ptr(int i0) const{
		return (_Tp*)(this->data + i0 * this->step());
	}

	template<typename _Tp>
	std::string Matrix_<_Tp>::toString(bool singleLine) const {
		std::stringstream buff;
		buff << "[";
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				for (size_t k = 0; k < channels; k++) {
					if (j < cols - 1)
						buff << (float)MAT_ELEM_M((*this), i, j, k) << ",";
					else
						buff << (float)MAT_ELEM_M((*this), i, j, k) << ";";
				}
			}
			if (!singleLine && i < rows - 1)
				buff << endl;
			else if (i == rows - 1)
				buff << "]";
		}
		return buff.str();
	};

	//template<typename _Tp>
	//std::ostream& operator<<(std::ostream& os, const Matrix_<_Tp>& rhs) {
	//	os.setf(ios::fixed);
	//	os.precision(8);
	//	os << "[";
	//	for (size_t i = 0; i < rhs.rows; i++) {
	//		for (size_t j = 0; j < rhs.cols; j++) {
	//			for (size_t k = 0; k < rhs.channels; k++) {
	//				if (j < rhs.cols - 1)
	//					os << MAT_ELEM_M(rhs, i, j, k) << ",";
	//				else
	//					os << MAT_ELEM_M(rhs, i, j, k) << ";";
	//			}
	//		}
	//		if (i < rhs.rows - 1)
	//			os << endl;
	//		else
	//			os << "]";
	//	}
	//	return os;
	//};

	template<typename _Tp>
	Matrix_<_Tp> Matrix_<_Tp>::subMat(const Range4i& rowRange,
		const Range4i& colRange) {
		assert(rowRange.end >= rowRange.start &&
			rowRange.end <= rows && rowRange.start >= 0);
		assert(colRange.end >= colRange.start &&
			colRange.end <= cols && colRange.start >= 0);
		Matrix_ submat;
		submat.create(rowRange.end - rowRange.start,
			colRange.end - colRange.start, channels,0);
		size_t cpySize = (colRange.end - colRange.start)*channels;
		for (size_t i = rowRange.start; i < rowRange.end; i++) {
			_Tp* ptrS = this->data + i*step() + colRange.start*channels;
			_Tp* ptrD = submat.data + (i-rowRange.start)*submat.step();
			memcpy(ptrD, ptrS, cpySize * sizeof(_Tp));
		}
		return submat;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::subMat(const Rect4i & roi) {
		Range4i row_range(roi.y, roi.y + roi.height);
		Range4i col_range(roi.x, roi.x + roi.width);
		return subMat(row_range, col_range);
	}

	template<typename _Tp>
	_Scalar<_Tp> Matrix_<_Tp>::dot(const Matrix_& m) {
		assert(rows == m.rows &&
			cols == m.cols &&
			channels == m.channels);
		_Scalar<_Tp> scRes = _Scalar<_Tp>::all(_Tp(0));

		auto process_single = [](const Matrix_<_Tp>& s0, 
			const Matrix_<_Tp>& s1)->_Tp {
			int i = 0;
			_Tp res = 0;
			if (s0.totalSizes() > 8) {
#			ifdef USE_OMP
#			pragma omp parallel for reduction(+:res)
#			endif
				for (; i < s0.totalSizes() - 8; i += 8) {
					res += *(s0.data + i + 0) * *(s1.data + i + 0);
					res += *(s0.data + i + 1) * *(s1.data + i + 1);
					res += *(s0.data + i + 2) * *(s1.data + i + 2);
					res += *(s0.data + i + 3) * *(s1.data + i + 3);
					res += *(s0.data + i + 4) * *(s1.data + i + 4);
					res += *(s0.data + i + 5) * *(s1.data + i + 5);
					res += *(s0.data + i + 6) * *(s1.data + i + 6);
					res += *(s0.data + i + 7) * *(s1.data + i + 7);
				}

#			ifdef USE_OMP
#			pragma omp parallel for reduction(+:res)
#			endif
			}
			for (; i < s0.totalSizes(); i++) {
				res += *(s0.data + i) * *(s1.data + i);
			}
			return res;
		};

		if (channels == 1) {
			scRes[0] = process_single(*this, m);
		}
		else {
			std::vector<Matrix_<_Tp> > mvs0,mvs1;
			split(*this, mvs0);
			split(m, mvs1);
#			ifdef USE_OMP
#			pragma omp parallel for
#			endif
			for (int i = 0; i < mvs0.size(); i++) {
				scRes[i] = process_single(mvs0[i], mvs1[i]);
			}
		}
		return scRes;
	};

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::hadamardProduct(const Matrix_<_Tp> & m)	{
		assert(rows == m.rows);
		assert(cols == m.cols);
		Matrix_<_Tp> result;

		auto process_single = [](const Matrix_<_Tp>& s0,
			const Matrix_<_Tp>& s1)->Matrix_<_Tp> {
			int i = 0;
			Matrix_<_Tp> res(s0.rows,s0.cols,s0.channels);
			if (s0.totalSizes() > 8) {
#			ifdef USE_OMP
#			pragma omp parallel for reduction(+:res)
#			endif
				for (; i < s0.totalSizes() - 8; i += 8) {
					*(res.data + i + 0) = *(s0.data + i + 0) * *(s1.data + i + 0);
					*(res.data + i + 1) = *(s0.data + i + 1) * *(s1.data + i + 1);
					*(res.data + i + 2) = *(s0.data + i + 2) * *(s1.data + i + 2);
					*(res.data + i + 3) = *(s0.data + i + 3) * *(s1.data + i + 3);
					*(res.data + i + 4) = *(s0.data + i + 4) * *(s1.data + i + 4);
					*(res.data + i + 5) = *(s0.data + i + 5) * *(s1.data + i + 5);
					*(res.data + i + 6) = *(s0.data + i + 6) * *(s1.data + i + 6);
					*(res.data + i + 7) = *(s0.data + i + 7) * *(s1.data + i + 7);
				}

#			ifdef USE_OMP
#			pragma omp parallel for reduction(+:res)
#			endif
			}
			for (; i < s0.totalSizes(); i++) {
				*(res.data + i) = *(s0.data + i) * *(s1.data + i);
			}
			return res;
		};

		if (channels == 1) {
			result = process_single(*this, m);
		}
		else {
			std::vector<Matrix_<_Tp> > mvs0, mvs1,mvs;
			split(*this, mvs0);
			split(m, mvs1);
#			ifdef USE_OMP
#			pragma omp parallel for
#			endif
			for (int i = 0; i < mvs0.size(); i++) {
				Matrix_<_Tp> m = process_single(mvs0[i], mvs1[i]);
				mvs.push_back(m);
			}
			merge(mvs, result);
		}

		return result;
	};

#ifdef HAVE_OPENCV
	template<typename _Tp>
	Matrix_<_Tp>::Matrix_(const cv::Mat& src, bool bCopy) {
		CV_Assert(!src.empty());
		int h = src.rows;
		int w = src.cols;
		int c = src.channels();
		this->rows = h;
		this->cols = w;
		this->channels = c;
		if (bCopy) {
			if (this->data)
				delete[] this->data;
			this->data = (uchar*)fastAlloc(w*h*c * sizeof(uchar));
			memset(this->data, 0, w*h*c * sizeof(uchar));
		}
		else {
			this->data = (uchar*)src.ptr<uchar>();
		}
	};

	template<typename _Tp>
	cv::Mat Matrix_<_Tp>::to_cvmat() {
		cv::Mat dst;
		if (channels == 1) {
			dst = cv::Mat(rows, cols, CV_8UC1);
		}
		else if (channels == 3) {
			dst = cv::Mat(rows, cols, CV_8UC3);
		}
		else {
			fprintf(stderr, "unsupport format image.");
			return cv::Mat();
		}
		memcpy(dst.ptr<uchar>(), this->data, this->totalSizes() * sizeof(uchar));
		return dst;
	};

	template<typename _Tp>
	void Matrix_<_Tp>::from_cvmat(const cv::Mat& src) {
		CV_Assert(!src.empty());
		this->rows = src.rows;
		this->cols = src.cols;
		this->channels = src.channels();
		if (data)
			delete[] data;
		size_t size = rows*cols*channels;
		data = (uchar*)fastAlloc(size * sizeof(uchar));
		memcpy(data, src.ptr<uchar>(), size * sizeof(uchar));
	};
#endif

	template<typename _Tp>
	Matrix_<_Tp>& Matrix_<_Tp>::operator =(const Matrix_& rhs) {
		if (this == &rhs) {
			return *this;
		}
		release();
		create(rhs.rows, rhs.cols, rhs.channels);
		size_t total = this->totalBytes();
		memcpy(this->data, rhs.data, total);
		return *this;
	}

	template<typename _Tp>
	inline bool Matrix_<_Tp>::operator==(const Matrix_ & rhs) {
		assert(rhs.checkValid());
		bool bres = true;
		if (rows == rhs.rows &&
			cols == rhs.cols &&
			channels == rhs.channels) {
			int res = memcmp(this->data, rhs.data, this->totalBytes());
			bres = res == 0 ? true : false;
		} else {
			bres = false;
		}
		return bres;
	}

	template<typename _Tp>
	inline bool Matrix_<_Tp>::operator!=(const Matrix_ & rhs) {
		return !(*this == rhs);
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator+(const Matrix_ & rhs) {
		assert(checkValid() && rhs.checkValid());
		assert(rows == rhs.rows &&
			cols == rhs.cols &&
			channels == rhs.channels);
		Matrix_<_Tp> dst(rhs.rows, rhs.cols, rhs.channels);
		int i = 0;
		if (totalSizes() > 8) {
			for (; i < rhs.totalSizes() - 8; i += 8) {
				*(dst.data + i + 0) = *(data + i + 0) + *(rhs.data + i + 0);
				*(dst.data + i + 1) = *(data + i + 1) + *(rhs.data + i + 1);
				*(dst.data + i + 2) = *(data + i + 2) + *(rhs.data + i + 2);
				*(dst.data + i + 3) = *(data + i + 3) + *(rhs.data + i + 3);
				*(dst.data + i + 4) = *(data + i + 4) + *(rhs.data + i + 4);
				*(dst.data + i + 5) = *(data + i + 5) + *(rhs.data + i + 5);
				*(dst.data + i + 6) = *(data + i + 6) + *(rhs.data + i + 6);
				*(dst.data + i + 7) = *(data + i + 7) + *(rhs.data + i + 7);
			}
		}

		for (; i < rhs.totalSizes(); i++) {
			*(dst.data + i) = *(data + i) + *(rhs.data + i);
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator+(const Scalar4d & scalar) {
		assert(channels < 5);
		Matrix_<_Tp> dst;
		dst.create(rows, cols, channels, 0);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				_Tp* ptrS = data + i*step() + j*channels;
				_Tp* ptrD = dst.data + i*step() + j*channels;
				for (size_t k = 0; k < channels; k++) {
					*(ptrD++) = *(ptrS++) + (_Tp)(scalar._data[k]);
				}
			}
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator-(const Matrix_ & rhs) {
		assert(checkValid() && rhs.checkValid());
		assert(rows == rhs.rows &&
			cols == rhs.cols &&
			channels == rhs.channels);
		Matrix_<_Tp> dst;
		dst.create(rhs.rows, rhs.cols, rhs.channels);
		int i = 0;
		if (totalSizes() > 8) {
			for (; i < rhs.totalSizes() - 8; i += 8) {
				*(dst.data + i + 0) = *(data + i + 0) - *(rhs.data + i + 0);
				*(dst.data + i + 1) = *(data + i + 1) - *(rhs.data + i + 1);
				*(dst.data + i + 2) = *(data + i + 2) - *(rhs.data + i + 2);
				*(dst.data + i + 3) = *(data + i + 3) - *(rhs.data + i + 3);
				*(dst.data + i + 4) = *(data + i + 4) - *(rhs.data + i + 4);
				*(dst.data + i + 5) = *(data + i + 5) - *(rhs.data + i + 5);
				*(dst.data + i + 6) = *(data + i + 6) - *(rhs.data + i + 6);
				*(dst.data + i + 7) = *(data + i + 7) - *(rhs.data + i + 7);
			}
		}

		for (; i < rhs.totalSizes(); i++) {
			*(dst.data + i) = *(data + i) - *(rhs.data + i);
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator-(const Scalar4d & scalar) {
		assert(channels < 5);
		Matrix_<_Tp> dst;
		dst.create(rows, cols, channels, 0);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				_Tp* ptrS = data + i*step() + j*channels;
				_Tp* ptrD = dst.data + i*step() + j*channels;
				for (size_t k = 0; k < channels; k++) {
					*(ptrD++) = *(ptrS++) - (_Tp)(scalar._data[k]);
				}
			}
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator/(const Scalar4d & scalar) {
		assert(channels < 5);
		Matrix_<_Tp> dst;
		dst.create(rows, cols, channels, 0);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				_Tp* ptrS = data + i*step() + j*channels;
				_Tp* ptrD = dst.data + i*step() + j*channels;
				for (size_t k = 0; k < channels; k++) {
					*(ptrD++) = (_Tp)(*(ptrS++) / scalar._data[k]);
				}
			}
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator*(const Scalar4d& scalar) {
		assert(channels < 5);
		Matrix_<_Tp> dst;
		dst.create(rows, cols, channels, 0);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				_Tp* ptrS = data + i*step() + j*channels;
				_Tp* ptrD = dst.data + i*step() + j*channels;
				for (int k = 0; k < channels; k++) {
					double v = scalar._data[k];
					*(ptrD++) = (_Tp)(*(ptrS++) * v);
				}
			}
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::operator*(const Matrix_ & rhs) {
		assert(cols == rhs.rows);
		assert(channels == 1);
		int nr = rows;
		int nc = rhs.cols;
		Matrix_<_Tp> dst;
		dst.create(nr, nc, 1, 0);
		for (size_t i = 0; i < nr; i++)	{
			for (size_t j = 0; j < nc; j++) {
				for (size_t idx = 0; idx < cols; idx++)	{
					dst.at(i, j) += this->at(i, idx) * rhs.at(idx, j);
				}
			}
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::rowRange(const Range4i & row_range) {
		assert(row_range.end >= row_range.start &&
			row_range.start >= 0 && row_range.end <= rows);
		return subMat(row_range,Range4i(0,cols));
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::colRange(const Range4i & col_range) {
		assert(col_range.end >= col_range.start &&
			col_range.start >= 0 && col_range.end <= cols);
		return subMat(Range_<int>(0, rows), col_range);
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::row(int row) {
		assert(0 <= row && row < rows);
		return subMat(Range_<int>(row, row + 1), Range_<int>(0, cols));
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::col(int col) {
		assert(0 <= col && col < cols);
		return subMat(Range_<int>(0, rows), Range_<int>(col, col + 1));
	}

	template<typename _Tp>
	inline void Matrix_<_Tp>::copyTo(Matrix_ & m, Rect4i roi)	{
		assert(roi.y >= 0 && (roi.y + roi.height) <= m.rows);
		assert(roi.x >= 0 && (roi.x + roi.width) <= m.cols);
		assert(roi.height == rows && roi.width == cols);

		size_t cpySize = roi.width*channels;
		for (size_t i = roi.y; i < (roi.y+roi.height); i++) {
			_Tp* ptrD = m.data + i*m.step() + roi.x*m.channels;
			_Tp* ptrS = data + (i-roi.y)*step();
			memcpy(ptrD, ptrS, cpySize * sizeof(_Tp));
		}
	}

	template<typename _Tp>
	template<typename _Tp2>
	inline void Matrix_<_Tp>::convertTo(Matrix_<_Tp2>& dst) {
		if (!dst.empty()) dst.release();
		int i = 0;
		dst.create(rows, cols, channels, 0);
		if (this->totalSizes() > 8) {
			for (; i < this->totalSizes() - 8; i += 8) {
				*(dst.data + i + 0) = (_Tp2)*(this->data + i + 0);
				*(dst.data + i + 1) = (_Tp2)*(this->data + i + 1);
				*(dst.data + i + 2) = (_Tp2)*(this->data + i + 2);
				*(dst.data + i + 3) = (_Tp2)*(this->data + i + 3);
				*(dst.data + i + 4) = (_Tp2)*(this->data + i + 4);
				*(dst.data + i + 5) = (_Tp2)*(this->data + i + 5); 
				*(dst.data + i + 6) = (_Tp2)*(this->data + i + 6);
				*(dst.data + i + 7) = (_Tp2)*(this->data + i + 7);
			}
		}
		for (; i < this->totalSizes(); i++)	{
			dst.at(i) = (_Tp2)this->at(i);
		}
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::t() {
		assert(channels == 1);
		Matrix_<_Tp> t;
		t.create(cols, rows, 1, 0);
		for (size_t i = 0; i < t.rows; i++) {
			for (size_t j = 0; j < t.cols; j++) {
				t.at(i, j) = this->at(j, i);
			}
		}
		return t;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::concatenate(const Matrix_ & m, int axes) {
		assert(axes >= 0 && axes < 2);
		Matrix_<_Tp> dst;
		switch (axes)
		{
			//concatenate at axes 0 =>row
			case 0:
			{
				assert(cols == m.cols && channels == m.channels);
				int w = cols;
				int h = rows + m.rows;
				int c = channels;
				dst.create(h, w, c);
				Rect4i roi0(0, 0, cols, rows);
				Rect4i roi1(0, rows, cols, m.rows);
				this->copyTo(dst, roi0);
				Matrix_(m).copyTo(dst, roi1);

			}break;
			//concatenate at axes 1 =>col
			case 1:
			{
				assert(rows == m.rows && channels == m.channels);
				int w = cols + m.cols;
				int h = rows;
				int c = channels;
				dst.create(h, w, c,0);
				Rect4i roi0(0, 0, cols, rows);
				this->copyTo(dst, roi0);
				size_t cpySize = m.step();
				for (size_t i = 0; i < rows; i++) {
					_Tp* ptrS = m.data + i*m.step();
					_Tp* ptrD = dst.data + i*dst.step() + cols*channels;
					memcpy(ptrD, ptrS, cpySize * sizeof(_Tp));
				}
			}break;
			default:
				break;
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::zeros(int r, int c, int cn) {
		Matrix_<_Tp> dst;
		dst.create(r, c, cn, 0);
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::ones(int r, int c, int cn) {
		Matrix_<_Tp> dst;
		dst.create(r, c, cn, 1);
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::eyes(int r, int c) {
		assert(r == c);
		Matrix_<_Tp> dst;
		dst.create(r, c, 1,0);
		for (size_t i = 0; i < r; i++) {
			dst.at(i, i) = _Tp(1.);
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::shuffle(const Matrix_ & m) {
		Matrix_<_Tp> dst(m.rows,m.cols,m.channels);
		std::vector<_Tp*> pVec;
		for (size_t i = 0; i < m.rows; i++) {
			_Tp* p = m.ptr(i);
			pVec.push_back(p);
		}

		//! shuffle data address of each rows
		auto randFunc = [&](int i)->int {
			return int(_rand<double>(0, 0xffff)) % i;
		};
		std::random_shuffle(pVec.begin(),
			pVec.end(), randFunc);

		size_t cpySizeOfPerRow = m.step();
		for (size_t i = 0; i < dst.rows; i++) {
			memcpy(dst.ptr(i), pVec[i], cpySizeOfPerRow * sizeof(_Tp));
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::rand(int r, int c, int cn) {
		Matrix_<_Tp> dst;
		dst.create(r, c, cn, 0);
		for (size_t i = 0; i < dst.totalSizes(); i++) {
			*(dst.data + i) = _rand<double>(0, 1);//std::rand() / (RAND_MAX + 1.0);
		}
		return dst;
	}

	template<typename _Tp>
	inline Matrix_<_Tp> Matrix_<_Tp>::randn(int r, int c, int cn, float esp) {
		Matrix_<_Tp> dst;
		dst.create(r, c, cn, 0);
		for (size_t i = 0; i < dst.totalSizes(); i++) {
			*(dst.data + i) =(_Tp) _randn<double>((double)-esp, (double)esp);
		}
		return dst;
	}

};
