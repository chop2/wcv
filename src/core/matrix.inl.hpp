#pragma once
#include "matrix.h"

namespace wcv {
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
		rows = h;
		cols = w;
		channels = c;
		this->data = (_Tp*)fastAlloc(w*h*c * sizeof(_Tp));
		memset(this->data, val, w*h*c * sizeof(_Tp));
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
	inline _Tp * Matrix_<_Tp>::ptr(int i0) {
		return (_Tp*)(this->data + i0 * step());
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

	template<typename _Tp>
	std::ostream& operator<<(std::ostream& os, const Matrix_<_Tp>& rhs) {
		os.setf(ios::fixed);
		os.precision(8);
		os << "[";
		for (size_t i = 0; i < rhs.rows; i++) {
			for (size_t j = 0; j < rhs.cols; j++) {
				for (size_t k = 0; k < rhs.channels; k++) {
					if (j < rhs.cols - 1)
						os << MAT_ELEM_M(rhs, i, j, k) << ",";
					else
						os << MAT_ELEM_M(rhs, i, j, k) << ";";
				}
			}
			if (i < rhs.rows - 1)
				os << endl;
			else
				os << "]";
		}
		return os;
	};

	template<typename _Tp>
	Matrix_<_Tp> Matrix_<_Tp>::subMat(const Range4i& rowRange,
		const Range4i& colRange) {
		assert(rowRange.end >= rowRange.start &&
			rowRange.end <= rows && rowRange.start >= 0);
		assert(colRange.end >= colRange.start &&
			colRange.end <= cols && colRange.start >= 0);
		Matrix_ submat(rowRange.end - rowRange.start,
			colRange.end - colRange.start, channels);
		size_t cpySize = (colRange.end - colRange.start)*channels;
		for (size_t i = rowRange.start; i < rowRange.end; i++) {
			uchar* ptrS = this->data + i*step() + colRange.start*channels;
			uchar* ptrD = submat.data + i*submat.step();
			memcpy(ptrD, ptrS, cpySize * sizeof(_Tp));
		}
		return submat;
	};

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

};