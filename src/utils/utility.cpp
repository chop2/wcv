#include "utility.h"


namespace wcv {

//	template<typename _Tp>
//	Matrix_<_Tp>::Matrix_(int h, int w, int c) {
//		create(h, w, c);
//	}
//
//	template<typename _Tp>
//	bool Matrix_<_Tp>::empty() const {
//		return (height*width*nchannels == 0) || (data == nullptr);
//	}
//	
//	template<typename _Tp>
//	void Matrix_<_Tp>::create(int h, int w, int c,int val) {
//		height = h;
//		width = w;
//		nchannels = c;
//		this->data = (uchar*)fastAlloc(w*h*c * sizeof(uchar));
//		memset(this->data, val, w*h*c * sizeof(uchar));
//	}
//	template<typename _Tp>
//	size_t Matrix_<_Tp>::totalSizes() const {
//		return size_t(height * width * nchannels);
//	}
//	
//	template<typename _Tp>
//	int Matrix_<_Tp>::step() const {
//		return width * nchannels;
//	}
//
//	template<typename _Tp>
//	bool Matrix_<_Tp>::checkValid() const {
//		bool bCheck = true;
//		if ((height < 0 || width < 0) ||
//			(height == 0 && width == 0) ||
//			data == nullptr)
//			bCheck = false;
//		return bCheck;
//	}
//
//#ifdef HAVE_OPENCV
//	template<typename _Tp>
//	Matrix_<_Tp>::Matrix_(const cv::Mat& src, bool bCopy) {
//		CV_Assert(!src.empty());
//		int h = src.rows;
//		int w = src.cols;
//		int c = src.channels();
//		this->height = h;
//		this->width = w;
//		this->nchannels = c;
//		if (bCopy) {
//			if (this->data)
//				delete[] this->data;
//			this->data = (uchar*)fastAlloc(w*h*c*sizeof(uchar));
//			memset(this->data, 0, w*h*c*sizeof(uchar));
//		} else {
//			this->data = (uchar*)src.ptr<uchar>();
//		}
//	}
//
//	template<typename _Tp>
//	cv::Mat Matrix_<_Tp>::to_cvmat() {
//		cv::Mat dst; 
//		if (nchannels == 1) {
//			dst = cv::Mat(height, width, CV_8UC1);
//		} else if (nchannels == 3) {
//			dst = cv::Mat(height, width, CV_8UC3);
//		} else {
//			fprintf(stderr, "unsupport format image.");
//			return cv::Mat();
//		}
//		memcpy(dst.ptr<uchar>(),this->data , this->totalSizes()*sizeof(uchar));
//		return dst;
//	}
//
//	template<typename _Tp>
//	void Matrix_<_Tp>::from_cvmat(const cv::Mat & src)	{
//		CV_Assert(!src.empty());
//		this->height = src.rows;
//		this->width = src.cols;
//		this->nchannels = src.channels();
//		if (data)
//			delete[] data;
//		size_t size = height*width*nchannels;
//		data = (uchar*)fastAlloc(size * sizeof(uchar));
//		memcpy(data, src.ptr<uchar>(), size * sizeof(uchar));
//	}
//#endif

	void * fastAlloc(size_t size) {
		uchar* udata = (uchar*)malloc(size + sizeof(void*) + MALLOC_ALIGN);
		if (!udata) {
			fprintf(stderr, "Failed to allocate %lu bytes", (unsigned long)size);
			return nullptr;
		}
		uchar** adata = alignPtr((uchar**)udata + 1, MALLOC_ALIGN);
		adata[-1] = udata;
		return adata;
	}

	void fastFree(void * _mm) {
		if (_mm) {
			uchar* udata = ((uchar**)_mm)[-1];
			free(udata);
		}
	}

	void * aligned_malloc(size_t size, int align)	{
		//align is a power of 2
		assert((align & (align - 1)) == 0);
		void* _mm = malloc(sizeof(void*) + size + align);
		if (!_mm) {
			fprintf(stderr, "out of memory error.");
			return nullptr;
		}

		void** _almm = (void**)_mm + 1;
		void** _almmd = (void**)(((size_t)_almm + align - 1) & -align);
		_almmd[-1] = _mm;
		return _almmd;
	}

	void aligned_free(void* _mm) {
		if (_mm) {
			free(((void**)_mm)[-1]);
		}
	}
};