#include "imgproc.h"
#include "../utils/utility.h"
#include <iostream>

#define CHECK_IN_ROI(x,y,size)			\
		((x >= 0 && x < size.width) &&	\
		(y >= 0 && y < size.height))

/**@brief get pixel element from single image*/
#define _PIXEL_ELEM_S(img,r,c)			\
		(*(img.data + r*img.step() + c))

#define _PIXEL_ELEM_M(img,r,c,k)		\
		(*(img.data + r*img.step() + c*img.channels + k))


namespace wcv {
	void threshold(const Image & src, Image & dst, int t,int max_val) {
		assert(src.checkValid() && src.channels == 1);
		if (dst.empty()) {
			dst.create(src.rows, src.cols, src.channels);
		}

		int i = 0;
		size_t size = src.cols*src.rows*src.channels;
		int step = src.cols * src.channels;
#ifdef USE_SSE
		__m128i s,tt,v;
		tt = _mm_set1_epi8 (short(t));
		v = _mm_set1_epi8(short(max_val));
		for (; i < size - 16; i += 16) {
			s = _mm_load_si128((__m128i*)(src.data + i));
			__m128i msk = _mm_cmpgt_epi8(s, tt);
			msk = _mm_and_si128(msk, v);
			_mm_store_si128((__m128i*)(dst.data + i), msk);
		}
#endif
		for (;  i< size; i++) {
			uchar& vs = *(src.data + i);
			uchar& vd = *(dst.data + i);
			vd = vs > t ? max_val : 0;
		}
	}

	void windowTrans(const Image & src, Image & dst, int l, int u) {
		assert(src.checkValid() && src.channels == 1);
		if (dst.empty())
			dst.create(src.rows, src.cols, src.channels);
		int step = src.cols * src.channels;

		for (size_t i = 0; i < src.rows; i++)	{
			for (size_t j = 0; j < src.cols; j++) {
				uchar& pData = *(src.data + i * step + j);
				uchar& pDst = *(dst.data + i * step + j);
				pDst = pData < l ? 0 : pDst;
				pDst = pData > u ? 255 : pDst;
			}
		}
	}

	void cvtColorGray(const Image & src, Image & dst) {
		assert(src.checkValid() && src.channels == 3);		
		if (dst.empty()) {
			dst.create(src.rows, src.cols, 1);
		}
		int i = 0,k = 0;
		float bRat = 0.114,gRat = 0.581,rRat = 0.299;
		for (; i < src.totalSizes(); i+=3,k++) {
			uchar& b = *(src.data + i + 0);
			uchar& g = *(src.data + i + 1);
			uchar& r = *(src.data + i + 2);
			*(dst.data + k) = b * bRat + g*gRat + r*rRat;
		}
	}

	void equalize(const Image & src, Image & dst) {
		assert(src.checkValid() && src.channels == 1);
		if (dst.empty()) {
			dst.create(src.rows, src.cols, src.channels);
		}
		int table[256] = { 0 };
		int acchist[256] = { 0 };
		int size_area = src.cols * src.rows;
		for (size_t i = 0; i < src.totalSizes(); i++) {
			table[*(src.data + i)]++;
		}
		acchist[0] = table[0];
		for (size_t i = 1; i < 256; i++) {
			acchist[i] = acchist[i - 1] + table[i];
		}

		for (size_t i = 0; i < src.totalSizes(); i++) {
			uchar& pix = *(src.data + i);
			*(dst.data + i) = (uchar)((float)acchist[pix] / size_area * 255.0);
		}
	}

	void translation(const Image & src, Image & dst, int dx, int dy, Scalar4i fillColor) {
		assert(src.checkValid());
		if (dst.empty()) {
			dst.create(src.rows, src.cols, src.channels);
		}

		for (size_t i = 0; i < dst.rows; i++) {
			for (size_t j = 0; j < dst.cols; j++) {
				int x0 = j - dx;
				int y0 = i - dy;
				for (size_t c = 0; c < src.channels; c++) {
					if ((x0 >= 0 && x0 < src.cols) &&
						y0 >= 0 && y0 < src.rows) {
						uchar& pix = *(src.data + y0 * src.step() + x0*src.channels + c);
						*(dst.data + i * dst.step() + j*dst.channels + c) = pix;
					}
					else {
						for (size_t c = 0; c < src.channels; c++) {
							*(dst.data + i * dst.step() + j*dst.channels + c) = fillColor[c];
						}
					}
				}
			}
		}
	}

	void mirror(const Image & src, Image & dst, bool bhor)	{
		assert(src.checkValid());
		if (dst.empty()) {
			dst.create(src.rows, src.cols, src.channels);
		}

		for (size_t i = 0; i < src.rows; i++)	{
			for (size_t j = 0; j < src.cols; j++) {
				int x0, y0;
				if (bhor) {
					x0 = src.cols - j;
					y0 = i;
				} else {
					x0 = j;
					y0 = src.rows - i;
				}

				for (size_t c = 0; c < src.channels; c++) {
					if ((x0 >= 0 && x0 <= src.cols) &&
						y0 >= 0 && y0 <= src.rows) {
						uchar& pix = *(src.data + y0 * src.step() + x0*src.channels + c);
						*(dst.data + i*dst.step() + j*dst.channels + c) = pix;
					}
				}
			}
		}
	}

	void tranpose(const Image & src, Image & dst) {
		assert(src.checkValid());
		if (dst.empty()) {
			dst.create(src.cols, src.rows, src.channels);
		}
		
		for (size_t i = 0; i < dst.rows; i++)	{
			for (size_t j = 0; j < dst.cols; j++) {
				int x0 = i;
				int y0 = j;
				for (size_t c = 0; c < src.channels; c++) {
					uchar& pix = *(src.data + y0 * src.step() + x0*src.channels + c);
					*(dst.data + i * dst.step() + j*dst.channels + c) = pix;
				}
			}
		}
	}

	void resize(const Image & src, Image & dst, Size4i size) {
		assert(src.checkValid());
		if (dst.empty()) {
			dst.create(size.height, size.width, src.channels);
		}
		int w = src.cols;
		int h = src.rows;
		int c = src.channels;

		float fx = float(size.width) / float(w);
		float fy = float(size.height) / float(h);

		for (size_t i = 0; i < dst.rows; i++) {
			for (size_t j = 0; j < dst.cols; j++) {
				int x0 = j / fx;
				int y0 = i / fy;
				for (size_t c = 0; c < src.channels; c++) {
					if (CHECK_IN_ROI(x0, y0, src.size())) {
						uchar& pix = (uchar)_PIXEL_ELEM_M(src, y0, x0, c);
						_PIXEL_ELEM_M(dst, i, j, c) = pix;
					}
				}
			}
		}
	}

	void rotate(const Image & src, Image & dst, float angle, Scalar4i fillColor) {
		assert(src.checkValid());
		float theta = angle / 180 * PI;
		float sin_theta = sin(theta);
		float cos_theta = cos(theta);
		float cols = src.cols;
		float rows = src.rows;
		float X1, Y1, X2, Y2, X3, Y3, X4, Y4;
		float XX1, YY1, XX2, YY2, XX3, YY3, XX4, YY4;
		X1 = -(cols - 1) / 2, Y1 = (rows - 1) / 2;
		X2 = (cols - 1) / 2, Y2 = (rows - 1) / 2;
		X3 = -(cols - 1) / 2, Y3 = -(rows - 1) / 2;
		X4 = (cols - 1) / 2, Y4 = -(rows - 1) / 2;

		XX1 = cos_theta * X1 + sin_theta * Y1;
		YY1 = -sin_theta* X1 + cos_theta * Y1;
		XX2 = cos_theta * X2 + sin_theta * Y2;
		YY2 = -sin_theta * X2 + cos_theta * Y2;
		XX3 = cos_theta * X3 + sin_theta * Y3;
		YY3 = -sin_theta * X3 + cos_theta * Y3;
		XX4 = cos_theta * X4 + sin_theta * Y4;
		YY4 = -sin_theta * X4 + cos_theta * Y4;

		float newWidth = std::max(abs(XX4 - XX1), abs(XX3 - XX2));
		float newHeight = std::max(abs(YY4 - YY1), abs(YY3 - YY2));
		if (dst.empty()) {
			dst.create(newHeight, newWidth, src.channels);
		}

		float f1 = -(newWidth - 1) * 0.5 * cos_theta - (newHeight - 1) * 0.5 * sin_theta + (cols - 1) * 0.5;
		float f2 = (newWidth - 1) * 0.5 * sin_theta - (newHeight - 1) * 0.5 * cos_theta + (rows - 1) * 0.5;

		for (size_t i = 0; i < dst.rows; i++)	{
			for (size_t j = 0; j < dst.cols; j++) {
				int i0 = -j * sin_theta + i * cos_theta + f2 + 0.5;
				int j0 = j * cos_theta + i * sin_theta + f1 + 0.5;
				for (size_t c = 0; c < src.channels; c++) {
					if (CHECK_IN_ROI(j0, i0, src.size())) {
						uchar& pix = _PIXEL_ELEM_M(src, i0, j0, c);
						_PIXEL_ELEM_M(dst, i, j, c) = pix;
					} else {
						_PIXEL_ELEM_M(dst, i, j, c) = fillColor[c];
					}
				}
			}
		}
	}

	void gaussianBlur(const Image & src_, Image & dst, Size4i kSize, double sigma) {
		assert(src_.checkValid());
		Image src = src_.clone();
		if (!dst.empty()) dst.release();
		dst.create(src.rows, src.cols, src.channels);
		auto process_single_channel = [&](const Image& s, 
			Size4i ks, double sgm)->Image {
			Mat64f img64f,dst64f;
			Image dst;
			s.convertTo(img64f);
			Mat64f kernel = getGaussianKernel2D(ks, sgm);
			templateOp(img64f, kernel, dst64f, SAME);
			dst64f.convertTo(dst);
			return dst;
		};

		if (src.channels == 1) {
			dst = process_single_channel(src, kSize, sigma);
		} else {
			vector<Image> mvs,mvs1;
			split(src, mvs);
#			ifdef USE_OMP
#			pragma omp parallel for
#			endif
			for (size_t i = 0; i < mvs.size(); i++)	{
				Image d = process_single_channel(mvs[i], kSize, sigma);
				mvs1.push_back(d);
			}
			merge(mvs1, dst);
		}
	}

	void medianBlur(const Image & src_, Image & dst, Size4i kSize) {
		assert(src_.checkValid() && 
			kSize.width % 2 != 0 &&
			kSize.height % 2 != 0);

		Image src_bd;
		int rows = src_.rows;
		int cols = src_.cols;
		int channels = src_.channels;
		int radius_h = kSize.height / 2;
		int radius_w = kSize.width / 2;
		copymakeBoarder(src_, kSize, SAME, src_bd);
		if (!dst.empty()) dst.release();
		dst.create(rows, cols, channels, 0);
		int hist[256] = { 0 };

		auto localhist = [&](Image& s, int i0, int j0, int c) {
			for (int i = -radius_h; i <= radius_h ; i++)	{
				for (int j = -radius_w; j <= radius_w; j++) {
					uchar& pix = s.at(i0 + i, j0 + j, c);
					hist[pix]++;
				}
			}
		};

		auto getMedianVal = [&](int hst[],int& med,int& mNum ) {
			float sum = 0;
			for (size_t i = 0; i < 256; i++) {
				mNum += hist[i];
				sum += hist[i] * 1. / kSize.area();
				if (sum > 0.5) {
					med = i;
					break;
				}
			}
		};

		int th = kSize.area() >> 1 + 1;

		for (size_t c = 0; c < channels; c++) {

			for (int i = radius_h; i < src_bd.rows - radius_h; i++) {
				int med = 0, mNum = 0;
				memset(hist, 0, 256 * sizeof(int));
				localhist(src_bd, i, radius_w, c);
				getMedianVal(hist, med, mNum);
				dst.at(i - radius_h, 0, c) = (uchar)(med);
				for (int j = radius_w + 1; j < src_bd.cols - radius_w; j++) {
					for (int k = -radius_h; k <= radius_h; k++) {
						uchar& pixL = src_bd.at(i + k, j - radius_w - 1, c);
						uchar& pixR = src_bd.at(i + k, j + radius_w, c);
						hist[pixL]--,hist[pixR]++;
						if (pixL <= med)
							mNum--;
						if (pixR <= med)
							mNum++;
					}

					while (mNum > th) {
						if (med == 0)
							break;
						mNum -= hist[med];
						med--;
					}

					while (mNum < th) {
						if (med == 255)
							break;
						med++;
						mNum += hist[med];
					}

					dst.at(i - radius_h, j - radius_w, c) = (uchar)med;
				}
			}

		}
	}

	void boxBlur(const Image & src_, Image & dst, Size4i kSize) {
		assert(src_.checkValid());
		Mat64f kernel;
		Image src = src_.clone();
		if (!dst.empty()) dst.release();
		dst.create(src.rows, src.cols, src.channels, 0);
		double val = 1. / kSize.area();
		kernel.create(kSize.height, kSize.width, 1, val);

		auto process_single_channel = [](Image& s,Mat64f& ker)->Image {
			Mat64f img64f, dst64f;
			Image dst;
			s.convertTo(img64f);
			templateOp(img64f, ker, dst64f, SAME);
			dst64f.convertTo(dst);
			return dst;
		};

		if (src.channels == 1) {
			dst = process_single_channel(src, kernel);
		} else {
			vector<Image> mvs,mvs1;
			split(src, mvs);
#			ifdef USE_OMP
#			pragma omp parallel for
#			endif
			for (size_t i = 0; i < mvs.size(); i++)	{
				Image d = process_single_channel(mvs[i], kernel);
				mvs1.push_back(d);
			}

			merge(mvs1, dst);
		}
	}

	void graySharp(const Image & src_, Image & dst, int thresh) {
		assert(src_.checkValid() && src_.channels == 1);
		Image src = src_.clone();
		if (!dst.empty()) dst.release();
		for (size_t i = 1; i < src.rows; i++) {
			for (size_t j = 1; j < src.cols; j++) {
				uchar& p0 = src.at(i - 1, j - 1);
				uchar& p1 = src.at(i, j - 1);
				uchar& p2 = src.at(i, j);
				int temp = abs(p0 - p1) + abs(p0 - p2);
				if (temp < 255) {
					if (temp >= (uchar)thresh) {
						p0 = temp;
					}
				} else {
					p0 = 255;
				}
			}
		}
		dst = src;
	}
};