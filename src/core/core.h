#pragma once
#include "core_base.hpp"

namespace wcv {
	typedef _Scalar<int> Scalar4i;
	typedef _Scalar<float> Scalar4f;
	typedef _Scalar<double> Scalar4d;
	
	typedef _Size<int> Size4i;
	typedef _Size<float> Size4f;
	typedef _Size<double> Size4d;

	typedef Matrix_<uchar> Image8u;
	typedef Image8u Image;
	typedef Matrix_<short> Image16s;
	typedef Matrix_<unsigned short> Image16u;
	typedef Matrix_<int> Image32s;
	typedef Matrix_<unsigned int> Image32u;
	typedef Matrix_<float> Image32f;
	typedef Matrix_<double> Image64f;

#	define Reg_cvt_image_fmt(tp1,tp2)								\
	Image##tp2## cvt##tp1##_##tp2##(const Image##tp1##& src) {		\
		Image##tp2## dst;											\
		dst.create(src.height, src.width, src.nchannels);			\
		for (size_t i = 0; i < src.totalSizes; i++)	{				\
			if(tp2 == "8u")											\
				dst.data[i] = (uchar)src.data[i];					\
			else if(tp2 == "16s")									\
				dst.data[i] = (short)src.data[i];					\
			else if(tp2 == "16u")									\
				dst.data[i] = (unsigned short)src.data[i];			\
			else if(tp2 == "32s")									\
				dst.data[i] = (int)src.data[i];						\
			else if(tp2 == "32u")									\
				dst.data[i] = (unsigned int)src.data[i];			\
			else if(tp2 == "32f")									\
				dst.data[i] = (float)src.data[i];					\
			else if(tp2 == "64f")									\
				dst.data[i] = (double)src.data[i];					\
		}															\
	}

	//! 8u <--> 16s
//#	define cvt_img_fmt<uchar, short> cvt8u_16s
//#	define cvt_img_fmt<short, uchar> cvt16s_8u
//	//! 8u <--> 16u
//#	define cvt_img_fmt<uchar,unsigned short> cvt8u_16u
//#	define cvt_img_fmt<unsigned short,uchar> cvt16u_8u
//	//! 8u <--> 32s
//#	define cvt_img_fmt<uchar,int> cvt8u_32s
//#	define cvt_img_fmt<int,uchar> cvt32s_8u
//	//! 8u <--> 32u
//#	define cvt_img_fmt<uchar,unsigned int> cvt8u_32u
//#	define cvt_img_fmt<unsigned int,uchar,> cvt32u_8u
//	//! 8u <--> 32f
//#	define cvt_img_fmt<uchar,float> cvt8u_32f
//#	define cvt_img_fmt<float,uchar> cvt32f_8u
//	//! 8u <--> 64f
//#	define cvt_img_fmt<uchar,double> cvt8u_64f
//#	define cvt_img_fmt<double,uchar> cvt64f_8u
}