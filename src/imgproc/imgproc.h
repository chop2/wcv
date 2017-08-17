#pragma once
#include "../core/core.h"
namespace wcv {

	/**@brief applay threshold
	$f(x) = 0,if f(x) < t else f(x) = max_val$
	@param src - input Image
	@param dst - output result
	@param t - input threshold
	*/
	void threshold(const Image& src,Image& dst,int t,int max_val = 255);

	/**@brief window transform 
	$f(x) = 0,if x < l
	f(x) = x,if l<= x <= u
	f(x) = 255,if x > u
	$
	@param src - input image
	@param dst - output image
	@param l - low value
	@param u - up value
	*/
	void windowTrans(const Image& src, Image& dst, int l, int u);

	/**@brief convert bgr to gray*/
	void cvtColorGray(const Image& src, Image& dst);

	/**@brief image equalize in gray level*/
	void equalize(const Image& src, Image& dst);

	//**********************geometry transform*************************************//
	/**@brief translation transform
	$x0 = x1 - tx,y0 = y1 - ty$
	*/
	void translation(const Image& src, Image& dst,int dx,int dy,Scalar4i fillColor = Scalar4i::all(0));

	/**@brief mirror transform 
	@param src - input image
	@param dst - output image
	@param bhor - horizontal or vertical
	*/
	void mirror(const Image& src, Image& dst, bool bhor);

	/**@brief image transpose transform*/
	void tranpose(const Image& src, Image& dst);

	/**@brief image resize transfrom
	@param src - input image
	@param dst - [out] output image
	@param size - [in] size of new image
	*/
	void resize(const Image& src, Image& dst, Size4i size);

	/**@brief image rotation
	formula:
	$x_{0} = x_{1}cos\theta + y1sin\theta + f1
	y_{0} = -x_{1}sin\theta + y1cos\theta + f2
	f1 = - \frac{newWidth - 1}{2}cos\theta - \frac{newHeight - 1}{2}sin\theta + \frac{cols - 1}{2}
	f2 = \frac{newWidth - 1}{2}sin\theta - \frac{newHeight - 1}{2}cos\theta + \frac{rows- 1}{2}
	$
	@param src - input image
	@param dst - output image
	@param angle - rotate angle,degree	
	@param fillColor - boarder color
	*/
	void rotate(const Image& src, Image& dst, float angle,Scalar4i fillColor=Scalar4i::all(0));

	/**@brief gaussian blur
	@param src - input image
	@param dst - output image
	@param kSize - blur size
	@param sigma - gaussian parameter
	*/
	void gaussianBlur(const Image& src, Image& dst, Size4i kSize, double sigma = 1.);

	/**@brief median blur
	@param src - input image
	@param dst - output image
	@param kSize - apply windows size
	*/
	void medianBlur(const Image& src, Image& dst, Size4i kSize);

	/**@brief box blur*/
	void boxBlur(const Image& src, Image& dst, Size4i kSize);

	/**@brief image sharp
	@param src - input gray image
	@param dst - output image
	@param thresh - input threshold,if gradien(f(x,y)) >= thresh,
			then f(x,y) = gradient(f(x,y)),otherwise f(x,y) will not changed.
	*/
	void graySharp(const Image& src, Image& dst,int thresh);

	/***************************Morphological operator********************************************/

	enum eMorphShape {
		MORP_SHAPE_RECT = 0,
		MORP_SHAPE_ELLIPSE,
		MORP_SHAPE_CROSS
	};

	/**@brief construct a element for morphologic
	@param size - size of element
	@param shape -  shape of element,see eMorphShape
	*/
	Image getStructElement(const Size4i& size,eMorphShape shape );

	/**@brief erode op
	@param src -  input image,binary image
	@param dst - output image
	@param elem - [in] structure element
	@param anchor - anchor of element
	@param iterator - iterator of erode operation
	*/
	void erode(const Image& src, Image& dst, const Image& elem, Point4i anchor = Point4i(-1, -1),int iterator = 1);

	/**@brief dilate op
	@param src -  input image,binary image
	@param dst - output image
	@param elem - [in] structure element
	@param anchor - anchor of element
	@param iterator - iterator of dilate operation
	*/
	void dilate(const Image& src, Image& dst, const Image& elem, Point4i anchor = Point4i(-1, -1), int iterator = 1);

};