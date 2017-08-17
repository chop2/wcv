#include <iostream>
#include <string>
#include "src/imgproc/imgproc.h"
#include "src/utils/utility.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "src/core/core.h"

using namespace std;
using namespace wcv;
using namespace cv;

void test_matrix() {
	float dat[] = {
		-1.5,0.5,1.5,
		-2.,0.,2.6,
		-1.,0.3,1.7 };
	float dat1[] = {
		1,2,3,4,5,6,7,8,9,
		2,3,4,5,6,7,8,9,10,
		3,4,5,6,7,8,9,10,11,
		4,5,6,7,8,9,10,11,12,
		5,6,7,8,9,10,11,12,13,
		6,7,8,9,10,11,12,13,14,
		7,8,9,10,11,12,13,14,15,
		8,9,10,11,12,13,14,15,16,
		9,10,11,12,13,14,15,16,17
	};
	//Mat32f k = Mat32f(3, 3, 3, (uchar*)(&dat1)), kk;

	// test<at>
	Mat32f k = Mat32f(3, 3, 1, &dat);
	float val = k.at(0);
	cout << "test \"at\": "<<val << endl;
	val = k.at(1, 2);
	cout << "test \"at\": " << val << endl;

	// test<col>
	Mat32f kk = k.col(2);
	cout << "test \"col\": " << endl << kk.toString(false) << endl;

	// test<row>
	kk = k.row(2);
	cout << "test \"row\": " << endl << kk.toString(false) << endl;

	// test<colRange>
	kk = k.colRange(Range4i(1, 3));
	cout << "test \"colRange\": " << endl << kk.toString(false) << endl;

	// test<rowRange>
	kk = k.rowRange(Range4i(0, 2));
	cout << "test \"rowRange\": " << endl << kk.toString(false) << endl;

	// test<concatenate>
	kk = k.concatenate(k, 1);
	cout << "test \"rowRange\": " << endl << kk.toString(false) << endl;

	// test<copyTo>
	kk.create(9, 9, 1, 0);
	k.copyTo(kk, Rect4i(3, 3, 3, 3));
	cout << "test \"copyTo\": " << endl << kk.toString(false) << endl;

	// test<eyes>
	kk = Matrix_<float>::eyes(3, 3);
	cout << "test \"eyes\": " << endl << kk.toString(false) << endl;

	// test<dot>
	Scalar4f s = k.dot(k);
	cout << "test \"eyes\": " << endl << s << endl;

	// test<hadamardProduct>
	kk = k.hadamardProduct(k);
	cout << "test \"hadamardProduct\": " << endl << kk.toString(false) << endl;

	// test<*>
	Scalar4d ss = Scalar4d::all(2);
	kk = k * ss;
	cout << "test \"*\": " << endl << kk.toString(false) << endl;

	// test<*>
	kk = k * k;
	cout << "test \"*\": " << endl << kk.toString(false) << endl;

	// test</>
	kk = k / 2;
	cout << "test \"/\": " << endl << kk.toString(false) << endl;

	// test<+>
	kk = k + 2;
	cout << "test \"+\": " << endl << kk.toString(false) << endl;

	// test<+>
	kk = k + k;
	cout << "test \"+\": " << endl << kk.toString(false) << endl;

	// test<->
	kk = k - 2;
	cout << "test \"-\": " << endl << kk.toString(false) << endl;

	// test<->
	kk = k - k;
	cout << "test \"-\": " << endl << kk.toString(false) << endl;

	// test<==>
	bool b = k == kk;
	cout << "test \"==\": " << endl << b << endl;
	
	// test<!=>
	b = k != k;
	cout << "test \"!=\": " << endl << b << endl;

	// test<ones()>
	kk = Mat32f::ones(3, 3, 1);
	cout << "test \"ones()\": " << endl << kk.toString(false) << endl;

	// test<zeros()>
	kk = Mat32f::zeros(3, 3, 1);
	cout << "test \"zeros()\": " << endl << kk.toString(false) << endl;

	// test<rand()>
	kk = Mat32f::rand(3, 3, 1);
	cout << "test \"rand()\": " << endl << kk.toString(false) << endl;

	// test<randn()>
	kk = Mat32f::randn(3, 3, 1);
	cout << "test \"randn()\": " << endl << kk.toString(false) << endl;

	// test<shuffle()>
	Mat32f ks = Mat32f(9, 9, 1, &dat1);
	kk = Mat32f::shuffle(ks);
	cout << "test \"shuffle()\": " << endl << kk.toString(false) << endl;

	// test<t()>
	kk = k.t();
	cout << "test \"t()\": " << endl << kk.toString(false) << endl;

	// test<converTo>
	Mat32s k32s;
	k.convertTo(k32s);
	cout << "test \"convertTo()\": " << endl << k32s.toString(false) << endl;
}

void test_matop() {
	Mat64f kernel = getDefaultGaussianKernel2D_3x3();
	//cout << kernel << endl;

	kernel = getGaussianKernel1D(7, 1);
	//cout << kernel << endl;

	kernel = getGaussianKernel2D(Size4i(7, 5), 1.);
	cout << kernel << endl;

	Scalar4d s = sum(kernel);
	cout << s << endl;
}

void test_imgproc() {
	Mat src = imread("3.jpg", 0);

	Image img,ds;
	img.from_cvmat(src);
	wcv::threshold(img, img, 50);
	wcv::not(img, img);
	cv::Mat dd = img.to_cvmat();
	//wcv::cvtColorGray(img, ds);
	//wcv::equalize(img, ds);
	//wcv::translation(img, ds, -10, -10,Scalar4i::all(255));
	//wcv::mirror(img, ds, false);
	//wcv::tranpose(img, ds);
	//wcv::resize(img, ds, Size4i(src.cols * 4, src.rows * 4));
	//wcv::rotate(img, ds, 30);
	//wcv::copymakeBoarder(img, Size4i(7, 7), wcv::MIRROR, ds);
	//wcv::gaussianblur(img, ds, Size4i(5, 5));
	//wcv::medianBlur(img, ds, Size4i(3, 3));
	//wcv::boxBlur(img, ds, Size4i(5, 5));
	//wcv::graySharp(img, ds, 128);

	//vector<Mat8u> mvs;
	//wcv::split(img, mvs);
	//Mat8u img;
	//wcv::merge(mvs, imgm);
	Image elem = getStructElement(Size4i(3, 3), MORP_SHAPE_CROSS);
	cout << elem << endl;
	wcv::erode(img, ds, elem);
	//wcv::dilate(img, ds, elem);

	Mat dst = ds.to_cvmat();
	imwrite("te.jpg", dst);
}

void main()
{
	test_matop();
	test_matrix();
	test_imgproc();
}

