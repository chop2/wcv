#include <iostream>
#include <string>
#include "src/imgproc/imgproc.h"
#include "src/utils/utility.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "src/core/matOp.hpp"

using namespace std;
using namespace wcv;
using namespace cv;

void main()
{
	Mat src = imread("2.jpg",0);
	//src.convertTo(src, CV_16S);

	Image img;
	img.from_cvmat(src);

	Image ds;
	//wcv::threshold(img, ds, 50);
	//wcv::cvtColorGray(img, ds);
	//wcv::equalize(img, ds);
	//wcv::translationTrans(img, ds, -10, -10,Scalar4i::all(255));
	//wcv::mirrorTrans(img, ds, false);
	//wcv::tranpose(img, ds);
	//wcv::resize(img, ds, Size4i(src.cols * 4, src.rows * 4));
	//wcv::rotate(img, ds, 30);
	//wcv::copymakeBoarder(img, Size4i(7, 7), wcv::MIRROR, ds);

	float dat[] = {
		-1.,0.,1.,
		-2.,0.,2.,
		-1.,0.,1. };
	double dat1[] = {
		1,2,3,4,5,6,7,
		2,3,4,5,6,7,8,
		3,4,5,6,7,8,9,
		4,5,6,7,8,9,10,
		5,6,7,8,9,10,11,
		6,7,8,9,10,11,12,
		7,8,9,10,11,12,13
	};
	Mat64f k = Mat64f(7, 7, 1, (uchar*)(&dat1)),kk;
	wcv::copymakeBoarder(k, Size4i(3, 3), wcv::MIRROR, kk);
	cout << k.toString(false) << endl << endl;;
	cout << kk.toString(false);

	Mat32f imgflt;
	cvt_img_fmt(img, imgflt);
	//wcv::templateOp(img, k, ds, SAME);

	cvt_img_fmt(imgflt, img);

	//cv::Mat dst;
	//thresh(src, dst,20);

	Mat dst = ds.to_cvmat();
	imwrite("te.jpg", dst);
}

